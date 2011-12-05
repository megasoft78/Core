/****************************************************************************
 *		directlight.cc: an integrator for direct lighting only
 *		This is part of the yafaray package
 *		Copyright (C) 2006  Mathias Wein (Lynx)
 *		Copyright (C) 2009  Rodrigo Placencia (DarkTide)
 *
 *		This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <sstream>

#include <yafraycore/irradianceCache.h>

__BEGIN_YAFRAY

class YAFRAYPLUGIN_EXPORT directIC_t : public mcIntegrator_t
{
public:
	directIC_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
	virtual bool preprocess();
	virtual color_t getRadiance(renderState_t &state, ray_t &ray) const;
	virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
	virtual void cleanup();
	static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);
};

directIC_t::directIC_t(bool transpShad, int shadowDepth, int rayDepth)
{
	type = SURFACE;
	causRadius = 0.25;
	causDepth = 10;
	nCausPhotons = 100000;
	nCausSearch = 100;
	trShad = transpShad;
	usePhotonCaustics = false;
	sDepth = shadowDepth;
	rDepth = rayDepth;
	intpb = 0;
	integratorName = "DirectIC";
	integratorShortName = "DIC";
}

bool directIC_t::preprocess()
{
	bool success = true;
	std::stringstream set;
	settings = "";

	if(trShad)
	{
		set << "ShadowDepth: [" << sDepth << "]";
	}
	if(!set.str().empty()) set << "+";
	set << "RayDepth: [" << rDepth << "]";

	background = scene->getBackground();
	lights = scene->lights;

	if(usePhotonCaustics)
	{
		success = createCausticMap();
		if(!set.str().empty()) set << "+";
		set << "Caustics:" << nCausPhotons << " photons. ";
	}

	if(useAmbientOcclusion)
	{
		if(!set.str().empty()) set << "+";
		set << "AO";
	}

	settings = set.str();

	//((qtFilm_t *)imageFilm)->drawTriangle(10, 10, 3);

	// setup cache tree
	if(useIrradianceCache)
	{
		icTree = new icTree_t(scene->getSceneBound(), 10, icMDivs);
	}

	return success;
}

/*! returns the incoming radiance from dir direction
  */
// in dir we return the length of hit
// ray must have from and dir
color_t directIC_t::getRadiance(renderState_t &state, ray_t &ray) const
{
	color_t result;
	ray.tmin = MIN_RAYDIST;
	ray.tmax = -1.0;
	surfacePoint_t hitpoint;
	if (scene->intersect(ray, hitpoint)) {
		BSDF_t matBSDF;
		hitpoint.material->initBSDF(state, hitpoint, matBSDF);
		vector3d_t wo = -ray.dir;
		if (! (matBSDF & BSDF_EMIT) ) {
			if ( matBSDF & (BSDF_DIFFUSE | BSDF_GLOSSY) ) {
				// Totally diffusive! not taking into account any glossyness
				result = estimateAllDirectLight(state, hitpoint, wo);
			}
		}
	} else {
		if (background) {
			result = (*background)(ray, state, false);
		}
		ray.tmax = std::numeric_limits<float>::max();
	}
	return result;
}

colorA_t directIC_t::integrate(renderState_t &state, diffRay_t &ray) const
{
	color_t col(0.0);
	float alpha = 0.0;
	void *o_udat = state.userdata;
	bool oldIncludeLights = state.includeLights;
	surfacePoint_t sp;

	// Shoot ray into scene
	if(scene->intersect(ray, sp)) // If it hits
	{
		// create new memory on stack for material setup
		unsigned char *newUserData[USER_DATA_SIZE];
		state.userdata = (void *)newUserData;
		const material_t *material = sp.material;
		BSDF_t bsdfs;
		material->initBSDF(state, sp, bsdfs);

		vector3d_t wo = -ray.dir;
		if(state.raylevel == 0) state.includeLights = true;

		// obtain material self emitance
		if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);

		if(bsdfs & BSDF_DIFFUSE)
		{
			// obtain direct illumination
			col += estimateAllDirectLight(state, sp, wo);
			col += estimateCausticPhotons(state, sp, wo);
			if(useAmbientOcclusion) col += sampleAmbientOcclusion(state, sp, wo);

			// check for an interpolated result
			if (useIrradianceCache) {
				if (ray.hasDifferentials) {
					icRec_t *icRecord = new icRec_t(icKappa, sp, &(icTree->stratHemi) ); // M, Kappa
					icRecord->setNup(wo);
					icRecord->setPixelArea(ray);
					if (!icTree->getIrradiance(icRecord)) {
						setICRecord(state, ray, icRecord);
						icTree->neighborClamp(icRecord);
						icTree->add(icRecord);
					}
					col += icRecord->irr * icRecord->material->eval(state, *icRecord, wo, icRecord->getNup(), BSDF_DIFFUSE) * M_1_PI;
				} else
					Y_INFO << "NO DIFFERENTIALS!!!" << std::endl;
			}
		}
		// Reflective?, Refractive?
		recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);
		float m_alpha = material->getAlpha(state, sp, wo);
		alpha = m_alpha + (1.f - m_alpha) * alpha;
	}
	else // Nothing hit, return background if any
	{
		if(background) col += (*background)(ray, state, false);
	}
	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;
	return colorA_t(col, alpha);
}

integrator_t* directIC_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
	bool transpShad=false;
	bool caustics=false;
	bool do_AO=false;
	int shadowDepth=5;
	int raydepth=5, cDepth=10;
	int search=100, photons=500000;
	int AO_samples = 32;
	double cRad = 0.25;
	double AO_dist = 1.0;
	color_t AO_col(1.f);
	bool do_IC=false;
	int IC_M=10;
	double IC_K=2.5;
	params.getParam("raydepth", raydepth);
	params.getParam("transpShad", transpShad);
	params.getParam("shadowDepth", shadowDepth);
	params.getParam("caustics", caustics);
	params.getParam("photons", photons);
	params.getParam("caustic_mix", search);
	params.getParam("caustic_depth", cDepth);
	params.getParam("caustic_radius", cRad);
	params.getParam("do_AO", do_AO);
	params.getParam("AO_samples", AO_samples);
	params.getParam("AO_distance", AO_dist);
	params.getParam("AO_color", AO_col);
	params.getParam("do_IC", do_IC);
	params.getParam("IC_M_Divs", IC_M);
	params.getParam("IC_Kappa", IC_K);

	directIC_t *inte = new directIC_t(transpShad, shadowDepth, raydepth);
	// caustic settings
	inte->usePhotonCaustics = caustics;
	inte->nCausPhotons = photons;
	inte->nCausSearch = search;
	inte->causDepth = cDepth;
	inte->causRadius = cRad;
	// AO settings
	inte->useAmbientOcclusion = do_AO;
	inte->aoSamples = AO_samples;
	inte->aoDist = (PFLOAT)AO_dist;
	inte->aoCol = AO_col;
	// IC settings
	inte->useIrradianceCache = do_IC;
	inte->icMDivs = IC_M;
	inte->icKappa = IC_K;
	return inte;
}

void directIC_t::cleanup() {
	icTree->saveToXml("dump.xml");
}

extern "C"
{

	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("directIC",directIC_t::factory);
	}

}

__END_YAFRAY

