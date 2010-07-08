/****************************************************************************
 *		pathtracer.cc: A rather simple MC path integrator
 *		This is part of the yafaray package
 *		Copyright (C) 2006  Mathias Wein (Lynx)
 *		Copyright (C) 2009  Rodrigo Placencia (DarkTide)
 *
 *		This library is free software; you can redistribute it and/or
 *		modify it under the terms of the GNU Lesser General Public
 *		License as published by the Free Software Foundation; either
 *		version 2.1 of the License, or (at your option) any later version.
 *
 *		This library is distributed in the hope that it will be useful,
 *		but WITHOUT ANY WARRANTY; without even the implied warranty of
 *		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *		Lesser General Public License for more details.
 *
 *		You should have received a copy of the GNU Lesser General Public
 *		License along with this library; if not, write to the Free Software
 *		Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
 
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/volume.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>

#include <yafraycore/scr_halton.h>
#include <yafraycore/photon.h>
#include <yafraycore/spectrum.h>

#include <utilities/mcqmc.h>
//#include <integrators/integr_utils.h>
#include <sstream>
#include <yafraycore/triangle.h>
#include <yafraycore/meshtypes.h>
#include <limits>

__BEGIN_YAFRAY

class YAFRAYPLUGIN_EXPORT pathIntegrator_t: public mcIntegrator_t
{
	public:
		pathIntegrator_t(bool transpShad=false, int shadowDepth=4);
		virtual bool preprocess();
		virtual colorA_t integrate(renderState_t &state, diffRay_t &ray/*, sampler_t &sam*/) const;
		static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);
		enum { NONE, PATH, PHOTON, BOTH };
	protected:
		bool traceCaustics; //!< use path tracing for caustics (determined by causticType)
		bool no_recursive;
		float invNPaths;
		int causticType;
};

pathIntegrator_t::pathIntegrator_t(bool transpShad, int shadowDepth)
{
	type = SURFACE;
	trShad = transpShad;
	sDepth = shadowDepth;
	causticType = PATH;
	rDepth = 6;
	maxBounces = 5;
	nPaths = 64;
	invNPaths = 1.f/64.f;
	no_recursive = false;
	integratorName = "PathTracer";
	integratorShortName = "PT";
}

bool pathIntegrator_t::preprocess()
{
	std::stringstream set;
	background = scene->getBackground();
	lights = scene->lights;
	
	if(trShad) set << "ShadowDepth: [" << sDepth << "]"; 

	if(!set.str().empty()) set << "+";
	set << "RayDepth: [" << rDepth << "]";

	bool success = true;
	traceCaustics = false;
	
	if(causticType == PHOTON || causticType == BOTH)
	{
		success = createCausticMap();
	}

	if(causticType == BOTH || causticType == PATH) traceCaustics = true;
	
	if(causticType == PATH)
	{
		if(!set.str().empty()) set << "+";
		set << "Caustics: Path";
	}
	else if(causticType == PHOTON)
	{
		if(!set.str().empty()) set << "+";
		set << "Caustics: Photon(" << nCausPhotons << ")";
	}
	else if(causticType == BOTH)
	{
		if(!set.str().empty()) set << "+";
		set << "Caustics: Path+Photon(" << nCausPhotons << ")";
	}
	
	settings = set.str();

	scene_t *new_scene = new scene_t;

	int max_to_generate = 32000;
	struct ptd {
		point3d_t p;
		bool t;
	};
	std::vector<ptd> points;
	points.resize(max_to_generate);

	struct disk
	{
		point3d_t p;
		int n_idx;
		float r;
		disk(const point3d_t &p, int n_idx, float r) : p(p), n_idx(n_idx), r(r) {}
	};

	std::vector<vector3d_t> normals;
	std::vector<disk> disks;

	int num_tri_points = 20;
	std::vector<point3d_t> tri_points;
	tri_points.resize(num_tri_points);

	scene_t::objDataArray_t &meshes = scene->getMeshes();
	for(scene_t::objDataArray_t::iterator itr = meshes.begin(); itr != meshes.end(); ++itr) {
		triangleObject_t *obj = itr->second.obj;
		int nt = obj->numPrimitives();
		const triangle_t **tris = new const triangle_t*[nt];
		itr->second.obj->getPrimitives(tris);
		for(int i = 0; i < nt; ++i) {
			const triangle_t *t = tris[i];
			point3d_t a, b, c;
			t->getVertices(a, b, c);
			normals.push_back(t->getNormal());
			
			vector3d_t ab = b - a, ac = c - a;
			float area = 0.5 * (ab ^ ac).length();

			// TODO: compute these from the area
			int nr_to_generate = 50;
			int nr_to_keep = 10;
			
			
			for(int j = 0; j < nr_to_generate; j++)
			{
				float u = ourRandom();
				float v = ourRandom();
				points[i].p = a + u * ab + v * ac; 
				points[i].t = false;
			}

			for(int j = 0; j < nr_to_keep; j++)
			{
				int ppoz = 0;
				float mind = std::numeric_limits<float>::max();

				for(int k = 0; k < nr_to_generate; k++)
				{
					if(!points[k].t)
					{
						// TODO: avoid computing the distances more than once
						for(int l = 0; l < j; l++) {
							float d = (points[k].p - disks[l].p).lengthSqr();
							if(d < mind) {
								mind = d;
								ppoz = k;
							}
						}
					}
				}

				point3d_t &p = points[ppoz].p;

				// TODO: determine this from the area
				float r = 10;
				disks.push_back(disk(p, i, r));

				vector3d_t &n = normals[i];
				for(int j = 0; j < num_tri_points; j++) {
					float beta_j = j * M_2PI / num_tri_points;
					float cos_beta_j = cos(beta_j);
					float sin_beta_j = sin(beta_j);

					float saj_p = (n.x * cos_beta_j + n.y * sin_beta_j);
					float nz_sq = n.z * n.z;
					float sinsq_alfa_j = nz_sq / ( saj_p * saj_p + nz_sq);
					float sin_alfa_j = sqrt(sinsq_alfa_j);
					float cos_alfa_j = sqrt(1 - sinsq_alfa_j);

					float r_sin_alfa_j = r * sin_alfa_j;
					
					tri_points[j] = point3d_t(
						p.x + r_sin_alfa_j * cos_beta_j,
						p.y + r_sin_alfa_j * sin_beta_j,
						p.z + r * cos_alfa_j
					);
				}
			}
		}
	}

	return success;
}

colorA_t pathIntegrator_t::integrate(renderState_t &state, diffRay_t &ray/*, sampler_t &sam*/) const
{
	static int calls=0;
	++calls;
	color_t col(0.0);
	CFLOAT alpha=0.0;
	surfacePoint_t sp;
	void *o_udat = state.userdata;
	float W = 0.f;
	//shoot ray into scene
	if(scene->intersect(ray, sp))
	{
		// if camera ray initialize sampling offset:
		if(state.raylevel == 0)
		{
			state.includeLights = true;
			//...
		}
		unsigned char userdata[USER_DATA_SIZE+7];
		userdata[0] = 0;
		state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
		BSDF_t bsdfs;
		const material_t *material = sp.material;
		material->initBSDF(state, sp, bsdfs);
		vector3d_t wo = -ray.dir;
		const volumeHandler_t *vol;
		color_t vcol(0.f);

		// contribution of light emitting surfaces		
		if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);
		
		if(bsdfs & BSDF_DIFFUSE)
		{
			col += estimateAllDirectLight(state, sp, wo);
			if(causticType == PHOTON || causticType == BOTH) col += estimateCausticPhotons(state, sp, wo);
		}
				
		// path tracing:
		// the first path segment is "unrolled" from the loop because for the spot the camera hit
		// we do things slightly differently (e.g. may not sample specular, need not to init BSDF anymore,
		// have more efficient ways to compute samples...)
		
		bool was_chromatic = state.chromatic;
		BSDF_t path_flags = no_recursive ? BSDF_ALL : (BSDF_DIFFUSE);
		
		if(bsdfs & path_flags)
		{
			color_t pathCol(0.0), wl_col;
			path_flags |= (BSDF_DIFFUSE | BSDF_REFLECT | BSDF_TRANSMIT);
			int nSamples = std::max(1, nPaths/state.rayDivision);
			for(int i=0; i<nSamples; ++i)
			{
				void *first_udat = state.userdata;
				unsigned char userdata[USER_DATA_SIZE+7];
				void *n_udat = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
				unsigned int offs = nPaths * state.pixelSample + state.samplingOffs + i; // some redunancy here...
				color_t throughput( 1.0 );
				color_t lcol, scol;
				surfacePoint_t sp1=sp, sp2;
				surfacePoint_t *hit=&sp1, *hit2=&sp2;
				vector3d_t pwo = wo;
				ray_t pRay;

				state.chromatic = was_chromatic;
				if(was_chromatic) state.wavelength = RI_S(offs);
				//this mat already is initialized, just sample (diffuse...non-specular?)
				float s1 = RI_vdC(offs);
				float s2 = scrHalton(2, offs);
				if(state.rayDivision > 1)
				{
					s1 = addMod1(s1, state.dc1);
					s2 = addMod1(s2, state.dc2);
				}
				// do proper sampling now...
				sample_t s(s1, s2, path_flags);
				scol = material->sample(state, sp, pwo, pRay.dir, s, W);
				
				scol *= W;
				throughput = scol;
				state.includeLights = false;

				pRay.tmin = MIN_RAYDIST;
				pRay.tmax = -1.0;
				pRay.from = sp.P;
				
				if(!scene->intersect(pRay, *hit)) continue; //hit background

				state.userdata = n_udat;
				const material_t *p_mat = hit->material;
				BSDF_t matBSDFs;
				p_mat->initBSDF(state, *hit, matBSDFs);
				pwo = -pRay.dir;
				lcol = estimateOneDirectLight(state, *hit, pwo, offs);
				if(matBSDFs & BSDF_EMIT) lcol += p_mat->emit(state, *hit, pwo);

				pathCol += lcol*throughput;
				
				bool caustic = false;
				
				for(int depth = 1; depth < maxBounces; ++depth)
				{
					int d4 = 4*depth;
					s.s1 = scrHalton(d4+3, offs); //ourRandom();//
					s.s2 = scrHalton(d4+4, offs); //ourRandom();//

					if(state.rayDivision > 1)
					{
						s1 = addMod1(s1, state.dc1);
						s2 = addMod1(s2, state.dc2);
					}

					s.flags = BSDF_ALL;
					
					scol = p_mat->sample(state, *hit, pwo, pRay.dir, s, W);
					scol *= W;
					
					if(scol.isBlack()) break;
					
					throughput *= scol;
					caustic = traceCaustics && (s.sampledFlags & (BSDF_SPECULAR | BSDF_GLOSSY | BSDF_FILTER));
					state.includeLights = caustic;

					pRay.tmin = MIN_RAYDIST;
					pRay.tmax = -1.0;
					pRay.from = hit->P;

					if(!scene->intersect(pRay, *hit2)) //hit background
					{
						if((caustic && background))
						{
							pathCol += throughput * (*background)(pRay, state);
						}
						break;
					}
					
					std::swap(hit, hit2);
					p_mat = hit->material;
					p_mat->initBSDF(state, *hit, matBSDFs);
					pwo = -pRay.dir;

					if(matBSDFs & BSDF_DIFFUSE) lcol = estimateOneDirectLight(state, *hit, pwo, offs);
					else lcol = color_t(0.f);

					if((matBSDFs & BSDF_VOLUMETRIC) && (vol=p_mat->getVolumeHandler(hit->N * pwo < 0)))
					{
						if(vol->transmittance(state, pRay, vcol)) throughput *= vcol;
					}
					
					if (matBSDFs & BSDF_EMIT && caustic) lcol += p_mat->emit(state, *hit, pwo);
					
					pathCol += lcol*throughput;
				}
				state.userdata = first_udat;
				
			}
			col += pathCol / nSamples;
		}
		//reset chromatic state:
		state.chromatic = was_chromatic;

		recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);

		CFLOAT m_alpha = material->getAlpha(state, sp, wo);
		alpha = m_alpha + (1.f-m_alpha)*alpha;
	}
	else //nothing hit, return background
	{
		if(background)
		{
			col += (*background)(ray, state, false);
		}
	}

	state.userdata = o_udat;
	return colorA_t(col, alpha);
}

integrator_t* pathIntegrator_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
	bool transpShad=false, noRec=false;
	int shadowDepth = 5;
	int path_samples = 32;
	int bounces = 3;
	int raydepth = 5;
	const std::string *cMethod=0;
	
	params.getParam("raydepth", raydepth);
	params.getParam("transpShad", transpShad);
	params.getParam("shadowDepth", shadowDepth);
	params.getParam("path_samples", path_samples);
	params.getParam("bounces", bounces);
	params.getParam("no_recursive", noRec);
	
	pathIntegrator_t* inte = new pathIntegrator_t(transpShad, shadowDepth);
	if(params.getParam("caustic_type", cMethod))
	{
		bool usePhotons=false;
		if(*cMethod == "photon"){ inte->causticType = PHOTON; usePhotons=true; }
		else if(*cMethod == "both"){ inte->causticType = BOTH; usePhotons=true; }
		else if(*cMethod == "none") inte->causticType = NONE;
		if(usePhotons)
		{
			double cRad = 0.25;
			int cDepth=10, search=100, photons=500000;
			params.getParam("photons", photons);
			params.getParam("caustic_mix", search);
			params.getParam("caustic_depth", cDepth);
			params.getParam("caustic_radius", cRad);
			inte->nCausPhotons = photons;
			inte->nCausSearch = search;
			inte->causDepth = cDepth;
			inte->causRadius = cRad;
		}
	}
	inte->rDepth = raydepth;
	inte->nPaths = path_samples;
	inte->invNPaths = 1.f / (float)path_samples;
	inte->maxBounces = bounces;
	inte->no_recursive = noRec;
	return inte;
}

extern "C"
{

	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("pathtracing",pathIntegrator_t::factory);
	}

}

__END_YAFRAY
