/****************************************************************************
 *		mcintegrator.h: A basic abstract integrator for MC sampling
 *		This is part of the yafray package
 *		Copyright (C) 2010  Rodrigo Placencia (DarkTide)
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

#include <core_api/mcintegrator.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <yafraycore/photon.h>
#include <yafraycore/scr_halton.h>
#include <yafraycore/spectrum.h>
#include <utilities/mcqmc.h>
#include <core_api/object3d.h>
#include <core_api/primitive.h>

__BEGIN_YAFRAY

#define allBSDFIntersect (BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
#define loffsDelta 4567 //just some number to have different sequences per light...and it's a prime even...

inline color_t mcIntegrator_t::estimateAllDirectLight(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	color_t col;
	unsigned int loffs = 0;
	for(std::vector<light_t *>::const_iterator l=lights.begin(); l!=lights.end(); ++l)
	{
		col += doLightEstimation(state, (*l), sp, wo, loffs);
		loffs++;
	}
	
	return col;
}

inline color_t mcIntegrator_t::estimateOneDirectLight(renderState_t &state, const surfacePoint_t &sp, vector3d_t wo, int n) const
{
	int lightNum = lights.size();
	
	if(lightNum == 0) return color_t(0.f); //??? if you get this far the lights must be >= 1 but, what the hell... :)
	
	Halton hal2(2);
	hal2.setStart(n-1);
	
	int lnum = std::min((int)(hal2.getNext() * (float)lightNum), lightNum - 1);
	
	return doLightEstimation(state, lights[lnum], sp, wo, lnum) * lightNum;
	//return col * nLights;
}

inline color_t mcIntegrator_t::doLightEstimation(renderState_t &state, light_t *light, const surfacePoint_t &sp, const vector3d_t &wo, const unsigned int  &loffs) const
{
	color_t col(0.f);
	bool shadowed;
	unsigned int l_offs = loffs * loffsDelta;
	const material_t *material = sp.material;
	ray_t lightRay;
	lightRay.from = sp.P;
	color_t lcol(0.f), scol;
	float lightPdf;

	// handle lights with delta distribution, e.g. point and directional lights
	if( light->diracLight() )
	{
		if( light->illuminate(sp, lcol, lightRay) )
		{
			// ...shadowed...
			lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
			shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
			if(!shadowed)
			{
				if(trShad) lcol *= scol;
				color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
				color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
				col += surfCol * lcol * std::fabs(sp.N*lightRay.dir) * transmitCol;
			}
		}
	}
	else // area light and suchlike
	{
		Halton hal2(2);
		Halton hal3(3);
		int n = light->nSamples();
		if(state.rayDivision > 1) n = std::max(1, n/state.rayDivision);
		float invNS = 1.f / (float)n;
		unsigned int offs = n * state.pixelSample + state.samplingOffs + l_offs;
		bool canIntersect=light->canIntersect();
		color_t ccol(0.0);
		lSample_t ls;

		hal2.setStart(offs-1);
		hal3.setStart(offs-1);

		for(int i=0; i<n; ++i)
		{
			// ...get sample val...
			ls.s1 = hal2.getNext();
			ls.s2 = hal3.getNext();
			
			if( light->illumSample (sp, ls, lightRay) )
			{
				// ...shadowed...
				lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
				shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
				
				if(!shadowed && ls.pdf > 1e-6f)
				{
					if(trShad) ls.col *= scol;
					color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
					ls.col *= transmitCol;
					color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
					if( canIntersect)
					{
						float mPdf = material->pdf(state, sp, wo, lightRay.dir, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
						if(mPdf > 1e-6f)
						{
							float l2 = ls.pdf * ls.pdf;
							float m2 = mPdf * mPdf;
							float w = l2 / (l2 + m2);
							ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) * w / ls.pdf;
						}
						else
						{
							ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) / ls.pdf;
						}
					}
					else ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) / ls.pdf;
				}
			}
		}

		col += ccol * invNS;

		if(canIntersect) // sample from BSDF to complete MIS
		{
			color_t ccol2(0.f);

			hal2.setStart(offs-1);
			hal3.setStart(offs-1);

			for(int i=0; i<n; ++i)
			{
				ray_t bRay;
				bRay.tmin = MIN_RAYDIST; bRay.from = sp.P;

				float s1 = hal2.getNext();
				float s2 = hal3.getNext();
				float W = 0.f;

				sample_t s(s1, s2, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
				color_t surfCol = material->sample(state, sp, wo, bRay.dir, s, W);
				if( s.pdf>1e-6f && light->intersect(bRay, bRay.tmax, lcol, lightPdf) )
				{
					shadowed = (trShad) ? scene->isShadowed(state, bRay, sDepth, scol) : scene->isShadowed(state, bRay);
					if(!shadowed && lightPdf > 1e-6f)
					{
						if(trShad) lcol *= scol;
						color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
						lcol *= transmitCol;
						float lPdf = 1.f/lightPdf;
						float l2 = lPdf * lPdf;
						float m2 = s.pdf * s.pdf;
						float w = m2 / (l2 + m2);
						ccol2 += surfCol * lcol * w * W;
					}
				}
			}
			col += ccol2 * invNS;
		}
	}
	return col;
}

bool mcIntegrator_t::createCausticMap()
{
	causticMap.clear();
	ray_t ray;
	std::vector<light_t *> causLights;
	
	for(unsigned int i=0;i<lights.size();++i)
	{
		if(lights[i]->shootsCausticP())
		{
			causLights.push_back(lights[i]);
		}
		
	}

	int numLights = causLights.size();
	progressBar_t *pb;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);

	if(numLights > 0)
	{
		float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
		float fNumLights = (float)numLights;
		float *energies = new float[numLights];
		for(int i=0;i<numLights;++i) energies[i] = causLights[i]->totalEnergy().energy();
		pdf1D_t *lightPowerD = new pdf1D_t(energies, numLights);
		
		Y_INFO << integratorName << ": Light(s) photon color testing for caustics map:" << yendl;
		color_t pcol(0.f);
		
		for(int i=0;i<numLights;++i)
		{
			pcol = causLights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
			lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
			Y_INFO << integratorName << ": Light [" << i+1 << "] Photon col:" << pcol << " | lnpdf: " << lightNumPdf << yendl;
		}
		
		delete[] energies;

		int pbStep;
		Y_INFO << integratorName << ": Building caustics photon map..." << yendl;
		pb->init(128);
		pbStep = std::max(1U, nCausPhotons / 128);
		pb->setTag("Building caustics photon map...");

		bool done=false;
		unsigned int curr=0;
		surfacePoint_t sp1, sp2;
		surfacePoint_t *hit=&sp1, *hit2=&sp2;
		renderState_t state;
		unsigned char userdata[USER_DATA_SIZE+7];
		state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
		while(!done)
		{
			state.chromatic = true;
			state.wavelength = RI_S(curr);
			s1 = RI_vdC(curr);
			s2 = scrHalton(2, curr);
			s3 = scrHalton(3, curr);
			s4 = scrHalton(4, curr);
			
			sL = float(curr) / float(nCausPhotons);
			
			int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
			
			if(lightNum >= numLights)
			{
				Y_ERROR << integratorName << ": lightPDF sample error! " << sL << "/" << lightNum << yendl;
				delete lightPowerD;
				return false;
			}
			
			color_t pcol = causLights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
			ray.tmin = MIN_RAYDIST;
			ray.tmax = -1.0;
			pcol *= fNumLights * lightPdf / lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
			if(pcol.isBlack())
			{
				++curr;
				done = (curr >= nCausPhotons);
				continue;
			}
			BSDF_t bsdfs = BSDF_NONE;
			int nBounces = 0;
			bool causticPhoton = false;
			bool directPhoton = true;
			const material_t *material = 0;
			const volumeHandler_t *vol = 0;
			
			while( scene->intersect(ray, *hit2) )
			{
				if(isnan(pcol.R) || isnan(pcol.G) || isnan(pcol.B))
				{
					Y_WARNING << integratorName << ": NaN (photon color)" << yendl;
					break;
				}
				color_t transm(1.f), vcol;
				// check for volumetric effects
				if(material)
				{
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(hit->Ng * ray.dir < 0)))
					{
						vol->transmittance(state, ray, vcol);
						transm = vcol;
					}

				}
				std::swap(hit, hit2);
				vector3d_t wi = -ray.dir, wo;
				material = hit->material;
				material->initBSDF(state, *hit, bsdfs);
				if(bsdfs & (BSDF_DIFFUSE | BSDF_GLOSSY))
				{
					//deposit caustic photon on surface
					if(causticPhoton)
					{
						photon_t np(wi, hit->P, pcol);
						causticMap.pushPhoton(np);
						causticMap.setNumPaths(curr);
					}
				}
				// need to break in the middle otherwise we scatter the photon and then discard it => redundant
				if(nBounces == causDepth) break;
				// scatter photon
				int d5 = 3*nBounces + 5;
				//int d6 = d5 + 1;

				s5 = scrHalton(d5, curr);
				s6 = scrHalton(d5+1, curr);
				s7 = scrHalton(d5+2, curr);

				pSample_t sample(s5, s6, s7, BSDF_ALL_SPECULAR | BSDF_GLOSSY | BSDF_FILTER | BSDF_DISPERSIVE, pcol, transm);
				bool scattered = material->scatterPhoton(state, *hit, wi, wo, sample);
				if(!scattered) break; //photon was absorped.
				pcol = sample.color;
				// hm...dispersive is not really a scattering qualifier like specular/glossy/diffuse or the special case filter...
				causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
								((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
				// light through transparent materials can be calculated by direct lighting, so still consider them direct!
				directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;
				// caustic-only calculation can be stopped if:
				if(!(causticPhoton || directPhoton)) break;

				if(state.chromatic && (sample.sampledFlags & BSDF_DISPERSIVE))
				{
					state.chromatic=false;
					color_t wl_col;
					wl2rgb(state.wavelength, wl_col);
					pcol *= wl_col;
				}
				ray.from = hit->P;
				ray.dir = wo;
				ray.tmin = MIN_RAYDIST;
				ray.tmax = -1.0;
				++nBounces;
			}
			++curr;
			if(curr % pbStep == 0) pb->update();
			done = (curr >= nCausPhotons);
		}
		pb->done();
		pb->setTag("Caustic photon map built.");
		Y_INFO << integratorName << ": Done." << yendl;
		Y_INFO << integratorName << ": Shot " << curr << " caustic photons from " << numLights <<" light(s)." << yendl;
		Y_INFO << integratorName << ": Stored caustic photons: " << causticMap.nPhotons() << yendl;
		
		delete lightPowerD;
		
		if(causticMap.nPhotons() > 0)
		{
			pb->setTag("Building caustic photons kd-tree...");
			causticMap.updateTree();
			Y_INFO << integratorName << ": Done." << yendl;
		}
		
		if(!intpb) delete pb;

	}
	else
	{
		Y_INFO << integratorName << ": No caustic source lights found, skiping caustic map building..." << yendl;
	}

	return true;
}

inline color_t mcIntegrator_t::estimateCausticPhotons(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	if(!causticMap.ready()) return color_t(0.f);
	
	foundPhoton_t *gathered = new foundPhoton_t[nCausSearch];//(foundPhoton_t *)alloca(nCausSearch * sizeof(foundPhoton_t));
	int nGathered = 0;
	
	float gRadiusSquare = causRadius * causRadius;
	
	nGathered = causticMap.gather(sp.P, gathered, nCausSearch, gRadiusSquare);
	
	gRadiusSquare = 1.f / gRadiusSquare; 
	
	color_t sum(0.f);
	
	if(nGathered > 0)
	{
		const material_t *material = sp.material;
		color_t surfCol(0.f);
		float k = 0.f;
		const photon_t *photon;

		for(int i=0; i<nGathered; ++i)
		{
			photon = gathered[i].photon;
			surfCol = material->eval(state, sp, wo, photon->direction(), BSDF_ALL);
			k = kernel(gathered[i].distSquare, gRadiusSquare);
			sum += surfCol * k * photon->color();
		}
		sum *= 1.f / ( float(causticMap.nPaths()) );
	}
	
	delete [] gathered;
	
	return sum;
}

inline void mcIntegrator_t::recursiveRaytrace(renderState_t &state, diffRay_t &ray, BSDF_t bsdfs, surfacePoint_t &sp, vector3d_t &wo, color_t &col, float &alpha) const
{
	const material_t *material = sp.material;
	spDifferentials_t spDiff(sp, ray);

    state.raylevel++;
    
	if(state.raylevel <= rDepth)
	{
		Halton hal2(2);
		Halton hal3(3);

		// dispersive effects with recursive raytracing:
		if( (bsdfs & BSDF_DISPERSIVE) && state.chromatic)
		{
			state.includeLights = false; //debatable...
			int dsam = 8;
			int oldDivision = state.rayDivision;
			int oldOffset = state.rayOffset;
			float old_dc1 = state.dc1, old_dc2 = state.dc2;
			if(state.rayDivision > 1) dsam = std::max(1, dsam/oldDivision);
			state.rayDivision *= dsam;
			int branch = state.rayDivision*oldOffset;
			float d_1 = 1.f/(float)dsam;
			float ss1 = RI_S(state.pixelSample + state.samplingOffs);
			color_t dcol(0.f), vcol(1.f);
			vector3d_t wi;
			const volumeHandler_t *vol;
			diffRay_t refRay;
			float W = 0.f;
			
			for(int ns=0; ns<dsam; ++ns)
			{
				state.wavelength = (ns + ss1)*d_1;
				state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
				state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
				if(oldDivision > 1)  state.wavelength = addMod1(state.wavelength, old_dc1);
				state.rayOffset = branch;
				++branch;
				sample_t s(0.5f, 0.5f, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_DISPERSIVE);
				color_t mcol = material->sample(state, sp, wo, wi, s, W);
				
				if(s.pdf > 1.0e-6f && (s.sampledFlags & BSDF_DISPERSIVE))
				{
					state.chromatic = false;
					color_t wl_col;
					wl2rgb(state.wavelength, wl_col);
					refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
					dcol += (color_t)integrate(state, refRay) * mcol * wl_col * W;
					state.chromatic = true;
				}
			}
			
			if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
			{
				vol->transmittance(state, refRay, vcol);
				dcol *= vcol;
			}
			
			col += dcol * d_1;
			state.rayDivision = oldDivision;
			state.rayOffset = oldOffset;
			state.dc1 = old_dc1;
			state.dc2 = old_dc2;
		}

		// glossy reflection with recursive raytracing:
		if( bsdfs & BSDF_GLOSSY )
		{
			state.includeLights = false;
			int gsam = 8;
			int oldDivision = state.rayDivision;
			int oldOffset = state.rayOffset;
			float old_dc1 = state.dc1, old_dc2 = state.dc2;
			if(state.rayDivision > 1) gsam = std::max(1, gsam/oldDivision);
			state.rayDivision *= gsam;
			int branch = state.rayDivision*oldOffset;
			unsigned int offs = gsam * state.pixelSample + state.samplingOffs;
			float d_1 = 1.f/(float)gsam;
			color_t gcol(0.f), vcol(1.f);
			vector3d_t wi;
			const volumeHandler_t *vol;
			diffRay_t refRay;

			hal2.setStart(offs);
			hal3.setStart(offs);

			for(int ns=0; ns<gsam; ++ns)
			{
				state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
				state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
				state.rayOffset = branch;
				++offs;
				++branch;

				float s1 = hal2.getNext();
				float s2 = hal3.getNext();
				
				float W = 0.f;

				sample_t s(s1, s2, BSDF_ALL_GLOSSY);
				color_t mcol = material->sample(state, sp, wo, wi, s, W);

				if(s.sampledFlags & BSDF_GLOSSY)
				{
					refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
					if(s.sampledFlags & BSDF_REFLECT) spDiff.reflectedRay(ray, refRay);
					else if(s.sampledFlags & BSDF_TRANSMIT) spDiff.refractedRay(ray, refRay, material->getMatIOR());
					gcol += (color_t)integrate(state, refRay) * mcol * W;
				}

				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					if(vol->transmittance(state, refRay, vcol)) gcol *= vcol;
				}
			}

			col += gcol * d_1;
			state.rayDivision = oldDivision;
			state.rayOffset = oldOffset;
			state.dc1 = old_dc1;
			state.dc2 = old_dc2;
		}

		//...perfect specular reflection/refraction with recursive raytracing...
		if(bsdfs & (BSDF_SPECULAR | BSDF_FILTER) && state.raylevel < 20)
		{
			state.includeLights = true;
			bool reflect=false, refract=false;
			vector3d_t dir[2];
			color_t rcol[2], vcol;
			material->getSpecular(state, sp, wo, reflect, refract, dir, rcol);
			const volumeHandler_t *vol;
			if(reflect)
			{
				diffRay_t refRay(sp.P, dir[0], MIN_RAYDIST);
				spDiff.reflectedRay(ray, refRay);
				color_t integ = integrate(state, refRay);
				
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					if(vol->transmittance(state, refRay, vcol)) integ *= vcol;
				}
				
				col += integ * rcol[0];
			}
			if(refract)
			{
				diffRay_t refRay(sp.P, dir[1], MIN_RAYDIST);
				spDiff.refractedRay(ray, refRay, material->getMatIOR());
				colorA_t integ = integrate(state, refRay);
				
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					if(vol->transmittance(state, refRay, vcol)) integ *= vcol;
				}
				
				col += (color_t)integ * rcol[1];
				alpha = integ.A;
			}
		}

	}
	--state.raylevel;
}

color_t mcIntegrator_t::sampleAmbientOcclusion(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	color_t col(0.f), surfCol(0.f), scol(0.f);
	bool shadowed;
	const material_t *material = sp.material;
	ray_t lightRay;
	lightRay.from = sp.P;
	
	int n = aoSamples;
	if(state.rayDivision > 1) n = std::max(1, n / state.rayDivision);

	unsigned int offs = n * state.pixelSample + state.samplingOffs;
	
	Halton hal2(2);
	Halton hal3(3);

	hal2.setStart(offs-1);
	hal3.setStart(offs-1);

	for(int i = 0; i < n; ++i)
	{
		float s1 = hal2.getNext();
		float s2 = hal3.getNext();
		
		if(state.rayDivision > 1)
		{
			s1 = addMod1(s1, state.dc1);
			s2 = addMod1(s2, state.dc2);
		}
		
		lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is still bad...
		lightRay.tmax = aoDist;
		
		float W = 0.f;
		
		sample_t s(s1, s2, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_REFLECT );
		surfCol = material->sample(state, sp, wo, lightRay.dir, s, W);
		
		if(material->getFlags() & BSDF_EMIT)
		{
			col += material->emit(state, sp, wo) * s.pdf;
		}
		
		shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
		
		if(!shadowed)
		{
			float cos = std::fabs(sp.N * lightRay.dir);
			if(trShad) col += aoCol * scol * surfCol * cos * W;
			else col += aoCol * surfCol * cos * W;
		}
	}
	
	return col / (float)n;
}

void mcIntegrator_t::setICRecord(renderState_t &state, diffRay_t &ray, icRec_t *record) const {
	if (!ray.hasDifferentials)
		Y_INFO << "ERROR: ray from mcIntegrator_t::createNewICRecord() should have differentials" << std::endl;
	float oldRayLength[record->getM()];
	color_t oldRad[record->getM()];
	// we set the projected pixel area on the surface point
	record->setPixelArea(ray);
	ray_t sRay; // ray from hitpoint to hemisphere sample direction
	sRay.from = record->P;
	color_t innerRotValues;
	color_t innerTransValuesU;
	color_t innerTransValuesV;
	color_t radiance;
	unsigned int offs = 8 * state.pixelSample /*+ state.samplingOffs*/;
	record->stratHemi->randomize();
	for (int k=0; k<record->getN(); k++) {
		innerRotValues.black();
		innerTransValuesU.black();
		innerTransValuesV.black();
		for (int j=0; j<record->getM(); j++) {
			// Calculate each incoming radiance of hemisphere at point icRecord
			sRay.dir = record->getSampleHemisphere(j, k, offs);
			radiance = getRadiance(state, sRay);
			// note: oldRad[j] and oldRayLength[j] means L_j,k-1 and r_j,k-1 respectively
			//       oldRad[j-1] and oldRayLength[j-1] means L_j-1,k and r_j-1,k respectively
			if (k>0) {
				if (j>0) {
					// cos2(theta_j-)sin(theta_j-) * (L_j,k - L_j-1,k) / min(r_j,k , r_j-1,k)
					float cosThetaMin = record->stratHemi->getCosThetaMinus(j);
					innerTransValuesU +=
							( (cosThetaMin*cosThetaMin) * (radiance - oldRad[j-1]) )/
							(fmin(sRay.tmax, oldRayLength[j-1])) ;
					// cos(theta_j)[cos(theta_j-) - cos(theta_j+)] * (L_j,k - L_j,k-1) / [sin(theta_j,k) * min(r_j,k , r_j-1,k)]
					innerTransValuesV +=
							( record->stratHemi->getCosTheta(j) * (cosThetaMin - record->stratHemi->getCosThetaPlus(j)) *
								(radiance - oldRad[j]) ) /
							(record->stratHemi->getSinTheta(j) * fmin(sRay.tmax, oldRayLength[j]));
				}
			}
			record->irr += radiance;
			record->changeSampleRadius(sRay.tmax);
			// copy new rays and irradiance values over old ones
			oldRad[j] = radiance;
			oldRayLength[j] = sRay.tmax;
			innerRotValues -= record->stratHemi->getTanTheta(j) * radiance;
		}
		vector3d_t vk(record->stratHemi->getVk(k));
		record->rotGrad[0] += vk * innerRotValues.R;
		record->rotGrad[1] += vk * innerRotValues.G;
		record->rotGrad[2] += vk * innerRotValues.B;
		vector3d_t uk(record->stratHemi->getUk(k));
		vector3d_t vkm(record->stratHemi->getVkMinus(k));
		record->transGrad[0] +=
				( (innerTransValuesU.R * M_2PI / (float)record->getN() ) * uk ) +
				( innerTransValuesV.R * vkm );
		record->transGrad[1] +=
				( (innerTransValuesU.G * M_2PI / (float)record->getN() ) * uk ) +
				( innerTransValuesV.G * vkm );
		record->transGrad[2] +=
				( (innerTransValuesU.B * M_2PI / (float)record->getN() ) * uk ) +
				( innerTransValuesV.B * vkm );
	}
	float k = M_PI / ((float)record->getM() * (float)record->getN());
	record->irr = record->irr * k;
	for (int i=0; i<3; i++) {
		record->rotGrad[i] = changeBasis(
				record->rotGrad[i] * k,
				record->NU,
				record->NV,
				record->getNup());
		record->transGrad[i] = changeBasis(
				record->transGrad[i],
				record->NU,
				record->NV,
				record->getNup());
	}
	// HEURISTICS
	record->clampRbyGradient();
	if (ray.hasDifferentials)
		record->clampRbyScreenSpace();
	record->clampGradient();
}

void mcIntegrator_t::cleanup() {
	// do nothing, if IC implemented, may call the xml dump saving function
}

bool mcIntegrator_t::createSSSMaps()
{
	// init and compute light pdf etc.
	ray_t ray;
	int maxBounces = this->nSSSDepth;
	unsigned int nPhotons=this->nSSSPhotons;
	int numLights = lights.size();
	float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
	float fNumLights = (float)numLights;
	float *energies = new float[numLights];
	for(int i=0;i<numLights;++i) 
		energies[i] = lights[i]->totalEnergy().energy();
	pdf1D_t *lightPowerD = new pdf1D_t(energies, numLights);
	for(int i=0;i<numLights;++i) 
		Y_INFO << "energy: "<< energies[i] <<" (dirac: "<<lights[i]->diracLight()<<")\n";
	delete[] energies;
	
	// init progressbar
	progressBar_t *pb;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);
	
	int pbStep;
	Y_INFO << integratorName << ": Building SSS photon map..." << yendl;
	pb->init(128);
	pbStep = std::max(1U, nPhotons / 128);
	pb->setTag("Building SSS photon map...");
	
	// prepare for shooting photons
	bool done=false;
	unsigned int curr=0;
	surfacePoint_t sp1, sp2;
	surfacePoint_t *hit=&sp1, *hit2=&sp2;
	renderState_t state;
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	//std::cout<<"INFO: caustic map " << cMap.nPhotons() << std::endl;
	
	while(!done)
	{
		// sampling the light to shoot photon
		s1 = RI_vdC(curr);
		s2 = scrHalton(2, curr);
		s3 = scrHalton(3, curr);
		s4 = scrHalton(4, curr);
		
		//sL = RI_S(curr);
		sL = float(curr) / float(nPhotons);
		//sL = float(cMap.nPhotons()) / float(nPhotons);
		int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
		if(lightNum >= numLights){ std::cout << "lightPDF sample error! "<<sL<<"/"<<lightNum<< "  " << curr << "/" << nPhotons << "\n"; delete lightPowerD; return false; }
		
		// shoot photon
		color_t pcol = lights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
		color_t pcol_t;
		ray.tmin = MIN_RAYDIST;
		ray.tmax = -1.0;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
		if(pcol.isBlack())
		{
			++curr;
			//done = (curr >= nPhotons) ? true : false;
			done = (inCount >= nPhotons) ? true : false;
			continue;
		}
		
		// find instersect point
		BSDF_t bsdfs = BSDF_NONE;
		int nBounces=0;
		const material_t *material = 0;
		
		bool isRefrectedOut = false;
		
		while( scene->intersect(ray, *hit2) )
		{
			if(isnan(pcol.R) || isnan(pcol.G) || isnan(pcol.B))
			{ std::cout << "NaN WARNING (photon color)" << std::endl; break; }
			color_t transm(1.f), vcol;
			// check for volumetric effects
			if(material)
			{
				if((bsdfs&BSDF_VOLUMETRIC) && material->volumeTransmittance(state, *hit, ray, vcol))
				{
					transm = vcol;
				}
			}
			std::swap(hit, hit2);
			vector3d_t wi = -ray.dir, wo;
			material = hit->material;
			material->initBSDF(state, *hit, bsdfs);
			if(bsdfs & BSDF_TRANSLUCENT)
			{
				// if photon intersect with SSS material, add this photon to cooresponding object's SSSMap and absorb it
				photon_t np(wi, hit->P, pcol);
				np.hitNormal = hit->N;
				const object3d_t* hitObj = hit->object;
				if(hitObj)
				{
					//std::cout << curr <<" bounces:" << nBounces << std::endl;
					std::map<const object3d_t*, photonMap_t*>::iterator it = SSSMaps.find(hitObj);
					if(it!=SSSMaps.end()){
						// exist SSSMap for this object
						SSSMaps[hitObj]->pushPhoton(np);
						SSSMaps[hitObj]->setNumPaths(curr);
					}
					else {
						// need create a new SSSMap for this object
						//std::cout << "new translucent is " << bsdfs << "   " << hitObj << std::endl;
						photonMap_t* sssMap_t = new photonMap_t();
						sssMap_t->pushPhoton(np);
						sssMap_t->setNumPaths(curr);
						SSSMaps[hitObj] = sssMap_t;
					}
				} 
				break;
				//cMap.pushPhoton(np);
				//cMap.setNumPaths(curr);
			}
			
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			if(nBounces == maxBounces) break;
			// scatter photon
			int d5 = 3*nBounces + 5;
			//int d6 = d5 + 1;
			if(d5+2 <= 50)
			{
				s5 = scrHalton(d5, curr);
				s6 = scrHalton(d5+1, curr);
				s7 = scrHalton(d5+2, curr);
			}
			else
			{
				s5 = ourRandom();
				s6 = ourRandom();
				s7 = ourRandom();
			}
			pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);
			bool scattered = material->scatterPhoton(state, *hit, wi, wo, sample);
			if(!scattered) break; //photon was absorped.
			
			//std::cout << curr << " not translucent objects:" << std::endl;
			
			pcol = sample.color;
			
			ray.from = hit->P;
			ray.dir = wo;
			ray.tmin = 0.001;
			ray.tmax = -1.0;
			++nBounces;
		}
		++curr;
		if(curr % pbStep == 0) pb->update();
		done = (curr >= nPhotons) ? true : false;
		//done = (cMap.nPhotons() >= nPhotons) ? true : false;
	}
	pb->done();
	pb->setTag("SSS photon map built.");
	
	delete lightPowerD;
	
	return true;
}

float mcIntegrator_t::sssScale = 10.f;

bool mcIntegrator_t::createSSSMapsByPhotonTracing()
{
	// init and compute light pdf etc.
	ray_t ray;
	int maxBounces = this->nSSSDepth;
	unsigned int nPhotons=this->nSSSPhotons;
	int numLights = lights.size();
	float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
	float fNumLights = (float)numLights;
	float *energies = new float[numLights];
	for(int i=0;i<numLights;++i) 
		energies[i] = lights[i]->totalEnergy().energy();
	pdf1D_t *lightPowerD = new pdf1D_t(energies, numLights);
	for(int i=0;i<numLights;++i) 
		Y_INFO << "energy: "<< energies[i] <<" (dirac: "<<lights[i]->diracLight()<<")\n";
	delete[] energies;
	
	// init progressbar
	progressBar_t *pb;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);
	
	int pbStep;
	Y_INFO << "SSS: shoot photons " << nPhotons << yendl;
	Y_INFO << integratorName << ": Building SSS photon map by photon tracing..." << yendl;
	pb->init(128);
	pbStep = std::max(1U, nPhotons / 128);
	pb->setTag("Building SSS photon map by photon tracing...");
	
	// prepare for shooting photons
	bool done=false;
	unsigned int curr=0, scatteCount=0, inCount=0, absorbCount=0;
	surfacePoint_t sp1, sp2;
	surfacePoint_t *hit=&sp1, *hit2=&sp2;
	renderState_t state;
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	//std::cout<<"INFO: caustic map " << cMap.nPhotons() << std::endl;
	
	while(!done)
	{
		// sampling the light to shoot photon
		s1 = RI_vdC(curr);
		s2 = scrHalton(2, curr);
		s3 = scrHalton(3, curr);
		s4 = scrHalton(4, curr);
		//sL = RI_S(curr);
		//sL = float(curr) / float(nPhotons);
		sL = float(inCount) / float(nPhotons);
		int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
		if(lightNum >= numLights){ std::cout << "lightPDF sample error! "<<sL<<"/"<<lightNum<< "  " << curr << "/" << nPhotons << "\n"; delete lightPowerD; return false; }
		
		// shoot photon
		color_t pcol = lights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
		ray.tmin = 0.001;
		ray.tmax = -1.0;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
		if(pcol.isBlack())
		{
			++curr;
			done = (curr >= nPhotons) ? true : false;
			
			continue;
		}
		
		// find instersect point
		BSDF_t bsdfs = BSDF_NONE;
		int nBounces=0;
		const material_t *material = 0;
		while( scene->intersect(ray, *hit2) )
		{
			if(isnan(pcol.R) || isnan(pcol.G) || isnan(pcol.B))
			{ std::cout << "NaN WARNING (photon color)" << std::endl; break; }
			color_t transm(1.f), vcol;
			// check for volumetric effects
			if(material)
			{
				if((bsdfs&BSDF_VOLUMETRIC) && material->volumeTransmittance(state, *hit, ray, vcol))
				{
					transm = vcol;
				}
			}
			std::swap(hit, hit2);
			vector3d_t wi = -ray.dir, wo;
			material = hit->material;
			material->initBSDF(state, *hit, bsdfs);
			
			// if the rey intersects with translucent objects.
			if(bsdfs & BSDF_TRANSLUCENT)
			{
				/*std::cout << "enter ray = " << curr << "  wi = " << wi  << std::endl;
				if (isRefrectedOut) {
					std::cout << "In  curr=" << curr << "  wi = " << wi << "  N=" << hit->N << " from=" << ray.from << "   pcol=" << pcol << std::endl;
					isRefrectedOut = false;
				}*/
				color_t diffuseC, sigma_s, sigma_a;
				float IOR;
				TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
				diffuseC = dat->difC;
				sigma_a = dat->sig_a;
				sigma_s = dat->sig_s;
				IOR = dat->IOR;
				float g = 0.f;
				
				float sig_a_ = sigma_a.col2bri();
				float sig_s_ = sigma_s.col2bri();
				float sig_t_ = sig_a_ + (1.f-g)*sig_s_;
				float sig_t_1 = 1.f/sig_t_;
				
				//Halton hal2(7);
				//hal2.setStart(curr);
				
				//std::cout << "random seed " << curr << std::endl;
				
				// if photon intersect with SSS material, get the refract direction and continue to trace this photon.
				if( refract(hit->N, wi, wo, IOR) )
				{
					inCount++;
					if (inCount % pbStep == 0) pb->update();
					
					const object3d_t* refObj = hit->object;
					bool refracOut = false;
					bool isStored = false;
					
					
					// get the refrace try
					float sc1 = ourRandom();//hal2.getNext(); 
					float sc2,sc3;
					float scatteDist = -1.f*log(1-sc1)*sig_t_1/sssScale;
					//float scatteDist = 1.f/(sig_t_1*sssScale);
					vector3d_t sdir = wo;
					ray.from = hit->P;
					ray.dir = wo;
					ray.tmin = MIN_RAYDIST/sssScale;
					ray.tmax = scatteDist;
					
					
					point3d_t scattePt = ray.from + scatteDist*ray.dir;
					float cosWo = ray.dir*(-1.f*hit->N);
					pcol *= diffuseC;
					photon_t np(wi, hit->P, pcol);
					np.hitNormal = hit->N;
					np.sourcePos = scattePt;
					np.sourceDepth = cosWo*scatteDist;
					
					//std::cout << "cosWo = " << cosWo << " scatterDist = " << scatteDist << "  depth = " << np.sourceDepth << std::endl;
					
					if(refObj)
					{
						//std::cout << curr <<" bounces:" << nBounces << std::endl;
						std::map<const object3d_t*, photonMap_t*>::iterator it = SSSMaps.find(refObj);
						if(it!=SSSMaps.end()){
							// exist SSSMap for this object
							SSSMaps[refObj]->pushPhoton(np);
							SSSMaps[refObj]->setNumPaths(curr);
						}
						else {
							// need create a new SSSMap for this object
							//std::cout << "new translucent is " << bsdfs << "   " << hitObj << std::endl;
							photonMap_t* sssMap_t = new photonMap_t();
							sssMap_t->pushPhoton(np);
							sssMap_t->setNumPaths(curr);
							SSSMaps[refObj] = sssMap_t;
						}
					}
					
					//std::cout << "first refracted = " << curr << "  wi = " << wi << "  N=" << hit->N << " wo=" << wo << "   pcol=" << pcol << std::endl;
					
					
					int scatNum = 0;
					
					while (!(refracOut = scene->intersect(ray, *hit2)))
					{
						//std::cout << "ray = " << curr << "  ray.dir = " << ray.dir << "  from=" << ray.from << "  scatteDist=" << ray.tmax << std::endl;
						// compute scatter point
						point3d_t scattePt = ray.from + scatteDist*ray.dir;
						pcol_t = pcol;
						pcol *= fExp(-1*sig_t_*scatteDist*sssScale); // power attuation
						
						if (pcol.energy() < 1e-6) {
							break;
						}
						
						// roulette whether scatter or absorb
						float s = ourRandom();//hal2.getNext();
						if (  s < sig_a_*sig_t_1 ) {
							// absorbed, then break
							//std::cout << "absorbed" << std::endl;
							absorbCount++;
							break;
						}
						else 
						{
							// scattered
							// store photon
							//std::cout << "scattered " << s << " " <<sig_a_*sig_t_1 << std::endl;
							//if (!isStored) {
//								// store photon here
//								//photon_t np(ray.dir, scattePt, pcol);
//								
//								float cosWo = ray.dir*(-1.f*hit->N);
//								photon_t np(wi, hit->P, pcol_t);
//								np.hitNormal = hit->N;
//								np.sourcePos = scattePt;
//								np.sourceDepth = cosWo*scatteDist;
//								
//								//std::cout << "cosWo = " << cosWo << " scatterDist = " << scatteDist << "  depth = " << np.sourceDepth << std::endl;
//								
//								if(refObj)
//								{
//									//std::cout << curr <<" bounces:" << nBounces << std::endl;
//									std::map<const object3d_t*, photonMap_t*>::iterator it = SSSMaps.find(refObj);
//									if(it!=SSSMaps.end()){
//										// exist SSSMap for this object
//										SSSMaps[refObj]->pushPhoton(np);
//										SSSMaps[refObj]->setNumPaths(curr);
//									}
//									else {
//										// need create a new SSSMap for this object
//										//std::cout << "new translucent is " << bsdfs << "   " << hitObj << std::endl;
//										photonMap_t* sssMap_t = new photonMap_t();
//										sssMap_t->pushPhoton(np);
//										sssMap_t->setNumPaths(curr);
//										SSSMaps[refObj] = sssMap_t;
//									}
//								}
//								isStored = true;
//								
//								//break;
//							}
							
							// get the scatter direction
							sc2 = scrHalton(2, scatteCount);
							sc3 = scrHalton(3, scatteCount);
							sdir = SampleSphere(sc2,sc3);
							
							sc1 = ourRandom();//hal2.getNext();
							scatteDist = -1.f*log(1-sc1)*sig_t_1/sssScale;
							ray.from = scattePt;
							ray.dir = sdir;
							ray.tmin = MIN_RAYDIST/sssScale;
							ray.tmax = scatteDist;
							
							scatteCount++;
							scatNum++;
							//if (scatNum >= nBounces) {
							//	break;
							//}
							
						}
					}
					
					if (refracOut) {
						
						//std::cout << "ray = " << curr << "  ray.dir = " << ray.dir << "  from=" << ray.from << "  scatteDist=" << ray.tmax << std::endl;
						
						// refract out
						// compute new direction and
						std::swap(hit, hit2);
						wi = -ray.dir;
						material = hit->material;
						
						if (-1*hit->N*wi<=0) {
							break;
						}
						
						if( refract(-1*hit->N, wi, wo, 1.0f/IOR) )
						{
							//std::cout << "Out curr=" << curr << "  wi = " << wi << "  N=" << hit->N*-1 << " wo=" << wo << "  from=" << ray.from << "   pcol=" << pcol << std::endl;
							
							vector3d_t lastTransmit = hit->P - ray.from;
							scatteDist = lastTransmit.length();
							
							pcol *= fExp(-1*sig_t_*scatteDist*sssScale);
							ray.from = hit->P;
							ray.dir = wo;
							ray.tmin = MIN_RAYDIST/sssScale;
							ray.tmax = -1.0;
							++nBounces;
							//if (!isStored) {
							//std::cout << "photon refacted out" << "  " << pcol << std::endl;
							//}
							isRefrectedOut = true;
							//break;
							continue;
						}
						else
						{
							break;
						}
					}
					else {
						// not refracted out, just because absorbed or power is too small
						break;
					}
				}
				else{
					// not refracted in to object
					break;
				}
			}
			
			if (isDirectLight) {
				break;
			}
			
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			if(nBounces == maxBounces) break;
			// scatter photon
			int d5 = 3*nBounces + 5;
			//int d6 = d5 + 1;
			if(d5+2 <= 50)
			{
				s5 = scrHalton(d5, curr);
				s6 = scrHalton(d5+1, curr);
				s7 = scrHalton(d5+2, curr);
			}
			else
			{
				s5 = ourRandom();
				s6 = ourRandom();
				s7 = ourRandom();
			}
			pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);
			bool scattered = material->scatterPhoton(state, *hit, wi, wo, sample);
			if(!scattered) break; //photon was absorped.
			
			//std::cout << curr << " not translucent objects:" << std::endl;
			
			pcol = sample.color;
			
			ray.from = hit->P;
			ray.dir = wo;
			ray.tmin = MIN_RAYDIST;
			ray.tmax = -1.0;
			++nBounces;
			
		}
		++curr;
		//if(curr % pbStep == 0) pb->update();
		//done = (curr >= nPhotons) ? true : false;
		//if (inCount % pbStep == 0) pb->update();
		done = (inCount >= nPhotons) ? true : false;
	}
	pb->done();
	pb->setTag("SSS photon map built.");
	
	std::cout << "In photons: " << inCount << std::endl;
	std::cout << "absorbed photons: " << absorbCount << std::endl;
	
	delete lightPowerD;
	
	return true;
}

float RD(float sig_s, float sig_a, float g, float IOR, float r)
{
	//float rd = 1.f/(4*M_PI);
	float rd = 0.25*M_1_PI;//1.f/(4*M_PI);
	float sig_s_ = (1.f-g)*sig_s;
	float sig_t_ = sig_a + sig_s_;
	float alpha_ = sig_s_/sig_t_;
	float sig_tr = sqrtf(3*sig_a*sig_t_);
	float z_r = 1.f/sig_t_;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	float z_v = z_r*(1+4.f*A/3.f);
	
	//float dr = sqrtf(r*r + z_r*z_r);
	//float dv = sqrtf(r*r + z_v*z_v);
	float dr = fSqrt(r*r + z_r*z_r);//sqrtf(r*r + z_r*z_r);
	float dv = fSqrt(r*r + z_v*z_v);//sqrtf(r*r + z_v*z_v);
	
	rd *= alpha_;
	//float real = z_r*(sig_tr+1/dr)*expf(-1.f*sig_tr*dr)/(dr*dr);
	//float vir = z_v*(sig_tr+1/dv)*expf(-1.f*sig_tr*dv)/(dv*dv);
	float real = z_r*(sig_tr+1/dr)*fExp(-1.f*sig_tr*dr)/(dr*dr);//expf(-1.f*sig_tr*dr)/(dr*dr);
	float vir = z_v*(sig_tr+1/dv)*fExp(-1.f*sig_tr*dv)/(dv*dv);//expf(-1.f*sig_tr*dv)/(dv*dv);
	
	rd *= (real+vir);
	
	return rd;
}

color_t dipole(const photon_t& inPhoton, const surfacePoint_t &sp, const vector3d_t &wo, float IOR, float g, const color_t &sigmaS, const color_t &sigmaA )
{
	color_t rd(0.f);
	const color_t Li = inPhoton.c;
	const vector3d_t wi = inPhoton.direction();
	const vector3d_t No = sp.N;
	const vector3d_t Ni = inPhoton.hitNormal;
	
	float cosWiN = wi*Ni;
	
	float Kr_i, Kt_i, Kr_o, Kt_o;
	fresnel(wi, Ni, IOR, Kr_i, Kt_i);
	fresnel(wo, No, IOR, Kr_o, Kt_o);
	
	vector3d_t v = inPhoton.pos-sp.P;
	float r  = v.length()*mcIntegrator_t::sssScale;
	
	// compute RD
	rd.R = M_1_PI*cosWiN*Li.R*Kt_i*Kt_o*RD(sigmaS.R, sigmaA.R, g, IOR, r);
	rd.G = M_1_PI*cosWiN*Li.G*Kt_i*Kt_o*RD(sigmaS.G, sigmaA.G, g, IOR, r);
	rd.B = M_1_PI*cosWiN*Li.B*Kt_i*Kt_o*RD(sigmaS.B, sigmaA.B, g, IOR, r);
	
	
	//std::cout << "Kt_i=" << Kt_i << "    Kt_o=" << Kt_o << std::endl;
	//std::cout << "r=" << r << "    rd=" << RD(sigmaS.R, sigmaA.R, g, IOR, r) << std::endl;
	//std::cout << "Li=" << Li << "   r=" << r  << "  rd=" << rd << std::endl << std::endl;
	
	return rd;
}

color_t dipole2(const photon_t& inPhoton, const surfacePoint_t &sp, const vector3d_t &wo, float IOR, float g, const color_t &sigmaS, const color_t &sigmaA )
{
	color_t rd(0.25*M_1_PI);
	const color_t Li = inPhoton.c;
	const vector3d_t wi = inPhoton.direction();
	const vector3d_t No = sp.N;
	const vector3d_t Ni = inPhoton.hitNormal;
	
	float cosWiN = wi*Ni;
	
	float Kr_i, Kt_i, Kr_o, Kt_o;
	fresnel(wi, Ni, IOR, Kr_i, Kt_i);
	fresnel(wo, No, IOR, Kr_o, Kt_o);
	
	vector3d_t v = inPhoton.pos-sp.P;
	float r  = v.length()*mcIntegrator_t::sssScale;
	
	
	//float rd = 0.25*M_1_PI;//1.f/(4*M_PI);
	color_t sig_s_ = (1.f-g)*sigmaS;
	color_t sig_t_ = sigmaA + sig_s_;
	color_t alpha_ = sig_s_/sig_t_;
	color_t sig_tr = colorSqrt(3*sigmaA*sig_t_);
	
	//std::cout << "sig_t=" << sig_t_ << "  alpha=" << alpha_ << "  sig_tr=" << sig_tr << std::endl;
	
	float z_r = inPhoton.sourceDepth*mcIntegrator_t::sssScale;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	float z_v = z_r*(1+4.f*A/3.f);
	
	point3d_t vSourcePos = inPhoton.sourcePos + inPhoton.hitNormal*(inPhoton.sourceDepth*(2+4.f*A/3.f));
	vector3d_t vecRS = sp.P - inPhoton.sourcePos;
	vector3d_t vecVS = sp.P - vSourcePos;
	
	float dr = vecRS.length()*mcIntegrator_t::sssScale;
	float dv = vecVS.length()*mcIntegrator_t::sssScale;
	
	//std::cout << "dr=" << dr << "  dv=" << dv << std::endl;
	
	//float dr = sqrtf(r*r + z_r*z_r);
	//float dv = sqrtf(r*r + z_v*z_v);
	//float dr = fSqrt(r*r + z_r*z_r);//sqrtf(r*r + z_r*z_r);
	//float dv = fSqrt(r*r + z_v*z_v);//sqrtf(r*r + z_v*z_v);
	
	rd *= alpha_;
	//float real = z_r*(sig_tr+1/dr)*expf(-1.f*sig_tr*dr)/(dr*dr);
	//float vir = z_v*(sig_tr+1/dv)*expf(-1.f*sig_tr*dv)/(dv*dv);
	color_t real = z_r*(sig_tr+1/dr)*colorExp(-1.f*sig_tr*dr)/(dr*dr);//expf(-1.f*sig_tr*dr)/(dr*dr);
	color_t vir = z_v*(sig_tr+1/dv)*colorExp(-1.f*sig_tr*dv)/(dv*dv);//expf(-1.f*sig_tr*dv)/(dv*dv);
	
	rd *= (real+vir);
	
	if (rd.energy() > 1.0) 
	{
		std::cout << "rd=" << rd << " out pos = " << sp.P << "   photon pos=" << inPhoton.sourcePos << " dr=" << dr << "  dv=" << dv << std::endl;
		std::cout << "old RD = " << RD(sigmaS.R, sigmaA.R, g, IOR, r) << std::endl;
	}
	
	rd = rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	
	// compute RD
	//rd.R = cosWiN*Li.R*Kt_i*Kt_o*RD(sigmaS.R, sigmaA.R, g, IOR, r)/M_PI;
	//rd.G = cosWiN*Li.G*Kt_i*Kt_o*RD(sigmaS.G, sigmaA.G, g, IOR, r)/M_PI;
	//rd.B = cosWiN*Li.B*Kt_i*Kt_o*RD(sigmaS.B, sigmaA.B, g, IOR, r)/M_PI;
	
	
	//std::cout << "Kt_i=" << Kt_i << "    Kt_o=" << Kt_o << std::endl;
	//std::cout << "r=" << r << "    rd=" << RD(sigmaS.R, sigmaA.R, g, IOR, r) << std::endl;
	//std::cout << "Li=" << Li << "   r=" << r  << "  rd=" << rd << std::endl << std::endl;
	
	return rd;
}

color_t dipole3(const photon_t& inPhoton, const surfacePoint_t &sp, const vector3d_t &wo, float IOR, float g, const color_t &sigmaS, const color_t &sigmaA )
{
	color_t rd(0.25*M_1_PI);
	const color_t Li = inPhoton.c;
	const vector3d_t wi = inPhoton.direction();
	const vector3d_t No = sp.N;
	const vector3d_t Ni = inPhoton.hitNormal;
	
	float cosWiN = wi*Ni;
	
	float Kr_i, Kt_i, Kr_o, Kt_o;
	fresnel(wi, Ni, IOR, Kr_i, Kt_i);
	fresnel(wo, No, IOR, Kr_o, Kt_o);
	
	vector3d_t v = inPhoton.pos-sp.P;
	float r  = v.length()*mcIntegrator_t::sssScale;
	
	
	//float rd = 0.25*M_1_PI;//1.f/(4*M_PI);
	color_t sig_s_ = (1.f-g)*sigmaS;
	color_t sig_t_ = sigmaA + sig_s_;
	color_t alpha_ = sig_s_/sig_t_;
	color_t sig_tr = colorSqrt(3*sigmaA*sig_t_);
	
	//std::cout << "sig_t=" << sig_t_ << "  alpha=" << alpha_ << "  sig_tr=" << sig_tr << std::endl;
	
	color_t z_r = 1.f/sig_t_/mcIntegrator_t::sssScale;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	color_t z_v = z_r*(1+1.333333333f*A);
	
	point3d_t rSourcePosR = inPhoton.pos + inPhoton.hitNormal*-1*z_r.R;
	point3d_t rSourcePosG = inPhoton.pos + inPhoton.hitNormal*-1*z_r.G;
	point3d_t rSourcePosB = inPhoton.pos + inPhoton.hitNormal*-1*z_r.B;
	
	
	point3d_t vSourcePosR = inPhoton.pos + inPhoton.hitNormal*z_v.R;
	point3d_t vSourcePosG = inPhoton.pos + inPhoton.hitNormal*z_v.G;
	point3d_t vSourcePosB = inPhoton.pos + inPhoton.hitNormal*z_v.B;
	
	
	
	color_t dr, dv;
	dr.R = (sp.P-rSourcePosR).length()*mcIntegrator_t::sssScale;
	dr.G = (sp.P-rSourcePosG).length()*mcIntegrator_t::sssScale;
	dr.B = (sp.P-rSourcePosB).length()*mcIntegrator_t::sssScale;
	
	dv.R = (sp.P-vSourcePosR).length()*mcIntegrator_t::sssScale;
	dv.G = (sp.P-vSourcePosG).length()*mcIntegrator_t::sssScale;
	dv.B = (sp.P-vSourcePosB).length()*mcIntegrator_t::sssScale;
	 
	
	//std::cout << "dr=" << dr << "  dv=" << dv << std::endl;
	
	//float dr = sqrtf(r*r + z_r*z_r);
	//float dv = sqrtf(r*r + z_v*z_v);
	//float dr = fSqrt(r*r + z_r*z_r);//sqrtf(r*r + z_r*z_r);
	//float dv = fSqrt(r*r + z_v*z_v);//sqrtf(r*r + z_v*z_v);
	
	z_r = z_r*mcIntegrator_t::sssScale;
	z_v = z_v*mcIntegrator_t::sssScale;
	
	//color_t dr = colorSqrt(r*r + z_r*z_r);
	//color_t dv = colorSqrt(r*r + z_v*z_v);
	
	rd *= alpha_;
	//float real = z_r*(sig_tr+1/dr)*expf(-1.f*sig_tr*dr)/(dr*dr);
	//float vir = z_v*(sig_tr+1/dv)*expf(-1.f*sig_tr*dv)/(dv*dv);
	color_t real = z_r*(sig_tr+1/dr)*colorExp(-1.f*sig_tr*dr)/(dr*dr);//expf(-1.f*sig_tr*dr)/(dr*dr);
	color_t vir = z_v*(sig_tr+1/dv)*colorExp(-1.f*sig_tr*dv)/(dv*dv);//expf(-1.f*sig_tr*dv)/(dv*dv);
	
	rd *= (real+vir);
	
	if (rd.energy() > 1.0) 
	{
		std::cout << "rd=" << rd << " out pos = " << sp.P << "   photon pos=" << inPhoton.sourcePos << " dr=" << dr << "  dv=" << dv << std::endl;
		std::cout << "old RD = " << RD(sigmaS.R, sigmaA.R, g, IOR, r) << std::endl;
	}
	
	rd = rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	
	// compute RD
	//rd.R = cosWiN*Li.R*Kt_i*Kt_o*RD(sigmaS.R, sigmaA.R, g, IOR, r)/M_PI;
	//rd.G = cosWiN*Li.G*Kt_i*Kt_o*RD(sigmaS.G, sigmaA.G, g, IOR, r)/M_PI;
	//rd.B = cosWiN*Li.B*Kt_i*Kt_o*RD(sigmaS.B, sigmaA.B, g, IOR, r)/M_PI;
	
	
	//std::cout << "Kt_i=" << Kt_i << "    Kt_o=" << Kt_o << std::endl;
	//std::cout << "r=" << r << "    rd=" << RD(sigmaS.R, sigmaA.R, g, IOR, r) << std::endl;
	//std::cout << "Li=" << Li << "   r=" << r  << "  rd=" << rd << std::endl << std::endl;
	
	return rd;
}

// dipole and quadpole

color_t dipoleAdnQuadpole(const photon_t& inPhoton, const surfacePoint_t &sp, const vector3d_t &wo, float IOR, float g, const color_t &sigmaS, const color_t &sigmaA )
{
	color_t rd(0.25*M_1_PI);
	color_t qd(1.f);
	
	const color_t Li = inPhoton.c;
	const vector3d_t wi = inPhoton.direction();
	const vector3d_t No = sp.N;
	const vector3d_t Ni = inPhoton.hitNormal;
	
	
	float gamma = acosf(dot(No, Ni));
	
	float cosWiN = wi*Ni;
	
	float Kr_i, Kt_i, Kr_o, Kt_o;
	fresnel(wi, Ni, IOR, Kr_i, Kt_i);
	fresnel(wo, No, IOR, Kr_o, Kt_o);
	
	vector3d_t v = inPhoton.pos-sp.P;
	float r  = v.length()*mcIntegrator_t::sssScale;
	
	color_t sig_s_ = (1.f-g)*sigmaS;
	color_t sig_t_ = sigmaA + sig_s_;
	color_t alpha_ = sig_s_/sig_t_;
	color_t sig_tr = colorSqrt(3*sigmaA*sig_t_);
	
	color_t z_r = 1.f/sig_t_/mcIntegrator_t::sssScale;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	color_t z_v = z_r*(1+1.333333333f*A);
	
	point3d_t rSourcePosR = inPhoton.pos + inPhoton.hitNormal*-1*z_r.R;
	point3d_t rSourcePosG = inPhoton.pos + inPhoton.hitNormal*-1*z_r.G;
	point3d_t rSourcePosB = inPhoton.pos + inPhoton.hitNormal*-1*z_r.B;
	point3d_t vSourcePosR = inPhoton.pos + inPhoton.hitNormal*z_v.R;
	point3d_t vSourcePosG = inPhoton.pos + inPhoton.hitNormal*z_v.G;
	point3d_t vSourcePosB = inPhoton.pos + inPhoton.hitNormal*z_v.B;
	
//	point3d_t vmSourcePosR = rSourcePosR + 2*(((sp.P-rSourcePosR)*No+0.66666667f*A/sig_t_.R/sssScale)*No);
//	point3d_t vmSourcePosG = rSourcePosG + 2*(((sp.P-rSourcePosG)*No+0.66666667f*A/sig_t_.G/sssScale)*No);
//	point3d_t vmSourcePosB = rSourcePosB + 2*(((sp.P-rSourcePosB)*No+0.66666667f*A/sig_t_.B/sssScale)*No);
//	point3d_t rmSourcePosR = vSourcePosR + 2*(((sp.P-vSourcePosR)*No+0.66666667f*A/sig_t_.R/sssScale)*No);
//	point3d_t rmSourcePosG = vSourcePosG + 2*(((sp.P-vSourcePosG)*No+0.66666667f*A/sig_t_.G/sssScale)*No);
//	point3d_t rmSourcePosB = vSourcePosB + 2*(((sp.P-vSourcePosB)*No+0.66666667f*A/sig_t_.B/sssScale)*No);
	
	point3d_t mInPosR = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*No+0.66666667f*A/sig_t_.R/mcIntegrator_t::sssScale)*No);
	point3d_t mInPosG = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*No+0.66666667f*A/sig_t_.G/mcIntegrator_t::sssScale)*No);
	point3d_t mInPosB = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*No+0.66666667f*A/sig_t_.B/mcIntegrator_t::sssScale)*No);
	
	color_t mr;
	mr.R = (sp.P-mInPosR).length()*mcIntegrator_t::sssScale;
	mr.G = (sp.P-mInPosG).length()*mcIntegrator_t::sssScale;
	mr.B = (sp.P-mInPosB).length()*mcIntegrator_t::sssScale;
	
	vector3d_t iToOR = ((sp.P-rSourcePosR)*No)*No;
	vector3d_t iToOG = ((sp.P-rSourcePosG)*No)*No;
	vector3d_t iToOB = ((sp.P-rSourcePosB)*No)*No;
	color_t xr;
	xr.R = iToOR.length()*mcIntegrator_t::sssScale;
	xr.G = iToOG.length()*mcIntegrator_t::sssScale;
	xr.B = iToOB.length()*mcIntegrator_t::sssScale;
	color_t xv = xr + 1.333333333f*A/sig_t_;
	
	z_r = z_r*mcIntegrator_t::sssScale;
	z_v = z_v*mcIntegrator_t::sssScale;
	
	color_t dr = colorSqrt(r*r + z_r*z_r);
	color_t dv = colorSqrt(r*r + z_v*z_v);
	
	color_t drm, dvm;
//	drm.R = (sp.P-rmSourcePosR).length()*sssScale;
//	drm.G = (sp.P-rmSourcePosG).length()*sssScale;
//	drm.B = (sp.P-rmSourcePosB).length()*sssScale;
//	dvm.R = (sp.P-vmSourcePosR).length()*sssScale;
//	dvm.G = (sp.P-vmSourcePosG).length()*sssScale;
//	dvm.B = (sp.P-vmSourcePosB).length()*sssScale;
	
	dvm = colorSqrt(mr*mr+z_r*z_r);
	drm = colorSqrt(mr*mr+z_v*z_v);
	
	
	rd *= alpha_;
	
	color_t real = z_r*(sig_tr+1/dr)*colorExp(-1.f*sig_tr*dr)/(dr*dr);
	color_t vir = z_v*(sig_tr+1/dv)*colorExp(-1.f*sig_tr*dv)/(dv*dv);
	
	rd *= (real+vir);
	
	qd = xr*(1+sig_tr*dr)*colorExp(-1*sig_tr*dr)*0.125*M_1_PI/(dr*dr*dr)
		+ xr*(1+sig_tr*dv)*colorExp(-1*sig_tr*dv)*0.125*M_1_PI/(dv*dv*dv)
		+ xv*(1+sig_tr*drm)*colorExp(-1*sig_tr*drm)*0.125*M_1_PI/(drm*drm*drm)
		+ xv*(1+sig_tr*dvm)*colorExp(-1*sig_tr*dvm)*0.125*M_1_PI/(dvm*dvm*dvm); 
	
	if (qd.energy() > 1.0) 
	{
		std::cout << "qd=" << qd << std::endl;
		std::cout << "  pcol = " << Li << "   sig_tr = " << sig_tr	<< std::endl;
		std::cout << "  In pos = " << inPhoton.pos << "  out pos = " << sp.P << "  No = " << No << "  Ni = " << Ni << std::endl;
		std::cout << "  rSourcePosB = " << rSourcePosB << "   vSourcePosB = " << vSourcePosB << std::endl;
		//std::cout << "  vmSourcePosB = " << vmSourcePosB << "   rmSourcePosB = " << rmSourcePosB << std::endl;
		//std::cout << "  (sp.P-rSourcePosB)*No = " << (sp.P-rSourcePosB)*No << std::endl;
		//std::cout << "  (0.66666667f*A/sig_t_.R/sssScale) = " << 0.66666667f*A/sig_t_.B/sssScale << std::endl;
		//std::cout << "  ((sp.P-rSourcePosB)*No+0.66666667f*A/sig_t_.B/sssScale)*No = " << ((sp.P-rSourcePosB)*No+0.66666667f*A/sig_t_.B/sssScale)*No << std::endl;
		std::cout << "  xr = " << xr << "   xv = " << xv << std::endl;
		std::cout << "  drm = " << drm << "    dvm = " << dvm << " dr= " << dr << "  dv= " << dv << std::endl << std::endl;
	}
	
	
	color_t result;
	
	result = rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	
	/*if (gamma <= 0.5*M_PI && gamma >=0) {
		result = 2*M_1_PI*(0.5*M_PI-gamma)*rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
		result = result + 2*M_1_PI*gamma*qd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	}
	else {
		result = rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	}*/
	
	return result;
}

color_t dipoleAdnQuadpole2(const photon_t& inPhoton, const surfacePoint_t &sp, const vector3d_t &wo, float IOR, float g, const color_t &sigmaS, const color_t &sigmaA )
{
	color_t rd(0.25*M_1_PI);
	color_t qd(1.f);
	
	const color_t Li = inPhoton.c;
	const vector3d_t wi = inPhoton.direction();
	const vector3d_t No = sp.N;
	const vector3d_t Ni = inPhoton.hitNormal;
	
	
	float gamma = acosf(dot(No, Ni));
	
	float cosWiN = wi*Ni;
	
	float Kr_i, Kt_i, Kr_o, Kt_o;
	fresnel(wi, Ni, IOR, Kr_i, Kt_i);
	fresnel(wo, No, IOR, Kr_o, Kt_o);
	
	vector3d_t v = inPhoton.pos-sp.P;
	float r  = v.length()*mcIntegrator_t::sssScale;
	
	color_t sig_s_ = (1.f-g)*sigmaS;
	color_t sig_t_ = sigmaA + sig_s_;
	color_t alpha_ = sig_s_/sig_t_;
	color_t sig_tr = colorSqrt(3*sigmaA*sig_t_);
	
	color_t z_r = 1.f/sig_t_/mcIntegrator_t::sssScale;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	color_t z_v = z_r*(1+1.333333333f*A);
	
	point3d_t rSourcePosR = inPhoton.pos + inPhoton.hitNormal*-1*z_r.R;
	point3d_t rSourcePosG = inPhoton.pos + inPhoton.hitNormal*-1*z_r.G;
	point3d_t rSourcePosB = inPhoton.pos + inPhoton.hitNormal*-1*z_r.B;
	point3d_t vSourcePosR = inPhoton.pos + inPhoton.hitNormal*z_v.R;
	point3d_t vSourcePosG = inPhoton.pos + inPhoton.hitNormal*z_v.G;
	point3d_t vSourcePosB = inPhoton.pos + inPhoton.hitNormal*z_v.B;
	
		// compute the intersect diection of the two faces
	
	vector3d_t refDir;
	vector3d_t intersectDir = Ni^No;
	if ( intersectDir.length() < 1e-6 )
	{
		if (Ni*No >= 0) {
			refDir = sp.P - inPhoton.pos;
		}
		else {
			refDir = No;
		}

	}
	else
		refDir = intersectDir^Ni;
	refDir.normalize();
	
	if ((sp.P-inPhoton.pos)*refDir < 0) {
		refDir *= -1.f;
	}
	
//	point3d_t vmSourcePosR = rSourcePosR + 2*(((sp.P-rSourcePosR)*refDir+0.66666667f*A/sig_t_.R/sssScale)*refDir);
//	point3d_t vmSourcePosG = rSourcePosG + 2*(((sp.P-rSourcePosG)*refDir+0.66666667f*A/sig_t_.G/sssScale)*refDir);
//	point3d_t vmSourcePosB = rSourcePosB + 2*(((sp.P-rSourcePosB)*refDir+0.66666667f*A/sig_t_.B/sssScale)*refDir);
//	point3d_t rmSourcePosR = vSourcePosR + 2*(((sp.P-vSourcePosR)*refDir+0.66666667f*A/sig_t_.R/sssScale)*refDir);
//	point3d_t rmSourcePosG = vSourcePosG + 2*(((sp.P-vSourcePosG)*refDir+0.66666667f*A/sig_t_.G/sssScale)*refDir);
//	point3d_t rmSourcePosB = vSourcePosB + 2*(((sp.P-vSourcePosB)*refDir+0.66666667f*A/sig_t_.B/sssScale)*refDir);
	
	point3d_t mInPosR = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*refDir+0.66666667f*A/sig_t_.R/mcIntegrator_t::sssScale)*refDir);
	point3d_t mInPosG = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*refDir+0.66666667f*A/sig_t_.G/mcIntegrator_t::sssScale)*refDir);
	point3d_t mInPosB = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*refDir+0.66666667f*A/sig_t_.B/mcIntegrator_t::sssScale)*refDir);
	
	color_t mr;
	mr.R = (sp.P-mInPosR).length()*mcIntegrator_t::sssScale;
	mr.G = (sp.P-mInPosG).length()*mcIntegrator_t::sssScale;
	mr.B = (sp.P-mInPosB).length()*mcIntegrator_t::sssScale;
	
	vector3d_t iToOR = ((sp.P-rSourcePosR)*refDir)*refDir;
	vector3d_t iToOG = ((sp.P-rSourcePosG)*refDir)*refDir;
	vector3d_t iToOB = ((sp.P-rSourcePosB)*refDir)*refDir;
	color_t xr;
	xr.R = iToOR.length()*mcIntegrator_t::sssScale;
	xr.G = iToOG.length()*mcIntegrator_t::sssScale;
	xr.B = iToOB.length()*mcIntegrator_t::sssScale;
	color_t xv = xr + 1.333333333f*A/sig_t_;
	
	z_r = z_r*mcIntegrator_t::sssScale;
	z_v = z_v*mcIntegrator_t::sssScale;
	
	color_t dr = colorSqrt(r*r + z_r*z_r);
	color_t dv = colorSqrt(r*r + z_v*z_v);
	
	color_t drm, dvm;
//	drm.R = (sp.P-rmSourcePosR).length()*sssScale;
//	drm.G = (sp.P-rmSourcePosG).length()*sssScale;
//	drm.B = (sp.P-rmSourcePosB).length()*sssScale;
//	dvm.R = (sp.P-vmSourcePosR).length()*sssScale;
//	dvm.G = (sp.P-vmSourcePosG).length()*sssScale;
//	dvm.B = (sp.P-vmSourcePosB).length()*sssScale;
	
	dvm = colorSqrt(mr*mr+z_r*z_r);
	drm = colorSqrt(mr*mr+z_v*z_v);
	
	
	rd *= alpha_;
	
	color_t real = z_r*(sig_tr+1/dr)*colorExp(-1.f*sig_tr*dr)/(dr*dr);
	color_t vir = z_v*(sig_tr+1/dv)*colorExp(-1.f*sig_tr*dv)/(dv*dv);
	
	rd *= (real+vir);
	
	qd = xr*(1+sig_tr*dr)*colorExp(-1*sig_tr*dr)*0.125*M_1_PI/(dr*dr*dr)
	+ xr*(1+sig_tr*dv)*colorExp(-1*sig_tr*dv)*0.125*M_1_PI/(dv*dv*dv)
	+ xv*(1+sig_tr*drm)*colorExp(-1*sig_tr*drm)*0.125*M_1_PI/(drm*drm*drm)
	+ xv*(1+sig_tr*dvm)*colorExp(-1*sig_tr*dvm)*0.125*M_1_PI/(dvm*dvm*dvm); 
	
//	color_t temp = z_r*(1+sig_tr*dr)*colorExp(-1*sig_tr*dr)*0.125*M_1_PI/(dr*dr*dr)
//				 + z_v*(1+sig_tr*dv)*colorExp(-1*sig_tr*dv)*0.125*M_1_PI/(dv*dv*dv);
//	
//	if (temp.energy() < qd.energy()) 
//	{
//		std::cout << "qd=" << qd << std::endl;
//		std::cout << "rd=" << temp << std::endl;
//		std::cout << "  pcol = " << Li << "   sig_tr = " << sig_tr	<< std::endl;
//		std::cout << "  In pos = " << inPhoton.pos << "  out pos = " << sp.P << "  No = " << No << "  Ni = " << Ni << std::endl;
//		std::cout << "  refDir = " << refDir << "   (sp.P-inPhoton.pos)*refDir = " << (sp.P-inPhoton.pos)*refDir << std::endl;
//		std::cout << "  mInPosR = " << mInPosR << "   mInPosG = " << mInPosG << "   mInPosB = " << mInPosB << std::endl;
//		//std::cout << "  vmSourcePosB = " << vmSourcePosB << "   rmSourcePosB = " << rmSourcePosB << std::endl;
//		//std::cout << "  (sp.P-rSourcePosB)*No = " << (sp.P-rSourcePosB)*No << std::endl;
//		//std::cout << "  (0.66666667f*A/sig_t_.R/sssScale) = " << 0.66666667f*A/sig_t_.B/sssScale << std::endl;
//		//std::cout << "  ((sp.P-rSourcePosB)*No+0.66666667f*A/sig_t_.B/sssScale)*No = " << ((sp.P-rSourcePosB)*No+0.66666667f*A/sig_t_.B/sssScale)*No << std::endl;
//		std::cout << "  xr = " << xr << "   xv = " << xv << std::endl;
//		std::cout << "  drm = " << drm << "    dvm = " << dvm << " dr= " << dr << "  dv= " << dv << std::endl << std::endl;
//	}

	color_t result;
	
	result = qd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	
//	if (gamma <= 0.5*M_PI && gamma >=0) {
//	 result = 2*M_1_PI*(0.5*M_PI-gamma)*rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
//	 result = result + 2*M_1_PI*gamma*qd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
//	 }
//	 else {
//	 result = rd*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
//	 }
	
	return result;
}

color_t RdQdRm(const photon_t& inPhoton, const surfacePoint_t &sp, const vector3d_t &wo, float IOR, float g, const color_t &sigmaS, const color_t &sigmaA )
{
	
	int m_n = 2;
	
	color_t rd(0.25*M_1_PI);
	color_t qd(1.f);
	color_t rm(0.0f);
	
	const color_t Li = inPhoton.c;
	const vector3d_t wi = inPhoton.direction();
	const vector3d_t No = sp.N;
	const vector3d_t Ni = inPhoton.hitNormal;
	
	
	float gamma = acosf(dot(No, Ni));
	
	float cosWiN = wi*Ni;
	
	float Kr_i, Kt_i, Kr_o, Kt_o;
	fresnel(wi, Ni, IOR, Kr_i, Kt_i);
	fresnel(wo, No, IOR, Kr_o, Kt_o);
	
	vector3d_t v = inPhoton.pos-sp.P;
	float r  = v.length()*mcIntegrator_t::sssScale;
	
	color_t sig_s_ = (1.f-g)*sigmaS;
	color_t sig_t_ = sigmaA + sig_s_;
	color_t alpha_ = sig_s_/sig_t_;
	color_t sig_tr = colorSqrt(3*sigmaA*sig_t_);
	
	color_t z_r = 1.f/sig_t_/mcIntegrator_t::sssScale;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	color_t z_v = z_r*(1+1.333333333f*A);
	
	point3d_t rSourcePosR = inPhoton.pos + inPhoton.hitNormal*-1*z_r.R;
	point3d_t rSourcePosG = inPhoton.pos + inPhoton.hitNormal*-1*z_r.G;
	point3d_t rSourcePosB = inPhoton.pos + inPhoton.hitNormal*-1*z_r.B;
	point3d_t vSourcePosR = inPhoton.pos + inPhoton.hitNormal*z_v.R;
	point3d_t vSourcePosG = inPhoton.pos + inPhoton.hitNormal*z_v.G;
	point3d_t vSourcePosB = inPhoton.pos + inPhoton.hitNormal*z_v.B;
	
	// compute the intersect diection of the two faces
	
	vector3d_t refDir;
	vector3d_t intersectDir = Ni^No;
	if ( intersectDir.length() < 1e-6 )
	{
		if (Ni*No >= 0) {
			refDir = sp.P - inPhoton.pos;
		}
		else {
			refDir = No;
		}
		
	}
	else
		refDir = intersectDir^Ni;
	refDir.normalize();
	
	if ((sp.P-inPhoton.pos)*refDir < 0) {
		refDir *= -1.f;
	}
	
	point3d_t mInPosR = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*refDir+0.66666667f*A/sig_t_.R/mcIntegrator_t::sssScale)*refDir);
	point3d_t mInPosG = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*refDir+0.66666667f*A/sig_t_.G/mcIntegrator_t::sssScale)*refDir);
	point3d_t mInPosB = inPhoton.pos + 2*(((sp.P-inPhoton.pos)*refDir+0.66666667f*A/sig_t_.B/mcIntegrator_t::sssScale)*refDir);
	
	color_t mr;
	mr.R = (sp.P-mInPosR).length()*mcIntegrator_t::sssScale;
	mr.G = (sp.P-mInPosG).length()*mcIntegrator_t::sssScale;
	mr.B = (sp.P-mInPosB).length()*mcIntegrator_t::sssScale;
	
	vector3d_t iToOR = ((sp.P-rSourcePosR)*refDir)*refDir;
	vector3d_t iToOG = ((sp.P-rSourcePosG)*refDir)*refDir;
	vector3d_t iToOB = ((sp.P-rSourcePosB)*refDir)*refDir;
	color_t xr;
	xr.R = iToOR.length()*mcIntegrator_t::sssScale;
	xr.G = iToOG.length()*mcIntegrator_t::sssScale;
	xr.B = iToOB.length()*mcIntegrator_t::sssScale;
	color_t xv = xr + 1.333333333f*A/sig_t_;
	
	z_r = z_r*mcIntegrator_t::sssScale;
	z_v = z_v*mcIntegrator_t::sssScale;
	
	color_t dr = colorSqrt(r*r + z_r*z_r);
	color_t dv = colorSqrt(r*r + z_v*z_v);
	
	color_t drm, dvm;	
	dvm = colorSqrt(mr*mr+z_r*z_r);
	drm = colorSqrt(mr*mr+z_v*z_v);
	
	
	//rd *= alpha_;
	
	color_t real = z_r*(sig_tr+1/dr)*colorExp(-1.f*sig_tr*dr)/(dr*dr);
	color_t vir = z_v*(sig_tr+1/dv)*colorExp(-1.f*sig_tr*dv)/(dv*dv);
	
	rd *= (real+vir);
	
	qd = z_r*(1+sig_tr*dr)*colorExp(-1*sig_tr*dr)*0.125*M_1_PI/(dr*dr*dr)
	+ z_v*(1+sig_tr*dv)*colorExp(-1*sig_tr*dv)*0.125*M_1_PI/(dv*dv*dv)
	+ xv*(1+sig_tr*drm)*colorExp(-1*sig_tr*drm)*0.125*M_1_PI/(drm*drm*drm)
	+ xr*(1+sig_tr*dvm)*colorExp(-1*sig_tr*dvm)*0.125*M_1_PI/(dvm*dvm*dvm); 
	
	// compute rm
	
	color_t thickness(0.0f);
	color_t l = z_r;
	
	thickness += z_r;
	
	thickness.R += fabs((sp.P - rSourcePosR)*No)*mcIntegrator_t::sssScale;
	thickness.G += fabs((sp.P - rSourcePosG)*No)*mcIntegrator_t::sssScale;
	thickness.B += fabs((sp.P - rSourcePosB)*No)*mcIntegrator_t::sssScale;
	
	for (int i=-1*m_n; i<=m_n; i++) {
		z_r = 2*i*(thickness+1.33333333f*A/sig_t_) + l;
		z_v = 2*i*(thickness+1.33333333f*A/sig_t_) - l - 1.3333333f*A/sig_t_;

		dr = colorSqrt(r*r + z_r*z_r);
		dv = colorSqrt(r*r + z_v*z_v);	
		
		rm += ( z_r*(1+sig_tr*dr)*colorExp(-1*sig_tr*dr)*0.25*M_1_PI/(dr*dr*dr) 
			   - z_v*(1+sig_tr*dv)*colorExp(-1*sig_tr*dv)*0.25*M_1_PI/(dv*dv*dv)  );
	}
	
	color_t result(0.0f);
	
	//result = rm*Li*cosWiN*Kt_i*Kt_o*M_1_PI;
	
	if (gamma <= 0.5*M_PI && gamma >=0) {
		result += 2*M_1_PI*(0.5*M_PI-gamma)*rd;
		result += 2*M_1_PI*gamma*qd;
	}
	else if ( gamma > 0.5*M_PI && gamma <= M_PI ) {
		result += 2*M_1_PI*(M_PI - gamma)*qd;
		result += 2*M_1_PI*(gamma - 0.5*M_PI)*rm;
	}
	else {
		result += rd;
	}
	
	result = result*Li*cosWiN*Kt_i*Kt_o;
	return result;
}

color_t mcIntegrator_t::estimateSSSMaps(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo ) const
{
	color_t sum(0.f);
	vector3d_t wi(0.0f);
	const object3d_t* hitObj = sp.object;
	
	std::map<const object3d_t*, photonMap_t*>::const_iterator it = SSSMaps.find(hitObj);
	if ( it == SSSMaps.end() ) {
		return sum;
	}
	photonMap_t* sssMap_t = it->second;
	
	float photonSum = 0;
	it = SSSMaps.begin();
	while (it!=SSSMaps.end())
	{
		photonSum += it->second->nPhotons();
		it++;
	}
	
	BSDF_t bsdfs;
	
	void *o_udat = state.userdata;
	unsigned char userdata[USER_DATA_SIZE];
	state.userdata = (void *)userdata;
	
	const material_t *material = sp.material;
	material->initBSDF(state, sp, bsdfs);
	
	color_t sigma_s, sigma_a, diffuseC;
	float IOR;
	TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
	diffuseC = dat->difC;
	sigma_a = dat->sig_a;
	sigma_s = dat->sig_s;
	IOR = dat->IOR;
	
	// sum all photon in translucent object
	const std::vector<const photon_t*>& photons = sssMap_t->getAllPhotons(sp.P);
	
	for (uint i=0; i<photons.size(); i++) {
		//sum += dipole(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);
		//sum += dipole2(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);
		//sum += dipole3(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);
		//sum += dipoleAdnQuadpole(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);
		//sum += dipoleAdnQuadpole2(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);
		sum += RdQdRm(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);

	}
	
	sum *= sssScale*sssScale/((float)sssMap_t->nPaths());
	
	/*if (sum.energy() > 1.f ) {
		std::cout << sp.P << std::endl;	
	}*/
	
	sum *= diffuseC;
	
	state.userdata = o_udat;
	
	return sum;
	
}

float phaseFunc ( const vector3d_t &wi, const vector3d_t &wo )
{
	float g = 0;
	float cosTheta = wi*wo;
	return (1+3*g*cosTheta)*0.25*M_1_PI;
}

color_t mcIntegrator_t::estimateSSSSingleScattering(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	float stepSize = 0.1f/sssScale;
	color_t singleS(0.0f);
	
	if (wo*sp.N < 0) {
		return singleS;
	}
	
	float t0 = 1e10f, t1 = -1e10f;
//	std::cout << "entry point is " << sp.P << std::endl;
//	std::cout << "dir  is " << -1*wo << std::endl;
	
	// get the material infomation
	void *o_udat = state.userdata;
	unsigned char userdata[USER_DATA_SIZE];
	state.userdata = (void *)userdata;
	
	BSDF_t bsdfs;
	const material_t *material = sp.material;
	material->initBSDF(state, sp, bsdfs);
	
	color_t diffuseC;
	color_t sigma_s, sigma_a, sigma_t;
	float IOR;
	TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
	diffuseC = dat->difC;
	sigma_a = dat->sig_a;
	sigma_s = dat->sig_s;
	sigma_t = sigma_s + sigma_a;
	IOR = dat->IOR;

	float Kr_o, Kt_o;
	fresnel(wo, sp.N, IOR, Kr_o, Kt_o);	
	
	// get the refracted direction
	vector3d_t refDir;
	if (!refract(sp.N, wo, refDir, IOR)) {
		return singleS;
	}
	
	
	ray_t ray;
	ray.from = sp.P;
	ray.dir = refDir;
	ray.tmin = MIN_RAYDIST;
	ray.tmax = -1;
	
	surfacePoint_t hit;
	
	if (!scene->intersect(ray, hit)) {
		return singleS;
	}
	
	float badDist = (hit.P-sp.P).length();
	
	bool ismeetBadFace = false;
	while ( hit.N * refDir < 0 )
	{
		ismeetBadFace = true;
		//std::cout << "not out " << hit.P << std::endl;
		ray.from = hit.P;
		ray.tmin = MIN_RAYDIST;
		ray.tmax = -1;
		if (!scene->intersect(ray, hit)) {
			return singleS;
		}
	}
//	if(ismeetBadFace){
//		std::cout << " bad dist = " << badDist << "  new dist = " << (hit.P-sp.P).length() << std::endl;
//		std::cout << std::endl;
//	}
	t0 = 0;
	t1 = (hit.P-sp.P).length();
	float dist = (t1-t0);			
	float pos = t0 + (*state.prng)()*stepSize;
	int samples = dist/stepSize + 1;
	float currentStep = stepSize;
	int stepLength = 1;
	color_t trTmp(1.f);
	
	color_t stepTau(0.f); 
	for (int stepSample = 0; stepSample < samples; stepSample += stepLength)
	{
		ray_t stepRay(sp.P + (ray.dir * pos), ray.dir, 0, currentStep, 0);
		
		stepTau += sigma_t*currentStep*sssScale;
		
		trTmp = colorExp(-1*stepTau);
		
		//std::cout << "\t scatter point is " << stepRay.from << std::endl;
		//std::cout << "trTmp is " << trTmp << std::endl;
		
		color_t radiance = getTranslucentInScatter(state, stepRay, currentStep);
		//if ( state.pixelNumber == 489949 )
		//	std::cout << "radiance is " << radiance << std::endl;
		
		singleS += trTmp * radiance * sigma_s * currentStep * Kt_o * sssScale;
		
		pos += currentStep;
		if(pos - t0 >= dist)
			break;
	}
	//singleS += 0.001;
	//if (singleS.energy() > 1) {
	//if ( state.pixelNumber == 489949 )
	//{
	//	std::cout << " index = " << state.pixelNumber << std::endl;
	//	std::cout << "singleS is " << singleS << " stepTau = " << stepTau << "    dist = " << dist << std::endl << std::endl;
	//}
	//}
	//std::cout << "singleS is " << singleS << std::endl << std::endl;
//	std::cout << "refracted dir  is " << refDir << std::endl;
//	std::cout << "exit point is " << hit.P << std::endl;
//	std::cout << "the length of ray is " << t1-t0 << std::endl << std::endl;
	
	// restore old render state data
	state.userdata = o_udat;
	
	singleS *= diffuseC;
	
	return singleS;
}


color_t mcIntegrator_t::estimateSSSSingleSImportantSampling(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	//float stepSize = 0.1f/sssScale;
	int	samples = this->nSingleScatterSamples;
	std::vector<float> stepSizes;
	color_t singleS(0.0f);
	
	if (wo*sp.N < 0) {
		return singleS;
	}
	
	float t0 = 1e10f, t1 = -1e10f;
	//	std::cout << "entry point is " << sp.P << std::endl;
	//	std::cout << "dir  is " << -1*wo << std::endl;
	
	// get the material infomation
	void *o_udat = state.userdata;
	unsigned char userdata[USER_DATA_SIZE];
	state.userdata = (void *)userdata;
	
	BSDF_t bsdfs;
	const material_t *material = sp.material;
	material->initBSDF(state, sp, bsdfs);
	
	color_t diffuseC;
	color_t sigma_s, sigma_a, sigma_t;
	float IOR;
	TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
	diffuseC = dat->difC;
	sigma_a = dat->sig_a;
	sigma_s = dat->sig_s;
	sigma_t = sigma_s + sigma_a;
	IOR = dat->IOR;
	
	float Kr_o, Kt_o;
	fresnel(wo, sp.N, IOR, Kr_o, Kt_o);	
	
	// get the refracted direction
	vector3d_t refDir;
	if (!refract(sp.N, wo, refDir, IOR)) {
		return singleS;
	}
	
	ray_t ray;
	ray.from = sp.P;
	ray.dir = refDir;
	ray.tmin = MIN_RAYDIST;
	ray.tmax = -1;
	
	surfacePoint_t hit;
	
	if (!scene->intersect(ray, hit)) {
		return singleS;
	}
	
	float badDist = (hit.P-sp.P).length();
	
	bool ismeetBadFace = false;
	while ( hit.N * refDir < 0 )
	{
		ismeetBadFace = true;
		//std::cout << "not out " << hit.P << std::endl;
		ray.from = hit.P;
		ray.tmin = MIN_RAYDIST;
		ray.tmax = -1;
		if (!scene->intersect(ray, hit)) {
			return singleS;
		}
	}
	//if(ismeetBadFace){
	//	std::cout << " bad dist = " << badDist << "  new dist = " << (hit.P-sp.P).length() << std::endl;
	//	std::cout << std::endl;
	//}
	
	// get the light transmit distance
	t0 = 0;
	t1 = (hit.P-sp.P).length();
	float dist = (t1-t0);
	//important sampling
	float sigT = sigma_t.energy();
	float range = 1.f - exp(-1*dist*sigT);
	if (range == 1) {
		range -= 1e-6;
	}
//	if (state.pixelNumber == 70280) {
//		std::cout << "sigma_t = " << sigma_t << std::endl;
//		std::cout << "dist = " << dist << std::endl;
//		std::cout << "range = " << range << std::endl;
//	}
	float lastSamplePos = 0.0f, currSamplePos;
	for (int i=1; i<=samples; i++) {
		currSamplePos = -1.f*log(1-range*(float)i/(float)samples)/sigT;
		stepSizes.push_back(currSamplePos - lastSamplePos);
		lastSamplePos = currSamplePos;
	}
	
	//std::cout << " dist = " << dist << std::endl;
	
	//float pos = t0 + (*state.prng)()*stepSizes[0];
	float pos = t0 + 0.5*stepSizes[0];
	
	//int samples = dist/stepSize + 1;
	float currentStep = stepSizes[0];
	
	color_t trTmp(1.f);
	color_t stepTau(0.f); 
	
//	if ((state.pixelNumber == 223580) || (state.pixelNumber == 223600) )
//	{
//		std::cout << state.pixelNumber << std::endl;
//		std::cout << "Normal = " << sp.N << std::endl;
//	}
	
	for (int stepSample = 0; stepSample < samples; stepSample++)
	{
		currentStep = stepSizes[stepSample];
		
		ray_t stepRay(sp.P + (ray.dir * pos), ray.dir, 0, currentStep, 0);
		
		stepTau += sigma_t*currentStep*sssScale;
		
		trTmp = colorExp(-1*stepTau);
		
		//std::cout << "\t scatter point is " << stepRay.from << std::endl;
		
		color_t inScatter = getTranslucentInScatter(state, stepRay, currentStep);
		
		singleS += trTmp * inScatter * sigma_s * currentStep * Kt_o * sssScale;
		
		//std::cout << state.pixelNumber << std::endl;
//		if (state.pixelNumber == 70280) {
//			std::cout << "trTmp = " << trTmp << "   inScatter = " << inScatter << "  currentStep == " << currentStep << "  kt_o = " << Kt_o << std::endl;
//			std::cout << "singleS = " << singleS << std::endl;
//		}
//		else if (state.pixelNumber == 223600) {
//			std::cout << "trTmp = " << trTmp << "   inScatter = " << inScatter << "  currentStep == " << currentStep << "  kt_o = " << Kt_o << std::endl;
//			std::cout << "singleS = " << singleS << std::endl;
//		}
		
		pos += currentStep;
		
		if(pos - t0 >= dist)
			break;
	}
	//if ((state.pixelNumber == 70280) || (state.pixelNumber == 223600) )
	//	std::cout << std::endl;
	//	std::cout << "refracted dir  is " << refDir << std::endl;
	//	std::cout << "exit point is " << hit.P << std::endl;
	//	std::cout << "the length of ray is " << t1-t0 << std::endl << std::endl;
	
	// restore old render state data
	state.userdata = o_udat;
	
//	if (state.pixelNumber == 70280) {
//		singleS = color_t(1.0,0.0,0.0);
//	}
	
	singleS *= diffuseC;
	
	return singleS;
}

color_t mcIntegrator_t::getTranslucentInScatter(renderState_t& state, ray_t& stepRay, float currentStep) const
{	
	void *o_udat = state.userdata;
	unsigned char userdata[USER_DATA_SIZE];
	state.userdata = (void *)userdata;
	color_t sigma_s, sigma_a, sigma_t;
	float IOR;
	
	color_t inScatter(0.f);
	color_t diffuseC;
	surfacePoint_t sp;
	sp.P = stepRay.from;
	
	ray_t lightRay;
	lightRay.from = sp.P;
	
	ray_t outRay;
	surfacePoint_t outHit;
	
	for(std::vector<light_t *>::const_iterator l=lights.begin(); l!=lights.end(); ++l)
	{
		color_t lcol(0.0);
		
		if( (*l)->diracLight() )
		{
			if( (*l)->illuminate(sp, lcol, lightRay) )
			{
				//std::cout << "\t\t the light ray is from " << lightRay.from << "  dir = " << lightRay.dir << std::endl; 
				
				// get the exit point;
				outRay.from = sp.P;
				outRay.dir = lightRay.dir;
				outRay.tmin = 0;
				outRay.tmax = -1;
				if (!scene->intersect(outRay, outHit)) {
					continue;
				}
				// get the material infomation
				
				BSDF_t bsdfs;
				const material_t *material = outHit.material;
				material->initBSDF(state, outHit, bsdfs);
			
				TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
				diffuseC = dat->difC;
				sigma_a = dat->sig_a;
				sigma_s = dat->sig_s;
				sigma_t = sigma_s + sigma_a;
				IOR = dat->IOR;
				
				point3d_t exitP = outHit.P;
				float cosWi = fabs(outHit.N*outRay.dir);
				float dist = (exitP - sp.P).length();//*cosWi/sqrtf((1.f-(1.f/(float)IOR)*(1.f/(float)IOR))*(1-cosWi*cosWi));
				
				
				float Kr_i, Kt_i;
				fresnel(outRay.dir, outHit.N, IOR, Kr_i, Kt_i);
				
				//std::cout << "\t\t the light ray hit point is " << outHit.P << " and real dist = " << (exitP - sp.P).length() << "  approximation = " << dist << std::endl; 

				// ...shadowed...
				lightRay.from = exitP;
				lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
				if (lightRay.tmax < 0.f) lightRay.tmax = 1e10; // infinitely distant light
				bool shadowed = scene->isShadowed(state, lightRay);
				if (!shadowed)
				{
					color_t lightTr(0.0f);
					
					color_t lightstepTau(0.f);
					
					lightstepTau = sigma_t * dist * sssScale;
					
					lightTr = colorExp(-1*lightstepTau);
					
					//std::cout << "\t\t tau = " << lightstepTau << " and light tau = " << lightTr << std::endl;//
					
					inScatter += (lightTr * lcol * diffuseC * Kt_i) * phaseFunc(lightRay.dir, -1*stepRay.dir);
					
					//std::cout << "\t\t lcol = " << lcol << " contribute= " << (lightTr * lcol * Kt_i) << std::endl;
				}
			}
		} // end if diractLight
		
		else // area light and suchlike
		{
			//std::cout << "area light " << std::endl;
			int n = (*l)->nSamples() >> 2; // samples / 4
			if (n < 1) n = 1;
			float iN = 1.f / (float)n; // inverse of n
			color_t ccol(0.0);
			color_t lightTr(0.0f);
			lSample_t ls;
			
			for(int i=0; i<n; ++i)
			{
				// ...get sample val...
				ls.s1 = (*state.prng)();
				ls.s2 = (*state.prng)();
				
				if((*l)->illumSample(sp, ls, lightRay))
				{
					//if ( state.pixelNumber == 489949 )
					//	std::cout << "\t sample " << i << " lightRay.tmax = " << lightRay.tmax << std::endl;
					// get the exit point;
					outRay.from = sp.P;
					outRay.dir = lightRay.dir;
					outRay.tmin = 0;
					outRay.tmax = -1;
					if (!scene->intersect(outRay, outHit)) {
						continue;
					}
					// get the material infomation
					BSDF_t bsdfs;
					const material_t *material = outHit.material;
					material->initBSDF(state, outHit, bsdfs);
					
					if (!(bsdfs & BSDF_TRANSLUCENT)) {
						continue;
					}
					
					TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
					sigma_a = dat->sig_a;
					sigma_s = dat->sig_s;
					sigma_t = sigma_s + sigma_a;
					IOR = dat->IOR;
					
					
					point3d_t exitP = outHit.P;
					float cosWi = fabs(outHit.N*outRay.dir);
					float dist = (exitP - sp.P).length();//*cosWi/sqrtf((1.f-(1.f/(float)IOR)*(1.f/(float)IOR))*(1-cosWi*cosWi));
					//if ( state.pixelNumber == 489949 )
					//	std::cout << "\t sigma_t="<< sigma_t << "   IOR=" << IOR << "  cosWi="<<cosWi << "   " << (1.f-(1.f/(float)IOR)*(1.f/(float)IOR)) << std::endl;
					
					float Kr_i, Kt_i;
					fresnel(outRay.dir, outHit.N, IOR, Kr_i, Kt_i);
					
					// ...shadowed...
					lightRay.from = exitP;
					lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
					lightRay.tmax -= (exitP - sp.P).length();
					if (lightRay.tmax < 0.f) lightRay.tmax = 1e10; // infinitely distant light
					bool shadowed = scene->isShadowed(state, lightRay);
					
					//std::cout << "sample " << i << " isshadowed = " << shadowed << std::endl;
					//if ( state.pixelNumber == 489949 )
					//	std::cout << "dist = " << dist << "   kt_i=" << Kt_i << std::endl;
					
					if(!shadowed) {
						ccol += ls.col / ls.pdf;
						
						color_t lightstepTau = sigma_t * dist * sssScale;
						
						lightTr += colorExp(-1*lightstepTau)*Kt_i*phaseFunc(lightRay.dir, -1*stepRay.dir);;
					}
				}
			} // end of area light sample loop
			
			lightTr *= iN;
			
			//if ( state.pixelNumber == 489949 )
			//	std::cout << "\t lightTr = " << lightTr << std::endl;
			
			ccol = ccol * iN;
			//if ( state.pixelNumber == 489949 )
			//	std::cout << "\t ccol = " << ccol << std::endl;
			
			inScatter += lightTr * ccol;
		} // end of area lights loop
		
	}
	
	inScatter *= phaseFunc(lightRay.dir, -1*stepRay.dir);
	
	//inScatter *= 30.f;
	
	state.userdata = o_udat;
	return inScatter;
}

color_t mcIntegrator_t::estimateSSSSingleScatteringPhotons(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	float stepSize = 1.f/sssScale;
	color_t singleS(0.0f);
	
	// get the volumetric photonmap
	const object3d_t* hitObj = sp.object;
	std::map<const object3d_t*, photonMap_t*>::const_iterator it = SSSMaps.find(hitObj);
	if ( it == SSSMaps.end() ) {
		return singleS;
	}
	photonMap_t* sssMap_t = it->second;
	
	float t0 = 1e10f, t1 = -1e10f;
	//	std::cout << "entry point is " << sp.P << std::endl;
	//	std::cout << "dir  is " << -1*wo << std::endl;
	
	// get the material infomation
	void *o_udat = state.userdata;
	unsigned char userdata[USER_DATA_SIZE];
	state.userdata = (void *)userdata;
	
	BSDF_t bsdfs;
	const material_t *material = sp.material;
	material->initBSDF(state, sp, bsdfs);
	
	color_t diffuseC;
	color_t sigma_s, sigma_a, sigma_t;
	float IOR;
	TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
	diffuseC = dat->difC;
	sigma_a = dat->sig_a;
	sigma_s = dat->sig_s;
	sigma_t = sigma_s + sigma_a;
	IOR = dat->IOR;
	
	float Kr_o, Kt_o;
	fresnel(wo, sp.N, IOR, Kr_o, Kt_o);	
	
	// get the refracted direction
	vector3d_t refDir;
	if (!refract(sp.N, wo, refDir, IOR)) {
		return singleS;
	}
	
	
	ray_t ray;
	ray.from = sp.P;
	ray.dir = refDir;
	ray.tmin = MIN_RAYDIST;
	ray.tmax = -1;
	
	surfacePoint_t hit;
	
	if (!scene->intersect(ray, hit)) {
		return singleS;
	}
	
	float badDist = (hit.P-sp.P).length();
	
	bool ismeetBadFace = false;
	while ( hit.N * refDir < 0 )
	{
		ismeetBadFace = true;
		//std::cout << "not out " << hit.P << std::endl;
		ray.from = hit.P;
		ray.tmin = MIN_RAYDIST;
		ray.tmax = -1;
		if (!scene->intersect(ray, hit)) {
			return singleS;
		}
	}
//	if(ismeetBadFace){
//		std::cout << " bad dist = " << badDist << "  new dist = " << (hit.P-sp.P).length() << std::endl;
//		std::cout << std::endl;
//	}
	t0 = 0;
	t1 = (hit.P-sp.P).length();
	float dist = (t1-t0);			
	float pos = t0 + (*state.prng)()*stepSize;
	int samples = dist/stepSize + 1;
	float currentStep = stepSize;
	int stepLength = 1;
	color_t trTmp(1.f);
	
	float singleSRadius = 4.f / sssScale;
	
	color_t stepTau(0.f); 
	for (int stepSample = 0; stepSample < samples; stepSample += stepLength)
	{
		ray_t stepRay(sp.P + (ray.dir * pos), ray.dir, 0, currentStep, 0);
		
		stepTau += sigma_t*currentStep*sssScale;
		
		trTmp = colorExp(-1*stepTau);
		
		//std::cout << "\t scatter point is " << stepRay.from << std::endl;
		
		//singleS += trTmp * getTranslucentInScatter(state, stepRay, currentStep) * sigma_s * currentStep * Kt_o * sssScale;
		// here it use sssMaps to evaluate the single scattering
		
		// copied from estimateCaustic
		{
			foundPhoton_t *gathered = (foundPhoton_t *)alloca(nCausSearch * sizeof(foundPhoton_t));
			int nGathered = 0;
			
			float gRadiusSquare = singleSRadius * singleSRadius;
			
			nGathered = sssMap_t->gather(sp.P, gathered, nCausSearch, gRadiusSquare);
			
			//std::cout << "sample Idx " << stepSample << std::endl;
			//std::cout << "find photon " << nGathered << std::endl;
			
			color_t sum(0.f);
			
			if(nGathered > 0)
			{
				float k = 0.75f*M_1_PI/(singleSRadius*singleSRadius*singleSRadius);
				const photon_t *photon;
				
				for(int i=0; i<nGathered; ++i)
				{
					photon = gathered[i].photon;
					
					color_t attenuation = colorExp(-1*(photon->pos-photon->sourcePos).length()*sssScale*sigma_t);
					
					sum += photon->color() * phaseFunc(-1*photon->direction(), -1*stepRay.dir );
				}
				sum *= k / sssMap_t->nPaths();
				sum = sum / sigma_s;
				//std::cout << " photon color " << sum << std::endl << std::endl;
			}
			
			singleS += trTmp * sum * sigma_s * currentStep * Kt_o * sssScale;
			
		}			
		
		pos += currentStep;
		
		if(pos - t0 >= dist)
			break;
	}
	//	std::cout << "refracted dir  is " << refDir << std::endl;
	//	std::cout << "exit point is " << hit.P << std::endl;
	//	std::cout << "the length of ray is " << t1-t0 << std::endl << std::endl;
	
	// restore old render state data
	state.userdata = o_udat;
	
	singleS *= diffuseC;
	
	return singleS;
}

__END_YAFRAY
