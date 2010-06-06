/****************************************************************************
 * 			common.cc: common methods for light integrators
 *      This is part of the yafray package
 *      Copyright (C) 2006  Mathias Wein
 *
 *      This library is free software; you can redistribute it and/or
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
#include <core_api/material.h>
#include <core_api/integrator.h>
#include <core_api/light.h>
#include <yafraycore/scr_halton.h>
#include <yafraycore/photon.h>
#include <utilities/mcqmc.h>
#include <utilities/sample_utils.h>
#include <yafraycore/spectrum.h>
#include <integrators/integr_utils.h>

__BEGIN_YAFRAY

struct TranslucentData_t
{
	color_t sig_s;
	color_t sig_a;
	float	IOR;
};

//! estimate direct lighting with multiple importance sampling using the power heuristic with exponent=2
/*! sp.material must be initialized with "initBSDF()" before calling this function! */
color_t estimateDirect_PH(renderState_t &state, const surfacePoint_t &sp, const std::vector<light_t *> &lights, scene_t *scene, const vector3d_t &wo, bool trShad, int sDepth)
{
	color_t col;
	bool shadowed;
	unsigned int l_offs = 0;
	const material_t *material = sp.material;
	ray_t lightRay;
	lightRay.from = sp.P;
	for(std::vector<light_t *>::const_iterator l=lights.begin(); l!=lights.end(); ++l)
	{
		color_t lcol(0.0), scol;
		float lightPdf;
		// handle lights with delta distribution, e.g. point and directional lights
		if( (*l)->diracLight() )
		{
			if( (*l)->illuminate(sp, lcol, lightRay) )
			{
				// ...shadowed...
				lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
				shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
				if(!shadowed)
				{
					if(trShad) lcol *= scol;
					color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
					color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay); // FIXME: add also to the other lightsources!
					col += surfCol * lcol * std::fabs(sp.N*lightRay.dir) * transmitCol;
				}
				//else
				//	col = color_t(1, 0, 0); // make areas visible which are not lit due to being shadowed
			}
		}
		else // area light and suchlike
		{
			Halton hal3(3);
			int n = (*l)->nSamples();
			if(state.rayDivision > 1) n = std::max(1, n/state.rayDivision);
			float invNS = 1.f / (float)n;
			unsigned int offs = n * state.pixelSample + state.samplingOffs + l_offs;
			bool canIntersect=(*l)->canIntersect();//false;
			l_offs += 4567; //just some number to have different sequences per light...and it's a prime even...
			color_t ccol(0.0);
			lSample_t ls;
			hal3.setStart(offs-1);
			for(int i=0; i<n; ++i)
			{
				// ...get sample val...
				ls.s1 = RI_vdC(offs+i);
				ls.s2 = hal3.getNext();
				if(state.rayDivision > 1)
				{
					ls.s1 = addMod1(ls.s1, state.dc1);
					ls.s2 = addMod1(ls.s2, state.dc2);
				}

				if( (*l)->illumSample (sp, ls, lightRay) )
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
				for(int i=0; i<n; ++i)
				{
					ray_t bRay;
					bRay.tmin = MIN_RAYDIST; bRay.from = sp.P;
					float s1 = scrHalton(3, offs+i);
					float s2 = scrHalton(4, offs+i);
					if(state.rayDivision > 1)
					{
						s1 = addMod1(s1, state.dc1);
						s2 = addMod1(s2, state.dc2);
					}
					float W = 0.f;
					sample_t s(s1, s2, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
					color_t surfCol = material->sample(state, sp, wo, bRay.dir, s, W);
					if( s.pdf>1e-6f && (*l)->intersect(bRay, bRay.tmax, lcol, lightPdf) )
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
							CFLOAT cos2 = std::fabs(sp.N*bRay.dir);
							ccol2 += surfCol * lcol * cos2 * w / s.pdf;
						}
					}
				}
				col += ccol2 * invNS;
			}
		}
	} //end light loop
	return col;
}

color_t estimatePhotons(renderState_t &state, const surfacePoint_t &sp, const photonMap_t &map, const vector3d_t &wo, int nSearch, PFLOAT radius)
{
	if(!map.ready()) return color_t(0.f);

	foundPhoton_t *gathered = (foundPhoton_t *)alloca(nSearch * sizeof(foundPhoton_t));
	int nGathered = 0;

	float gRadiusSquare = radius * radius;

	nGathered = map.gather(sp.P, gathered, nSearch, gRadiusSquare);

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
		sum *= 1.f / ( float(map.nPaths()) );
	}
	return sum;
}

bool createCausticMap(const scene_t &scene, const std::vector<light_t *> &all_lights, photonMap_t &cMap, int depth, int count, progressBar_t *pb, std::string intName)
{
	cMap.clear();
	ray_t ray;
	int maxBounces = depth;
	unsigned int nPhotons=count;
	std::vector<light_t *> lights;

	for(unsigned int i=0;i<all_lights.size();++i)
	{
		if(all_lights[i]->shootsCausticP())
		{
			lights.push_back(all_lights[i]);
		}

	}

	int numLights = lights.size();

	if(numLights > 0)
	{
		float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
		float fNumLights = (float)numLights;
		float *energies = new float[numLights];
		for(int i=0;i<numLights;++i) energies[i] = lights[i]->totalEnergy().energy();
		pdf1D_t *lightPowerD = new pdf1D_t(energies, numLights);

		Y_INFO << intName << ": Light(s) photon color testing for caustics map:" << yendl;
		color_t pcol(0.f);
		for(int i=0;i<numLights;++i)
		{
			pcol = lights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
			lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
			Y_INFO << intName << ": Light [" << i+1 << "] Photon col:" << pcol << " | lnpdf: " << lightNumPdf << yendl;
		}

		delete[] energies;

		int pbStep;
		Y_INFO << intName << ": Building caustics photon map..." << yendl;
		pb->init(128);
		pbStep = std::max(1U, nPhotons / 128);
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

			sL = float(curr) / float(nPhotons);

			int lightNum = lightPowerD->DSample(sL, &lightNumPdf);

			if(lightNum >= numLights)
			{
				Y_ERROR << intName << ": lightPDF sample error! " << sL << "/" << lightNum << yendl;
				delete lightPowerD;
				return false;
			}

			color_t pcol = lights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
			ray.tmin = MIN_RAYDIST;
			ray.tmax = -1.0;
			pcol *= fNumLights * lightPdf / lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
			if(pcol.isBlack())
			{
				++curr;
				done = (curr >= nPhotons);
				continue;
			}
			BSDF_t bsdfs = BSDF_NONE;
			int nBounces = 0;
			bool causticPhoton = false;
			bool directPhoton = true;
			const material_t *material = 0;
			const volumeHandler_t *vol = 0;

			while( scene.intersect(ray, *hit2) )
			{
				if(isnan(pcol.R) || isnan(pcol.G) || isnan(pcol.B))
				{
					Y_WARNING << intName << ": NaN (photon color)" << yendl;
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
						cMap.pushPhoton(np);
						cMap.setNumPaths(curr);
					}
				}
				// need to break in the middle otherwise we scatter the photon and then discard it => redundant
				if(nBounces == maxBounces) break;
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
			done = (curr >= nPhotons);
		}
		pb->done();
		pb->setTag("Caustic photon map built.");
		Y_INFO << intName << ": Done." << yendl;
		Y_INFO << intName << ": Shot " << curr << " caustic photons from " << numLights <<" light(s)." << yendl;
		Y_INFO << intName << ": Stored caustic photons: " << cMap.nPhotons() << yendl;

		delete lightPowerD;

		if(cMap.nPhotons() > 0)
		{
			pb->setTag("Building caustic photons kd-tree...");
			cMap.updateTree();
			Y_INFO << intName << ": Done." << yendl;
		}
	}
	else
	{
		Y_INFO << intName << ": No caustic source lights found, skiping caustic map building..." << yendl;
	}
	return true;
}


bool createSSSMaps( const scene_t &scene, const std::vector<light_t *> &lights, std::map<const object3d_t*, photonMap_t*> &SSSMaps, int depth, int count, progressBar_t *pb, std::string intName )
{
	// init and compute light pdf etc.
	ray_t ray;
	int maxBounces = depth;
	unsigned int nPhotons=count;
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
	int pbStep;
	Y_INFO << intName << ": Building SSS photon map..." << yendl;
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
		while( scene.intersect(ray, *hit2) )
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

color_t estimateSSSMaps(renderState_t &state, const surfacePoint_t &sp, const std::map<const object3d_t*, photonMap_t*> &SSSMaps, const vector3d_t &wo )
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
	
	color_t sigma_s, sigma_a;
	float IOR;
	TranslucentData_t* dat = (TranslucentData_t*)state.userdata;
	sigma_a = dat->sig_a;
	sigma_s = dat->sig_s;
	IOR = dat->IOR;
	
	//std::cout << "sigma_a = " << sigma_a.R << "  sigma_s = " << sigma_s.R << "  IOR = " << IOR << std::endl;
	
	// sum all photon in translucent object
	const std::vector<const photon_t*>& photons = sssMap_t->getAllPhotons(sp.P);
	
	//std::cout << "Sample " << state.pixelNumber << "    Get photons number is " << photons.size() << std::endl;
	
	for (uint i=0; i<photons.size(); i++) {
		sum += dipole(*photons[i],sp,wo,IOR,0.f,sigma_s,sigma_a);
	}
	sum *= 100.f/photonSum;//(float)sssMap_t->nPhotons();
	//sum *= 1.f/photons.size();
	std::cout << "sum = " << sum << "" << photonSum / 10.f << std::endl;
	
	state.userdata = o_udat;
	
	return sum;
}

float RD(float sig_s, float sig_a, float g, float IOR, float r)
{
	float rd = 1.f/(4*M_PI);
	float sig_s_ = (1.f-g)*sig_s;
	float sig_t_ = sig_a + sig_s_;
	float alpha_ = sig_s_/sig_t_;
	float sig_tr = sqrtf(3*sig_a*sig_t_);
	float z_r = 1.f/sig_t_;
	float Fdr = -1.440/(IOR*IOR)+0.710/IOR+0.668+0.0636*IOR;
	float A = (1+Fdr)/(1-Fdr);
	float z_v = z_r*(1+4.f*A/3.f);
	float dr = sqrtf(r*r + z_r*z_r);
	float dv = sqrtf(r*r + z_v*z_v);
	
	rd *= alpha_;
	float real = z_r*(sig_tr+1/dr)*expf(-1.f*sig_tr*dr)/(dr*dr);
	float vir = z_v*(sig_tr+1/dv)*expf(-1.f*sig_tr*dv)/(dv*dv);
	
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
	float r  = v.length()*20.f;
	
	// compute RD
	rd.R = cosWiN*Li.R*Kt_i*Kt_o*RD(sigmaS.R, sigmaA.R, g, IOR, r)/M_PI;
	rd.G = cosWiN*Li.G*Kt_i*Kt_o*RD(sigmaS.G, sigmaA.G, g, IOR, r)/M_PI;
	rd.B = cosWiN*Li.B*Kt_i*Kt_o*RD(sigmaS.B, sigmaA.B, g, IOR, r)/M_PI;
	
	
	//std::cout << "Kt_i=" << Kt_i << "    Kt_o=" << Kt_o << std::endl;
	//std::cout << "r=" << r << "    rd=" << RD(sigmaS.R, sigmaA.R, g, IOR, r) << std::endl;
	//std::cout << "Li=" << Li << "  rd=" << rd << std::endl << std::endl;
	
	return rd;
}

__END_YAFRAY
