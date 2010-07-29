/****************************************************************************
 * 			myDiffuse.cc: a test of diffuse materials
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
 
#include <core_api/material.h>
#include <core_api/environment.h>
#include <core_api/scene.h>
#include <materials/maskmat.h>
#include <materials/microfacet.h>
//#include <yafraycore/spectrum.h>


/*=============================================================
just a test material for using material class
=============================================================*/

__BEGIN_YAFRAY

#define C_TRANSLUCENT 	0
#define C_GLOSSY		1
#define C_DIFFUSE		2

struct TranslucentData_t
{
	color_t difC;
	color_t sig_s;
	color_t sig_a;
	float IOR;
	float mTransl, mDiffuse, mGlossy, pDiffuse;
	
	void *stack;
};

class translucentMat_t: public nodeMaterial_t
{
	public:
		translucentMat_t(color_t diffuseC, color_t specC, color_t glossyC, color_t siga, color_t sigs, float ior, float mT, float mD, float mG, float exp);
		virtual ~translucentMat_t();
		virtual void initBSDF(const renderState_t &state, const surfacePoint_t &sp, unsigned int &bsdfTypes)const;
		virtual color_t eval(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wl, BSDF_t bsdfs)const;
		virtual color_t sample(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, vector3d_t &wi, sample_t &s)const;
		virtual color_t emit(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo)const;
		virtual float pdf(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wi, BSDF_t bsdfs)const;
		//virtual void getSpecular(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo,
		//						 bool &refl, bool &refr, vector3d_t *const dir, color_t *const col)const;
		static material_t* factory(paraMap_t &params, std::list< paraMap_t > &eparans, renderEnvironment_t &env);
	protected:
		shaderNode_t* diffuseS; //!< shader node for diffuse colo
		shaderNode_t* glossyS;
		shaderNode_t* glossyRefS;
		shaderNode_t* bumpS;
	
		color_t diffuseCol;
		color_t specRefCol;
		color_t gloss_color;
	
		float with_diffuse;
		float translucency, diffusity, glossity;
		float pDiffuse;
		float exponent;
	
		BSDF_t cFlags[3];
		int nBSDF;
	
		// parameters for translucent property	
		color_t sigma_s;
		color_t sigma_a;
		float	IOR;
};

translucentMat_t::translucentMat_t(color_t diffuseC, color_t specC, color_t glossyC, color_t siga, color_t sigs, float ior, float mT, float mD, float mG, float exp): diffuseCol(diffuseC),specRefCol(specC),gloss_color(glossyC),
								sigma_a(siga),sigma_s(sigs),IOR(ior),
								translucency(mT), diffusity(mD), glossity(mG), exponent(exp),
								diffuseS(0), glossyS(0), glossyRefS(0), bumpS(0)
{	
//	gloss_color = color_t(1.0f);
//	
//	diffusity = 0.001;
//	glossity = 1.0;
//	translucency = 0.9f;
//	exponent = 800;
	
	cFlags[C_TRANSLUCENT] = (BSDF_TRANSLUCENT);
	cFlags[C_GLOSSY] = (BSDF_GLOSSY | BSDF_REFLECT);
	
	if(diffusity>0)
	{
		cFlags[C_DIFFUSE] = BSDF_DIFFUSE | BSDF_REFLECT;
		with_diffuse = true;
		nBSDF = 3;
	}
	else
	{
		cFlags[C_DIFFUSE] = BSDF_NONE;
		nBSDF = 2;
	}
	
	bsdfFlags = cFlags[C_TRANSLUCENT] | cFlags[C_GLOSSY] | cFlags[C_DIFFUSE];
}

translucentMat_t::~translucentMat_t()
{
}

void translucentMat_t::initBSDF(const renderState_t &state, const surfacePoint_t &sp, unsigned int &bsdfTypes)const
{
	TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
	
	dat->stack = (char*)state.userdata + sizeof(TranslucentData_t);
	nodeStack_t stack(dat->stack);
	if(bumpS) evalBump(stack, state, sp, bumpS);
	
	std::vector<shaderNode_t *>::const_iterator iter, end=allViewindep.end();
	for(iter = allViewindep.begin(); iter!=end; ++iter) (*iter)->eval(stack, state, sp);
	
	dat->difC = diffuseS?diffuseS->getColor(stack):diffuseCol;
	dat->sig_s = this->sigma_s;
	dat->sig_a = this->sigma_a;
	dat->IOR = this->IOR;
	
	dat->mDiffuse = this->diffusity;
	dat->mGlossy = glossyRefS ? glossyRefS->getScalar(stack) : this->glossity;
	dat->mTransl = this->translucency;
	
	dat->pDiffuse = std::min(0.6f , 1.f - (dat->mGlossy/(dat->mGlossy + (1.f-dat->mGlossy)*dat->mDiffuse)) );
	
	bsdfTypes=bsdfFlags;
}

color_t translucentMat_t::eval(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wl, BSDF_t bsdfs)const
{
	if( !(bsdfs & BSDF_DIFFUSE) || ((sp.Ng*wl)*(sp.Ng*wo)) < 0.f ) 
		return color_t(0.f);
	
	TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
	color_t col(0.f);
	bool diffuse_flag = bsdfs & BSDF_DIFFUSE;
	
	nodeStack_t stack(dat->stack);
	vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
	
	float wiN = std::fabs(wl * N);
	float woN = std::fabs(wo * N);
	
	float Kr, Kt;
	fresnel(wl, N, IOR, Kr, Kt);
	
	float mR = (1.0f - Kt*dat->mTransl);
	
	if(bsdfs & BSDF_GLOSSY)
	{
		vector3d_t H = (wo + wl).normalize(); // half-angle
		float cos_wi_H = std::max(0.f, wl*H);
		float glossy;
		
		{
			//glossy = Blinn_D(H*N, exponent) * SchlickFresnel(cos_wi_H, dat->mGlossy) / ASDivisor(cos_wi_H, woN, wiN);
			glossy = mR*Blinn_D(H*N, exponent) * SchlickFresnel(cos_wi_H, dat->mGlossy) / ASDivisor(cos_wi_H, woN, wiN);
		}
		
		col = glossy*(glossyS ? glossyS->getColor(stack) : gloss_color);
		//col = glossy*gloss_color;
	}
	
	if(with_diffuse && diffuse_flag)
	{
		//col += diffuseReflect(wiN, woN, dat->mGlossy, dat->mDiffuse, (diffuseS ? diffuseS->getColor(stack) : diff_color)) * ((orenNayar)?OrenNayar(wi, wo, N):1.f);
		col += mR * diffuseReflect(wiN, woN, dat->mGlossy, dat->mDiffuse, (diffuseS ? diffuseS->getColor(stack) : diffuseCol));
	}
	
	return col;
}

color_t translucentMat_t::sample(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, vector3d_t &wi, sample_t &s)const
{
	TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
	float cos_Ng_wo = sp.Ng*wo;
	float cos_Ng_wi;
	vector3d_t N = (cos_Ng_wo<0) ? -sp.N : sp.N;
	vector3d_t Hs(0.f);
	s.pdf = 0.f;
	float Kr, Kt;
	float wiN = 0.f , woN = 0.f;
	
	fresnel(wi, N, IOR, Kr, Kt);
	
	// missing! get components
	nodeStack_t stack(dat->stack);
	bool use[3] = {false, false, false};
	float sum = 0.f, accumC[3], val[3], width[3];
	int cIndex[3]; // entry values: 0 := specular part, 1 := glossy part, 2:= diffuse part;
	int rcIndex[3]; // reverse fmapping of cIndex, gives position of spec/glossy/diff in val/width array
	accumC[0] = Kt*dat->mTransl;
	accumC[1] = (1.f - Kt*dat->mTransl)*(1.f - dat->pDiffuse);
	accumC[2] = (1.f - Kt*dat->mTransl)*(dat->pDiffuse);
	
	int nMatch = 0, pick = -1;
	for(int i = 0; i < nBSDF; ++i)
	{
		if((s.flags & cFlags[i]) == cFlags[i])
		{
			use[i] = true;
			width[nMatch] = accumC[i];
			cIndex[nMatch] = i;
			rcIndex[i] = nMatch;
			sum += width[nMatch];
			val[nMatch] = sum;
			++nMatch;
		}
	}
	if(!nMatch || sum < 0.00001){ return color_t(0.f); }
	else if(nMatch==1){ pick=0; width[0]=1.f; }
	else
	{
		float inv_sum = 1.f/sum;
		for(int i=0; i<nMatch; ++i)
		{
			val[i] *= inv_sum;
			width[i] *= inv_sum;
			if((s.s1 <= val[i]) && (pick<0 ))	pick = i;
		}
	}
	if(pick<0) pick=nMatch-1;
	float s1;
	if(pick>0) s1 = (s.s1 - val[pick-1]) / width[pick];
	else 	   s1 = s.s1 / width[pick];
	
	color_t scolor(0.f);
	switch(cIndex[pick])
	{
		case C_TRANSLUCENT: // specular reflect
			
			break;
		case C_GLOSSY: // glossy
			Blinn_Sample(Hs, s1, s.s2, exponent);
			break;
		case C_DIFFUSE: // lambertian
		default:
			wi = SampleCosHemisphere(N, sp.NU, sp.NV, s1, s.s2);
			cos_Ng_wi = sp.Ng*wi;
			if(cos_Ng_wo*cos_Ng_wi < 0) return color_t(0.f);
	}
	
	wiN = std::fabs(wi * N);
	woN = std::fabs(wo * N);
	
	if(cIndex[pick] != C_TRANSLUCENT)
	{
		// evaluate BSDFs and PDFs...
		if(use[C_GLOSSY])
		{
			PFLOAT glossy;
			PFLOAT cos_wo_H;
			if(cIndex[pick] != C_GLOSSY)
			{
				vector3d_t H = (wi+wo).normalize();
				Hs = vector3d_t(H*sp.NU, H*sp.NV, H*N);
				cos_wo_H = wo*H;
			}
			else
			{
				vector3d_t H = Hs.x*sp.NU + Hs.y*sp.NV + Hs.z*N;
				cos_wo_H = wo*H;
				if ( cos_wo_H < 0.f )
				{
					H = reflect_plane(N, H);
					cos_wo_H = wo*H;
				}
				// Compute incident direction by reflecting wo about H
				wi = reflect_dir(H, wo);
				cos_Ng_wi = sp.Ng*wi;
				if(cos_Ng_wo*cos_Ng_wi < 0) return color_t(0.f);
			}
			
			wiN = std::fabs(wi * N);
			
			{
				s.pdf += Blinn_Pdf(Hs.z, cos_wo_H, exponent) * width[rcIndex[C_GLOSSY]];
				glossy = Blinn_D(Hs.z, exponent) * SchlickFresnel(cos_wo_H, dat->mGlossy) / ASDivisor(cos_wo_H, woN, wiN);
			}
			scolor = (CFLOAT)glossy*(1.f-Kt*dat->mTransl)*(glossyS ? glossyS->getColor(stack) : gloss_color);
		}
		
		if(use[C_DIFFUSE])
		{
			scolor += (1.f-Kt*dat->mTransl)*diffuseReflect(wiN, woN, dat->mGlossy, dat->mDiffuse, (diffuseS ? diffuseS->getColor(stack) : diffuseCol));
			s.pdf += wiN * width[rcIndex[C_DIFFUSE]];
		}
	}
	s.sampledFlags = cFlags[cIndex[pick]];
	
	return scolor;
}

color_t translucentMat_t::emit(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo)const
{
	return color_t(0.f);
}

float translucentMat_t::pdf(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wi, BSDF_t bsdfs)const
{	
	TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
	bool transmit = ((sp.Ng * wo) * (sp.Ng * wi)) < 0.f;
	if(transmit) return 0.f;
	vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
	float pdf = 0.f;
	CFLOAT Kr, Kt;
	
	fresnel(wi, N, IOR, Kr, Kt);
	
	float accumC[3], sum=0.f, width;
	accumC[0] = Kt*dat->mTransl;
	accumC[1] = (1.f - Kt*dat->mTransl)*(1.f - dat->pDiffuse);
	accumC[2] = (1.f - Kt*dat->mTransl)*(dat->pDiffuse);
	
	int nMatch=0;
	for(int i=0; i<nBSDF; ++i)
	{
		if((bsdfs & cFlags[i]) == cFlags[i])
		{
			width = accumC[i];
			sum += width;
			if(i == C_GLOSSY)
			{
				vector3d_t H = (wi+wo).normalize();
				PFLOAT cos_wo_H = wo*H;
				PFLOAT cos_N_H = N*H;
				pdf += Blinn_Pdf(cos_N_H, cos_wo_H, exponent) * width;
			}
			else if(i == C_DIFFUSE)
			{
				pdf += std::fabs(wi*N) * width;
			}
			++nMatch;
		}
	}
	if(!nMatch || sum < 0.00001) return 0.f;
	return pdf / sum;
}

//void translucentMat_t::getSpecular(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo,
//						 bool &refl, bool &refr, vector3d_t *const dir, color_t *const col)const
//{
//	PFLOAT cos_Ng_wo = sp.Ng*wo, cos_Ng_wi;
//	
//	vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
//	
//	float Kr, Kt;
//	fresnel(wo, N, IOR, Kr, Kt);
//	
//	refr = false;
//	dir[0] = reflect_plane(N, wo);
//	//col[0] = (mirColS ? mirColS->getColor(stack) : specRefCol) * Kr;
//	col[0] = specRefCol * Kr;
//	refl = true;
//}

material_t* translucentMat_t::factory(paraMap_t &params, std::list< paraMap_t > &eparans, renderEnvironment_t &env)
{
	color_t col(1.0f);
	color_t glossyC(1.0f);
	color_t specC(1.0f);
	color_t siga(0.01f);
	color_t sigs(1.0f);
	float ior = 1.3;
	float mT=0.9, mG=1.0, mD=0.001f;
	float exp = 800;
	const std::string *name=0;
	params.getParam("color", col);
	params.getParam("glossy_color", glossyC);
	params.getParam("specular_color", specC);
	params.getParam("sigmaA", siga);
	params.getParam("sigmaS", sigs);
	params.getParam("IOR", ior);
	params.getParam("diffuse_reflect", mD);
	params.getParam("glossy_reflect", mG);
	params.getParam("sss_transmit", mT);
	params.getParam("exponent", exp);
	
	translucentMat_t* mat = new translucentMat_t(col,specC,glossyC,siga,sigs,ior, mT, mD, mG, exp);
	
	std::vector<shaderNode_t *> roots;
	std::map<std::string, shaderNode_t *> nodeList;
	std::map<std::string, shaderNode_t *>::iterator actNode;
	
	// Prepare our node list
	nodeList["diffuse_shader"] = NULL;
	nodeList["glossy_shader"] = NULL;
	nodeList["glossy_reflect_shader"] = NULL;
	nodeList["bump_shader"] = NULL;
	
	if(mat->loadNodes(eparans, env))
	{
		for(actNode = nodeList.begin(); actNode != nodeList.end(); actNode++)
		{
			if(params.getParam(actNode->first, name))
			{
				std::map<std::string,shaderNode_t *>::const_iterator i = mat->shader_table.find(*name);
				
				if(i!=mat->shader_table.end())
				{
					actNode->second = i->second;
					roots.push_back(actNode->second);
				}
				else Y_WARNING << "Glossy: Shader node " << actNode->first << " '" << *name << "' does not exist!" << yendl;
			}
		}
	}
	else Y_ERROR << "Glossy: loadNodes() failed!" << yendl;
	
	mat->diffuseS = nodeList["diffuse_shader"];
	mat->glossyS = nodeList["glossy_shader"];
	mat->glossyRefS = nodeList["glossy_reflect_shader"];
	mat->bumpS = nodeList["bump_shader"];
	
	// solve nodes order
	if(!roots.empty())
	{
		std::vector<shaderNode_t *> colorNodes;
		
		mat->solveNodesOrder(roots);
		
		if(mat->diffuseS) mat->getNodeList(mat->diffuseS, colorNodes);
		if(mat->glossyS) mat->getNodeList(mat->glossyS, colorNodes);
		if(mat->glossyRefS) mat->getNodeList(mat->glossyRefS, colorNodes);
		mat->filterNodes(colorNodes, mat->allViewdep, VIEW_DEP);
		mat->filterNodes(colorNodes, mat->allViewindep, VIEW_INDEP);
		if(mat->bumpS) mat->getNodeList(mat->bumpS, mat->bumpNodes);
	}
	mat->reqMem = mat->reqNodeMem + sizeof(TranslucentData_t);
	
	return mat;
}

extern "C"
{
	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("translucent", translucentMat_t::factory);
	}
}

__END_YAFRAY
