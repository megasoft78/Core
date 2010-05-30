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
//#include <yafraycore/spectrum.h>


/*=============================================================
just a test material for using material class
=============================================================*/

__BEGIN_YAFRAY


struct TranslucentData_t
{
	color_t sig_s;
	color_t sig_a;
	float	IOR;
};

class translucentMat_t: public nodeMaterial_t
{
	public:
		translucentMat_t(color_t diffuseC, color_t siga, color_t sigs, float ior);
		virtual ~translucentMat_t();
		virtual void initBSDF(const renderState_t &state, const surfacePoint_t &sp, unsigned int &bsdfTypes)const;
		virtual color_t eval(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wl, BSDF_t bsdfs)const;
		virtual color_t sample(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, vector3d_t &wi, sample_t &s)const;
		virtual color_t emit(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo)const;
		virtual float pdf(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wi, BSDF_t bsdfs)const;
		static material_t* factory(paraMap_t &params, std::list< paraMap_t > &eparans, renderEnvironment_t &env);
	protected:
		color_t diffuseCol;
		color_t sigma_s;
		color_t sigma_a;
		float	IOR;
};

translucentMat_t::translucentMat_t(color_t diffuseC,color_t siga, color_t sigs, float ior):diffuseCol(diffuseC),sigma_a(siga),sigma_s(sigs),IOR(ior)
{
	bsdfFlags = BSDF_TRANSLUCENT;
	diffuseCol = diffuseC;
	sigma_a = siga;
	sigma_s = sigs;
	IOR = ior;
	//std::cout << sigma_a << " " << sigma_s << " " << IOR << std::endl;
}

translucentMat_t::~translucentMat_t()
{
}

void translucentMat_t::initBSDF(const renderState_t &state, const surfacePoint_t &sp, unsigned int &bsdfTypes)const
{
	bsdfTypes=bsdfFlags;
	
	TranslucentData_t *dat = (TranslucentData_t *)state.userdata;
	//memset(dat, 0, sizeof(TranslucentData_t));

	dat->sig_s = this->sigma_s;
	dat->sig_a = this->sigma_a;
	dat->IOR = this->IOR;
	
}

color_t translucentMat_t::eval(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wl, BSDF_t bsdfs)const
{
	color_t scolor(0.f);
	
	// face forward:
	vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
	//if(!(bsdfs & bsdfFlags )) 
	//	return scolor;
	if(N*wl < 0.0) 
		return scolor;
	//if (bsdfs & BSDF_TRANSLUCENT) 
	{
		scolor = diffuseCol;
	}
	
	float Kr, Kt;
	fresnel(wl, N, IOR, Kr, Kt);
	
	return scolor*Kr;
}

color_t translucentMat_t::sample(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, vector3d_t &wi, sample_t &s)const
{
	color_t scolor(0.f);
	PFLOAT cos_Ng_wo = sp.Ng*wo, cos_Ng_wi;
	
	vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
	
	wi = SampleCosHemisphere(N, sp.NU, sp.NV, s.s1, s.s2);
	s.pdf = std::fabs(wi*N);
	
	cos_Ng_wi = sp.Ng*wi;
	if(cos_Ng_wo*cos_Ng_wi > 0)
		scolor = diffuseCol;
	return scolor;
}

color_t translucentMat_t::emit(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo)const
{
	return color_t(0.f);
}

float translucentMat_t::pdf(const renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo, const vector3d_t &wi, BSDF_t bsdfs)const
{	
	vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
	return std::fabs(wi*N);
}

material_t* translucentMat_t::factory(paraMap_t &params, std::list< paraMap_t > &eparans, renderEnvironment_t &env)
{
	color_t col(1.0f);
	color_t siga(0.01f);
	color_t sigs(1.0f);
	float ior = 1.3;
	params.getParam("color", col);
	params.getParam("sigmaA", siga);
	params.getParam("sigmaS", sigs);
	params.getParam("IOR", ior);
	
	return new translucentMat_t(col,siga,sigs,ior);
}

extern "C"
{
	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("translucent", translucentMat_t::factory);
	}
}

__END_YAFRAY
