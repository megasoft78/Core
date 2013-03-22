/****************************************************************************
 *
 * 		irradianceCache.h: icTree, icRecord and hemisphere definitions.
 *      This is part of the yafaray package
 *      Copyright (C) 2010  George Laskowsky Ziguilinsky
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
 *
 */

#ifndef Y_IRRADIANCECACHE_H
#define Y_IRRADIANCECACHE_H

#include <core_api/vector3d.h>
#include <core_api/bound.h>
#include <core_api/color.h>
#include <core_api/surface.h>
#include <core_api/material.h>
#include <core_api/scene.h>
#include <yafraycore/octree.h>
#include <utilities/mathOptimizations.h>
#include <utilities/mcqmc.h>

#include <limits>
#include <ctime>

__BEGIN_YAFRAY

inline vector3d_t changeBasis(const vector3d_t &vec, const vector3d_t &nx, const vector3d_t &ny, const vector3d_t &nz) {
	vector3d_t dir(nx * vec.x);
	dir += ny * vec.y;
	dir += nz * vec.z;
	return dir;
}

//! Stratified hemisphere for getting random direction vectors distributed proportionally to the cosine term.
/*!
  * There is no kind of copy constructor or assignment operator, please, don't do that, just use and discard!
  * TODO: checks for boundries in methods! (0<=j<M, 0<=k<N)
  */
struct YAFRAYCORE_EXPORT stratifiedHemisphere {
	stratifiedHemisphere(int nm);
	stratifiedHemisphere(const stratifiedHemisphere &strat);
	~stratifiedHemisphere();

	stratifiedHemisphere &operator=(const stratifiedHemisphere &strat);

	inline void randomize() { hal.setStart((unsigned int)time(0)); }

    //!< get random direction sample from section j,k in local coordinate system
	vector3d_t getDirection(int j, int k, unsigned int r);
	vector3d_t getDirection(int j, int k, float s1, float s2);

    //!< get vector v_k: base-plane vector in the direction (pi/2, phi_k + pi/2)
	inline const vector3d_t &getVk(int k) const { return vk[k]; }

	//!< get vector v_k: base-plane vector in the direction (pi/2, phi_k + pi/2)
    inline const vector3d_t &getVkMinus(int k) const { return vkMinus[k]; }

    //!< get vector u_k: base-plane vector in the direction (pi/2, phi_k)
	inline const vector3d_t &getUk(int k) const { return uk[k]; }
    //!< get tan(theta_j)
	inline float getTanTheta(int j) const { return tanTheta[j]; }
	inline float getSinTheta(int j) const { return sinTheta[j]; }
	inline float getSinThetaMinus(int j) const { return sinThetaMinus[j]; }
	inline float getCosTheta(int j) const { return cosTheta[j]; }
	inline float getCosThetaMinus(int j) const { return cosThetaMinus[j]; }
	inline float getCosThetaPlus(int j) const { return cosThetaPlus[j]; }
	inline int getM() const { return M; }
	inline int getN() const { return N; }

private:
	int M; //!< number of divisions along theta
	int N; //!< number of divisions along phi

	void calcVks();
	void calcVkMinuses();
	void calcUks();
	void calcTanThetas();
	void calcSinThetas();
	void calcSinThetaMinuses();
	void calcCosThetas();
	void calcCosThetaMinuses();
	void calcCosThetaPluses();

	vector3d_t *vk;
	vector3d_t *vkMinus;
	vector3d_t *uk;
	float *tanTheta;
	float *sinTheta;
	float *sinThetaMinus;
	float *cosTheta;
	float *cosThetaMinus;
	float *cosThetaPlus;
	Halton hal; //!< random number generator
};


//! Record of Irradiance Cache
struct YAFRAYCORE_EXPORT icRec_t : public surfacePoint_t
{
	icRec_t(float kappa, stratifiedHemisphere *strat);
	icRec_t(float kappa, const surfacePoint_t &sp, stratifiedHemisphere *strat);
	// METHODS
	vector3d_t		getSampleHemisphere(int j, int k, unsigned int r); //!< compute indirect light with direct lighting of first bounce
	vector3d_t		getSampleHemisphere(int j, int k, float s1, float s2); //!< compute indirect light with direct lighting of first bounce
	float			getWeight(const icRec_t &record) const;
	inline float	getRadius() const { return rClamp; } //!< return the radius of the sample "action" area
	bound_t			getBound() const; //!< get the bounding box of the sample sphere
	inline int		getM() const { return stratHemi->getM(); } //!< return the number of division if hemisphere along theta
	inline int		getN() const { return stratHemi->getN(); } //!< return the number of division if hemisphere along phi
	void			changeSampleRadius(float newr); //!< change radius if it is a new minimum
	void			setPixelArea(const diffRay_t &dir); //!< calculates the projected pixel area on the surface position of the sample
	void			setNup(const vector3d_t &wo);
	const vector3d_t& getNup() const { return Nup; }
	bool			inFront(const icRec_t &record) const;
	void			clampRbyGradient();
	void			clampRbyScreenSpace();
	void			clampRbyNeighbor();
	void			clampGradient();
	void			setRNeighbor(float r);
	// VARIABLES
	color_t			irr; //!< cached irradiance
	vector3d_t		rotGrad[3];
	vector3d_t		transGrad[3];
	stratifiedHemisphere *stratHemi; //!< sampling hemisphere at point location
	float			r; //!< minimum distance of all rays from hemisphere sampling
	float			rMin; //!< min radius based on screen space (1-3 times projected pixel area)
	float			rMax; //!< max radius based on screen space (20 times projected pixel area)
	float			rNeighbor; //!< saves the neighbor clamped radius
	// ToDo: adaptative sampling
	static const float NORMALIZATION_TERM; //!< from T&L weight function normalization term 1/sqrt(1-cos10°) for 10°
private:
	float			rClamp; //!< radius clamped by projected pixel area and gradients
	vector3d_t		Nup; //!< normal vector on the side of the hitting ray
	float			pArea; //!< projected pixel area over the sample
	float			kappa; //!< overall changing accuaracy constant
};


struct YAFRAYCORE_EXPORT icTree_t : public octree_t<icRec_t *>
{
	icTree_t(const bound_t &bound, int levels, int m):
			octree_t<icRec_t *>(bound, levels), stratHemi(m), totalRecords(0) {}
	//! Get irradiance estimation at point p. Return false if there isn't a cached irradiance sample near.
	bool getIrradiance(icRec_t *record);
	//! Add a new cached irradiance sample
	void add(icRec_t *record);
	//! Perform neighbor clamping on record
	void neighborClamp(icRec_t *record);
	//! Store the entire tree (with IC records data) into an xml file named fileName
	void saveToXml(const std::string &fileName);
	unsigned int getTotalRecords() const { return totalRecords; }
	stratifiedHemisphere stratHemi;
private:
	void recursiveFindNear(octNode_t<icRec_t *> *node, bound_t &nodeBound, const icRec_t *record,
						   std::vector<icRec_t *> &nearRecs, float &minR);
	struct icLookup_t {
		icLookup_t(const icRec_t *rec): record(rec),totalWeight(0.f) {}
		bool operator()(const point3d_t &p, const icRec_t *record);
		std::vector<color_t> radSamples;
		const icRec_t *record;
		float totalWeight;
	};
	unsigned int totalRecords;
};

__END_YAFRAY

#endif // Y_IRRADIANCECACHE_H
