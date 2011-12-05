/****************************************************************************
 *
 * 		irradianceCache.cc: icTree, icRecord and hemisphere types and
 *      operators implementation.
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

#include <yafraycore/irradianceCache.h>

#if HAVE_XML
#include <libxml/parser.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#define MY_ENCODING "ISO-8859-1"
#endif

__BEGIN_YAFRAY

// stratifiedHemisphere METHODS
// ***********************************************************************
stratifiedHemisphere::stratifiedHemisphere(int nm):M(nm), N(M_PI * M), hal(2) {
	hal.setStart(time(0));
	vk = new vector3d_t[N];
	vkMinus = new vector3d_t[N];
	uk = new vector3d_t[N];
	tanTheta = new float[M];
	sinTheta = new float[M];
	sinThetaMinus = new float[M];
	cosTheta = new float[M];
	cosThetaMinus = new float[M];
	cosThetaPlus = new float[M];

	calcSinThetas();
	calcCosThetas();
	calcTanThetas();
	calcSinThetaMinuses();
	calcCosThetaMinuses();
	calcCosThetaPluses();
	calcVks();
	calcUks();
	calcVkMinuses();
}

stratifiedHemisphere::stratifiedHemisphere(const stratifiedHemisphere &strat): hal(2) {
	hal.setStart(time(0));
	if (this != &strat) {
		if (&strat == NULL) {
			M = N = 0;
		} else {
			M = strat.getM();
			N = strat.getN();
			vk = new vector3d_t[N];
			vkMinus = new vector3d_t[N];
			uk = new vector3d_t[N];
			tanTheta = new float[M];
			sinTheta = new float[M];
			sinThetaMinus = new float[M];
			cosTheta = new float[M];
			cosThetaMinus = new float[M];
			cosThetaPlus = new float[M];
			calcSinThetas();
			calcCosThetas();
			calcTanThetas();
			calcSinThetaMinuses();
			calcCosThetaMinuses();
			calcCosThetaPluses();
			calcVks();
			calcUks();
			calcVkMinuses();
		}
	}
}

stratifiedHemisphere::~stratifiedHemisphere() {
	delete[] vk;
	delete[] vkMinus;
	delete[] uk;
	delete[] tanTheta;
	delete[] sinTheta;
	delete[] sinThetaMinus;
	delete[] cosTheta;
	delete[] cosThetaMinus;
	delete[] cosThetaPlus;
}

stratifiedHemisphere & stratifiedHemisphere::operator=(const stratifiedHemisphere &strat) {
	hal.setStart(time(0));
	if (&strat != NULL && this != &strat) {
		if (M != strat.getM()) {
			M = strat.getM();
			delete tanTheta;
			delete sinTheta;
			delete sinThetaMinus;
			delete cosTheta;
			delete cosThetaMinus;
			delete cosThetaPlus;
			tanTheta = new float[M];
			sinTheta = new float[M];
			sinThetaMinus = new float[M];
			cosTheta = new float[M];
			cosThetaMinus = new float[M];
			cosThetaPlus = new float[M];
			calcSinThetas();
			calcCosThetas();
			calcTanThetas();
			calcSinThetaMinuses();
			calcCosThetaMinuses();
			calcCosThetaPluses();
		}
		if (N != strat.getN()) {
			N = strat.getN();
			delete vk;
			delete vkMinus;
			delete uk;
			vk = new vector3d_t[N];
			vkMinus = new vector3d_t[N];
			uk = new vector3d_t[N];
			calcVks();
			calcUks();
			calcVkMinuses();
		}
	}
	return *this;
}

vector3d_t stratifiedHemisphere::getDirection(int j, int k, unsigned int r) {
	if (j<0 || j>M || k<0 || k>N)
		Y_INFO << "ERROR(stratifiedHemisphere.getDirection): j, k out of bound" << std::endl;
	float s1 = RI_vdC(j+k+r);
	float s2 = hal.getNext();
	float tmp = ((float)j+s1)/(float)M;
	float sinTheta = fSqrt(tmp);
	float phi = M_2PI*((float)k+s2)/N;
	return vector3d_t(sinTheta * fCos(phi),
					  sinTheta * fSin(phi),
					  fSqrt(1 - tmp));
}

vector3d_t stratifiedHemisphere::getDirection(int j, int k, float s1, float s2) {
	if (j<0 || j>M || k<0 || k>N)
		Y_INFO << "ERROR(stratifiedHemisphere.getDirection): j, k out of bound" << std::endl;
	float tmp = ((float)j+s1)/(float)M;
	float sinTheta = fSqrt(tmp);
	float phi = M_2PI*(k+s2)/N;
	return vector3d_t(sinTheta * fCos(phi),
					  sinTheta * fSin(phi),
					  fSqrt(1 - tmp));
}

void stratifiedHemisphere::calcUks() {
	for (int k=0; k<N; k++) {
		float phi = M_2PI * ((float)k + 0.5f) / (float)N;
		uk[k] = vector3d_t(fCos(phi), fSin(phi), 0.f);
	}
}

void stratifiedHemisphere::calcVks() {
	for (int k=0; k<N; k++) {
		float phi = M_2PI * ((float)k + 0.5f) / (float)N;
		vk[k] = vector3d_t(-fSin(phi), fCos(phi), 0.f);
	}
}

void stratifiedHemisphere::calcVkMinuses() {
	for (int k=0; k<N; k++) {
		float phi = M_2PI * (float)k / (float)N;
		vkMinus[k] = vector3d_t(-fSin(phi), fCos(phi), 0.f);
	}
}

void stratifiedHemisphere::calcTanThetas() {
	for (int j=0; j<M; j++) {
		//tanTheta[j] = fTan(fAsin( fSqrt( ((float)j+0.5f) / (float) M ) ));
		tanTheta[j] = fSqrt( ((float)j+0.5f) / ((float)M - (float)j - 0.5) );
	}
}

void stratifiedHemisphere::calcSinThetas() {
	for (int j=0; j<M; j++) {
		sinTheta[j] = fSqrt( ((float)j + 0.5f) / (float)M );
	}
}

void stratifiedHemisphere::calcSinThetaMinuses() {
	for (int j=0; j<M; j++) {
		sinThetaMinus[j] = fSqrt( (float)j / (float)M );
	}
}

void stratifiedHemisphere::calcCosThetas() {
	for (int j=0; j<M; j++) {
		cosTheta[j] = fSqrt( 1.0f - ((float)j + 0.5f) / (float)M );
	}
}

void stratifiedHemisphere::calcCosThetaMinuses() {
	for (int j=0; j<M; j++) {
		cosThetaMinus[j] = fSqrt( 1.0f - (float)j / (float)M );
	}
}

void stratifiedHemisphere::calcCosThetaPluses() {
	for (int j=0; j<M; j++) {
		cosThetaPlus[j] = fSqrt( 1.0f - ((float)j + 1.0f) / (float)M );
	}
}

// icREC_t METHODS
// ***********************************************************************

const float icRec_t::NORMALIZATION_TERM = 8.113140441;


//! number of total sections are nSamples = pi*m^2
icRec_t::icRec_t(float kappa, stratifiedHemisphere *strat):
		stratHemi(strat), kappa(kappa) {
	r = std::numeric_limits<float>::max();
}

icRec_t::icRec_t(float kappa, const surfacePoint_t &sp, stratifiedHemisphere *strat):
		surfacePoint_t(sp), stratHemi(strat), kappa(kappa) {
	r = std::numeric_limits<float>::max();
}

vector3d_t icRec_t::getSampleHemisphere(int j, int k, unsigned int r) {
	return changeBasis( stratHemi->getDirection(j, k, r), NU, NV, Nup );
}

vector3d_t icRec_t::getSampleHemisphere(int j, int k, float s1, float s2) {
	return changeBasis( stratHemi->getDirection(j, k, s1, s2), NU, NV, Nup );
}

void icRec_t::changeSampleRadius(float newr) {
	// we use minimal distance radius (without clamping for now)
	if (newr < r) {
		r = newr;
	}
}

float icRec_t::getWeight(const icRec_t &record) const {
	float dot = Nup * record.getNup();
	// if record is pointing to the other side, better not to count his contribution
	if (dot<0.f) {
		return 0.f;
	}
	float epNor = fSqrt(1.f - dot) * NORMALIZATION_TERM;
	float epPos = (P - record.P).length() * 2.f / rClamp;
	float weight = 1.f - kappa * std::max(epPos, epNor);
	return weight;
}

bound_t icRec_t::getBound() const {
	return bound_t(P - rClamp, P + rClamp);
}

void icRec_t::setPixelArea(const diffRay_t &ray) {
	spDifferentials_t diff(*this, ray);
	pArea = fSqrt(diff.projectedPixelArea());
	rMin = 3.f * pArea;
	rMax = 20.f * pArea;
}

void icRec_t::setNup(const vector3d_t &wo) {
	Nup = FACE_FORWARD(Ng, N, wo);
}

bool icRec_t::inFront(const icRec_t &record) const {
	float di = (P - record.P) * ((Nup + record.getNup() )/2.0f);
	if (di < -0.01f) // small negative value, ¿it works?
		return true;
	return false;
}

void icRec_t::clampRbyGradient() {
	rClamp = std::min(r, std::min(irr.R/transGrad[0].length(),
								  std::min( irr.G/transGrad[1].length(), irr.B/transGrad[2].length())) );
}

void icRec_t::clampRbyScreenSpace() {
	rClamp = std::min(std::max(rClamp, rMin), rMax);
}

void icRec_t::clampRbyNeighbor() {
	rClamp = std::min(rClamp, rNeighbor);
}

void icRec_t::clampGradient() {
	for (int i=0; i<3; i++) {
		transGrad[i] =
				transGrad[i] * std::min(1.f, r/rMin);
	}
}

void icRec_t::setRNeighbor(float r) {
	rNeighbor = r;
}

//
// ***********************************************************************

void icTree_t::add(icRec_t *rec) {
	lock.writeLock();
	const bound_t &bound = rec->getBound();
	recursiveAdd(&root, treeBound, rec, bound,
				 2*rec->getRadius()*M_SQRT3 ); // 2*r*sqrt(3) = (bound.a - bound.g).length
	lock.unlock();
	totalRecords++;
}

bool icTree_t::icLookup_t::operator()(const point3d_t &p, const icRec_t *sample) { // point p isn't used
	if (!record->inFront(*sample)) {
		float weight = sample->getWeight(*record);
		if (weight > 0.f) { // TODO: see if weight > 0 is correct or should be a small number
			// get weighted irradiance sample = E_i(p) * w_i(p)
			// E_i(p) = E_i + (n_i x n) * drotE_i
			color_t rotGradResult, transGradResult;
			vector3d_t NCross = record->getNup() ^ sample->getNup(); // should be normalized already
			rotGradResult.R = NCross * sample->rotGrad[0];
			rotGradResult.G = NCross * sample->rotGrad[1];
			rotGradResult.B = NCross * sample->rotGrad[2];
			vector3d_t posDif = record->P - sample->P;
			transGradResult.R = posDif * sample->transGrad[0];
			transGradResult.G = posDif * sample->transGrad[1];
			transGradResult.B = posDif * sample->transGrad[2];
			radSamples.push_back( (sample->irr + rotGradResult + transGradResult) * weight );
			totalWeight += weight;
		}
	} else {
	//	Y_INFO << "In front!" << std::endl;
	}
	return true; // when could it be false? example?
}

void icTree_t::recursiveFindNear(octNode_t<icRec_t *> *node, bound_t &nodeBound,
								 const icRec_t *record, std::vector<icRec_t *> &nearRecs,
								 float &minR) {
	for (unsigned int i = 0; i < node->data.size(); ++i) {
		// pass the "in front" test
		if (!record->inFront( *(node->data[i]) ) ) {
			float distance = (record->P - (node->data[i])->P).length();
			// if both radius overlaps
			if ( distance <= (record->r + (node->data[i])->r) ) {
				float rSum = (node->data[i])->r + distance;
				// checks for triangule inequality
				if ( rSum < record->r ){
					minR = std::min( minR, rSum );
					// add pointer to record to nearRecs
					nearRecs.push_back(node->data[i]);
				}
			}
		}
	}
	bound_t dataBound(record->P-record->r, record->P+record->r);
	// check on all the childrens that the records radius overlap
	point3d_t center = nodeBound.center();
	// Determine which children the item overlaps
	bool over[8];
	over[1] = over[3] = over[5] = over[7] = (dataBound.a.x <= center.x);
	over[0] = over[2] = over[4] = over[6] = (dataBound.g.x  > center.x);
	if(dataBound.a.y > center.y)  over[2] = over[3] = over[6] = over[7] = false;
	if(dataBound.g.y <= center.y) over[0] = over[1] = over[4] = over[5] = false;
	if(dataBound.a.z > center.z)  over[4] = over[5] = over[6] = over[7] = false;
	if(dataBound.g.z <= center.z) over[0] = over[1] = over[2] = over[3] = false;
	for (int child = 0; child < 8; ++child)
	{
		// don't do anything if radius not overlap or child node doent exist
		if (!over[child] || !node->children[child]) continue;
		// compute childbound and keep searching in child
		bound_t childBound;
		childBound.a.x = (child & 1) ? nodeBound.a.x : center.x;
		childBound.g.x = (child & 1) ? center.x : nodeBound.g.x;
		childBound.a.y = (child & 2) ? nodeBound.a.y : center.y;
		childBound.g.y = (child & 2) ? center.y : nodeBound.g.y;
		childBound.a.z = (child & 4) ? nodeBound.a.z : center.z;
		childBound.g.z = (child & 4) ? center.z : nodeBound.g.z;
		recursiveFindNear(node->children[child], childBound, record, nearRecs, minR);
	}
}

void icTree_t::neighborClamp(icRec_t *record) {
	// perform lock in tree, so we can't add other records at the same time
	lock.readLock();
	// if record is outside the scene don't do anything
	if (!treeBound.includes(record->P)) return;
	// create the vector to store all the neighbors
	std::vector<icRec_t *> nearRecs;
	// set neighbor clamped radius equal to the original distance to surfaces
	float minR = record->r;
	// search for all the near records
	recursiveFindNear(&root, treeBound, record, nearRecs, minR);
	// perform neighbor clamp on record
	record->setRNeighbor(minR);
	record->clampRbyNeighbor();
	// perform neighbor clamp to neighbors
	for (unsigned int i=0; i<nearRecs.size(); i++) {
		nearRecs[i]->setRNeighbor(minR + (nearRecs[i]->P - record->P).length());
		nearRecs[i]->clampRbyNeighbor();
	}
	// TODO: neighbor clamp algorith
	lock.unlock();
}

bool icTree_t::getIrradiance(icRec_t *record) {
	icLookup_t lookupProc(record);
	lookup(record->P, lookupProc); // ads weighted radiance values to vector
	// if there is no good irradiance sample return false
	if (lookupProc.radSamples.size() == 0) {
		return false;
	}
	// calculate Sum(E_i(p) * w_i(p))
	for (unsigned int i=0; i<lookupProc.radSamples.size(); i++) {
		record->irr += lookupProc.radSamples[i];
	}
	record->irr = record->irr / lookupProc.totalWeight; // E(p) = Sum(E_i(p) * w_i(p)) / Sum(w_i(p))
	return true;
}

void icTree_t::saveToXml(const std::string &fileName) {
	int rc;
	xmlTextWriterPtr writer;
	// Create a new XmlWriter for uri
	writer = xmlNewTextWriterFilename(fileName.c_str(), 1);
	if (writer == NULL) {
		Y_INFO << "testXmlwriterFilename: Error creating the xml writer" << std::endl;
		return;
	}
	// Start the document
	rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
	// Create root element named ICtree
	xmlTextWriterStartElement(writer, BAD_CAST "ICtree");
	xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "boundMin",
										   "%f,%f,%f", treeBound.a.x, treeBound.a.y, treeBound.a.z );
	xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "boundMax",
										   "%f,%f,%f", treeBound.g.x, treeBound.g.y, treeBound.g.z );
	octNode_t<icRec_t *> *nodes[maxDepth+1];
	int sibling[maxDepth+1];
	int level = 0;
	nodes[0] = &root;
	sibling[0] = 8; // end condition
	sibling[1] = 0;

	// search for the next valid node
	do {
		//  go down one level
		level++;
		nodes[level] = nodes[level-1]->children[sibling[level]];
		sibling[level+1] = 0;
		// if not a valid node
		if (!nodes[level]) {
			bool invalid = true;
			// go up until can move to next sibling
			do {
				sibling[level]++;
				level--;
				// The first one is for invalid nodes.
				if (invalid) {
					invalid = false;
				} else {
					// Close ICNode node
					rc = xmlTextWriterEndElement(writer);
					if (rc < 0) {
						Y_INFO << "testXmlwriterFilename: Error at closing ICNode \n" << std::endl;
						return;
					}
				}

			} while (sibling[level+1]==8);
		}
		else // find a valid one
		{
			// PROCESS DATA
			int size = nodes[level]->data.size();
			for (int i=0; i<size; i++) {
				icRec_t *record = nodes[level]->data[i];
				// Create Record node
				xmlTextWriterStartElement(writer, BAD_CAST "ICRecord");
				// Save members
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "pos", "%f,%f,%f", record->P.x, record->P.y, record->P.z );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "r", "%f", record->r );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "rMin", "%f", record->rMin );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "rMax", "%f", record->rMax );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "rNeighbor", "%f", record->rNeighbor );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "rClamp", "%f", record->getRadius() );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "irr", "%f,%f,%f", record->irr.R, record->irr.G, record->irr.B );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "nUp", "%f,%f,%f", record->getNup().x, record->getNup().y, record->getNup().z );
				xmlTextWriterStartElement(writer, BAD_CAST "RotGrad");
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "r", "%f,%f,%f", record->rotGrad[0].x, record->rotGrad[0].y, record->rotGrad[0].z );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "g", "%f,%f,%f", record->rotGrad[1].x, record->rotGrad[1].y, record->rotGrad[1].z );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "b", "%f,%f,%f", record->rotGrad[2].x, record->rotGrad[2].y, record->rotGrad[2].z );
				xmlTextWriterEndElement(writer);
				xmlTextWriterStartElement(writer, BAD_CAST "TransGrad");
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "r", "%f,%f,%f", record->transGrad[0].x, record->transGrad[0].y, record->transGrad[0].z );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "g", "%f,%f,%f", record->transGrad[1].x, record->transGrad[1].y, record->transGrad[1].z );
				xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "b", "%f,%f,%f", record->transGrad[2].x, record->transGrad[2].y, record->transGrad[2].z );
				xmlTextWriterEndElement(writer);
				// Close Record node
				xmlTextWriterEndElement(writer);
			}
			// Create ICNode node
			rc = xmlTextWriterStartElement(writer, BAD_CAST "ICNode");
			if (rc < 0) {
				Y_INFO << "testXmlwriterFilename: Error at xmlTextWriterStartElement\n" << std::endl;
				return;
			}
			xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "dataSize", "%d", size );
			xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "childId", "%d", sibling[level] );
		}
	} while (level>=0);
	// We already close ICtree element in the last loop
	xmlTextWriterEndDocument(writer);
	xmlFreeTextWriter(writer);
}


__END_YAFRAY
