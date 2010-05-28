#ifndef Y_PHOTONKDTREE_H
#define Y_PHOTONKDTREE_H

#include <yafray_config.h>

#include <utilities/y_alloc.h>
#include <core_api/bound.h>
#include <algorithm>
#include <vector>

__BEGIN_YAFRAY

namespace kdtree {

#define PHOTONKD_MAX_STACK 64
#define PHOTON_NON_REC_LOOKUP 1

	//*******************************
	// added by ronnie
	// for include boundingbox infos
	//*******************************
	template <class T>
	struct kdBoundNode
	{
		~kdBoundNode()
		{
			if (!IsLeaf()) {
				if (data!=NULL) {
					delete data;
				}
			}
		}
		void createLeaf(const T *d, bound_t &boundBox)
		{
			flags = 3;
			nodeBound = boundBox;
			data = d;
		}
		void createInterior(int axis, PFLOAT d, bound_t &boundBox)
		{
			division = d;
			nodeBound = boundBox;
			flags = (flags & ~3) | axis;
			data = NULL;
		}
		PFLOAT 	SplitPos() const { return division; }
		int 	SplitAxis() const { return flags & 3; }
		int 	nPrimitives() const { return flags >> 2; }
		bool 	IsLeaf() const { return (flags & 3) == 3; }
		u_int32	getRightChild() const { return (flags >> 2); }
		void 	setRightChild(u_int32 i) { flags = (flags&3) | (i << 2); }
		
		PFLOAT division;
		const T *data;
		bound_t nodeBound;
		int		photonNum;
		u_int32	flags;
	};

template<class NodeData> struct ComparePhotonNode
{
	ComparePhotonNode(int a) { axis = a; }
	int axis;
	bool operator()(const NodeData *d1,	const NodeData *d2) const
	{
		return d1->pos[axis] == d2->pos[axis] ? (d1 < d2) : d1->pos[axis] < d2->pos[axis];
	}
};

template <class T>
class photonKdTree
{
	public:
		photonKdTree(const std::vector<T> &dat);
		~photonKdTree(){ if(nodes) y_free(nodes); }
		template<class LookupProc> void lookup(const point3d_t &p, const LookupProc &proc, PFLOAT &maxDistSquared) const;
		double lookupStat()const{ return double(Y_PROCS)/double(Y_LOOKUPS); } //!< ratio of photons tested per lookup call
		int GetPhotons( const point3d_t& p, std::vector<const T*> &array, float threshold );
	protected:
		template<class LookupProc> void recursiveLookup(const point3d_t &p, const LookupProc &proc, PFLOAT &maxDistSquared, int nodeNum) const;
		void recursiveGetPhotons( const point3d_t& p, std::vector<const T*> &array, int nodeNum, float threshold );
		void recursiveSumPhotons(int nodeNum);
		struct KdStack
		{
			const kdBoundNode<T> *node; //!< pointer to far child
			PFLOAT s; 		//!< the split val of parent node
			int axis; 		//!< the split axis of parent node
		};
		void buildTree(u_int32 start, u_int32 end, bound_t &nodeBound, const T **prims);
		void buildTree2(u_int32 start, u_int32 end, bound_t &nodeBound, const T **prims, int axis=0);
		kdBoundNode<T> *nodes;
		u_int32 nElements, nextFreeNode;
		bound_t treeBound;
		mutable unsigned int Y_LOOKUPS, Y_PROCS;
};


template<class T>
photonKdTree<T>::photonKdTree(const std::vector<T> &dat)
{
	Y_LOOKUPS=0; Y_PROCS=0;
	nextFreeNode = 0;
	nElements = dat.size();
	
	if(nElements == 0)
	{
		Y_ERROR << "photonKdTree: Empty vector!" << yendl;
		return;
	}
	
	nodes = (kdBoundNode<T> *)y_memalign(64, 4*nElements*sizeof(kdBoundNode<T>)); //actually we could allocate one less...2n-1
	const T **elements = new const T*[nElements];
	
	for(u_int32 i=0; i<nElements; ++i) elements[i] = &dat[i];
	
	treeBound.set(dat[0].pos, dat[0].pos);
	
	for(u_int32 i=1; i<nElements; ++i) treeBound.include(dat[i].pos);
	
	Y_INFO << "photonKdTree: Starting recusive tree build for "<<nElements<<" elements..." << yendl;
	
	buildTree(0, nElements, treeBound, elements);
	
	recursiveSumPhotons(0);
	
	
	std::cout << "test: " << nodes[0].photonNum << std::endl;
	Y_INFO << "photonKdTree: Tree built." << yendl;
	
	delete[] elements;
}

template<class T>
void photonKdTree<T>::buildTree(u_int32 start, u_int32 end, bound_t &nodeBound, const T **prims)
{
	if(end - start == 1)
	{
		nodes[nextFreeNode].createLeaf(prims[start],nodeBound);
		nodes[nextFreeNode].photonNum = 1;
		nextFreeNode++;
		return;
	}
	int splitAxis = nodeBound.largestAxis();
	int splitEl = (start+end)/2;
	std::nth_element(&prims[start], &prims[splitEl],
					&prims[end], ComparePhotonNode<T>(splitAxis));
	u_int32 curNode = nextFreeNode;
	PFLOAT splitPos = prims[splitEl]->pos[splitAxis];
	nodes[curNode].createInterior(splitAxis, splitPos,nodeBound);
	// add boundbox and summarry photon here
	nodes[curNode].photonNum = end-start;
	
	++nextFreeNode;
	bound_t boundL = nodeBound, boundR = nodeBound;
	switch(splitAxis){
		case 0: boundL.setMaxX(splitPos); boundR.setMinX(splitPos); break;
		case 1: boundL.setMaxY(splitPos); boundR.setMinY(splitPos); break;
		case 2: boundL.setMaxZ(splitPos); boundR.setMinZ(splitPos); break;
	}
	//<< recurse below child >>
	buildTree(start, splitEl, boundL, prims);
	//<< recurse above child >>
	nodes[curNode].setRightChild (nextFreeNode);
	buildTree(splitEl, end, boundR, prims);
}

template<class T>
void photonKdTree<T>::buildTree2(u_int32 start, u_int32 end, bound_t &nodeBound, const T **prims, int axis)
{
	if(end - start == 1)
	{
		nodes[nextFreeNode].createLeaf(prims[start], nodeBound);
		nodes[nextFreeNode].photonNum = 1;
		nextFreeNode++;
		return;
	}

	int splitAxis = axis;
	int splitEl = (start+end)/2;
	std::nth_element(&prims[start], &prims[splitEl],
					&prims[end], ComparePhotonNode<T>(splitAxis));
	u_int32 curNode = nextFreeNode;
	PFLOAT splitPos = prims[splitEl]->pos[splitAxis];
	nodes[curNode].createInterior(splitAxis, splitPos,nodeBound);
	nodes[curNode].photonNum = end-start;
	
	++nextFreeNode;
	bound_t boundL = nodeBound, boundR = nodeBound;
	switch(splitAxis){
		case 0: boundL.setMaxX(splitPos); boundR.setMinX(splitPos); break;
		case 1: boundL.setMaxY(splitPos); boundR.setMinY(splitPos); break;
		case 2: boundL.setMaxZ(splitPos); boundR.setMinZ(splitPos); break;
	}
	//<< recurse below child >>
	buildTree2(start, splitEl, boundL, prims, (axis+1)%3);
	//<< recurse above child >>
	nodes[curNode].setRightChild (nextFreeNode);
	buildTree2(splitEl, end, boundR, prims, (axis+1)%3);
}


template<class T> template<class LookupProc> 
void photonKdTree<T>::lookup(const point3d_t &p, const LookupProc &proc, PFLOAT &maxDistSquared) const
{
#if NON_REC_LOOKUP > 0
	++Y_LOOKUPS;
	KdStack stack[KD_MAX_STACK];
	const kdBoundNode<T> *farChild, *currNode = nodes;
	
	int stackPtr = 1;
	stack[stackPtr].node = 0; // "nowhere", termination flag
	
	while (true)
	{
		while( !currNode->IsLeaf() )
		{
			int axis = currNode->SplitAxis();
			PFLOAT splitVal = currNode->SplitPos();
			
			if( p[axis] < splitVal ) //need traverse left first
			{
				farChild = &nodes[currNode->getRightChild()];
				currNode++;
			}
			else //need traverse right child first
			{
				farChild = currNode+1;
				currNode = &nodes[currNode->getRightChild()];
			}
			++stackPtr;
			stack[stackPtr].node = farChild;
			stack[stackPtr].axis = axis;
			stack[stackPtr].s = splitVal;
		}

		// Hand leaf-data kd-tree to processing function
		vector3d_t v = currNode->data->pos - p;
		PFLOAT dist2 = v.lengthSqr();

		if (dist2 < maxDistSquared)
		{
			++Y_PROCS;
			proc(currNode->data, dist2, maxDistSquared);
		}
		
		if(!stack[stackPtr].node) return; // stack empty, done.
		//radius probably lowered so we may pop additional elements:
		int axis = stack[stackPtr].axis;
		dist2 = p[axis] - stack[stackPtr].s;
		dist2 *= dist2;

		while(dist2 > maxDistSquared)
		{
			--stackPtr;
			if(!stack[stackPtr].node) return;// stack empty, done.
			axis = stack[stackPtr].axis;
			dist2 = p[axis] - stack[stackPtr].s;
			dist2 *= dist2;
		}
		currNode = stack[stackPtr].node;
		--stackPtr;
	}
#else
	recursiveLookup(p, proc, maxDistSquared, 0);
	++Y_LOOKUPS;
	if(Y_LOOKUPS == 159999)
	{
		Y_INFO << "pointKd-Tree:average photons tested per lookup:" << double(Y_PROCS)/double(Y_LOOKUPS) << yendl;
	}
#endif
}

template<class T> template<class LookupProc> 
void photonKdTree<T>::recursiveLookup(const point3d_t &p, const LookupProc &proc, PFLOAT &maxDistSquared, int nodeNum) const
{
	const kdBoundNode<T> *currNode = &nodes[nodeNum];
	if(currNode->IsLeaf())
	{
		vector3d_t v = currNode->data->pos - p;
		PFLOAT dist2 = v.lengthSqr();
		if (dist2 < maxDistSquared)
			proc(currNode->data, dist2, maxDistSquared);
			++Y_PROCS;
		return;
	}
	int axis = currNode->SplitAxis();
	PFLOAT dist2 = p[axis] - currNode->SplitPos();
	dist2 *= dist2;
	if(p[axis] < currNode->SplitPos())
	{
		recursiveLookup(p, proc, maxDistSquared, nodeNum+1);
		if (dist2 < maxDistSquared)
			recursiveLookup(p, proc, maxDistSquared, currNode->getRightChild());
	}
	else
	{
		recursiveLookup(p, proc, maxDistSquared, currNode->getRightChild());
		if (dist2 < maxDistSquared)
			recursiveLookup(p, proc, maxDistSquared, nodeNum+1);
	}
}
	
	template<class T>
	int photonKdTree<T>::GetPhotons( const point3d_t& p, std::vector<const T*> &array, float threshold )
	{
		array.clear();
		recursiveGetPhotons(p,array,0,threshold);
		return array.size();
	}
	
	template<class T>
	void photonKdTree<T>::recursiveSumPhotons(int nodeNum)
	{
		kdBoundNode<T> *currNode = &nodes[nodeNum];
		if(currNode->IsLeaf())
			return;
		
		// compute left
		recursiveSumPhotons(nodeNum+1);
		
		// compute right
		recursiveSumPhotons(currNode->getRightChild());
		
		// compute current
		// pos
		T* dat = new T();
		float weight = (float)nodes[nodeNum+1].photonNum/(float)currNode->photonNum;
		dat->pos = nodes[nodeNum+1].data->pos*weight + nodes[currNode->getRightChild()].data->pos*(1-weight);
		
		// c
		dat->c = nodes[nodeNum+1].data->c + nodes[currNode->getRightChild()].data->c;
		
		// dir
		vector3d_t left(nodes[nodeNum+1].data->dir);
		vector3d_t right(nodes[currNode->getRightChild()].data->dir);
		vector3d_t newDir = left*weight + right*(1.f-weight);
		dat->dir = newDir.normalize();;//nodes[nodeNum+1].data->dir*weight + nodes[currNode->getRightChild()].data->dir*(1-weight);
		
		//normal
		dat->hitNormal = nodes[nodeNum+1].data->hitNormal*weight + nodes[currNode->getRightChild()].data->hitNormal*(1-weight);
		dat->hitNormal.normalize();
		currNode->data = dat;
	}
	
	template<class T>
	void photonKdTree<T>::recursiveGetPhotons( const point3d_t& p, std::vector<const T*> &array, int nodeNum, float threshold )
	{
		kdBoundNode<T> *currNode = &nodes[nodeNum];
		if(currNode->IsLeaf())
		{
			array.push_back(currNode->data);
			return;
		}
		
		if( !currNode->nodeBound.includes(p) )
		{
			// compute the angle of node to p
			vector3d_t ptoCenter = p - currNode->nodeBound.center();
			vector3d_t boxDiag = currNode->nodeBound.g-currNode->nodeBound.a;
			
			//std::cout << "ptoCenter is " << ptoCenter << std::endl;
			//std::cout << "boxDiag is " << boxDiag << std::endl;
			
			float diagDis = boxDiag.length();
			float ptocDis = ptoCenter.length();
			ptoCenter.normalize();
			boxDiag.normalize();
			
			float cosAng1 = ptoCenter*boxDiag;
			boxDiag.x *= -1.f;
			float cosAng2 = ptoCenter*boxDiag;
			
			float cosAng = fabs(cosAng1)<fabs(cosAng2)?cosAng1:cosAng2;
			float sinAng = fSqrt(1-cosAng*cosAng);
			
			float aperture = sinAng*diagDis/ptocDis;
			
			if (ptocDis > 0.1 && aperture <= threshold) {
				array.push_back(currNode->data);
				return;
			}
		}
		// compute left
		recursiveGetPhotons(p,array,nodeNum+1,threshold);
		
		// compute right
		recursiveGetPhotons(p,array,currNode->getRightChild(),threshold);
	}
} // namespace::kdtree

__END_YAFRAY

#endif // Y_PHOTONKDTREE_H
