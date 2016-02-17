#ifndef BVHCollisionDetector_h
#define BVHCollisionDetector_h

#include <vector>
#include <pmmintrin.h>


template<class Real, int Dims, class Vector, class SpatialObject, class Config>
class BVHCollisionDetector
{
private:
	
	// Hehe.. compile time goodness!! >=)
	template<typename T, typename U> struct SameType { enum { value = false }; };
	template<typename T> struct SameType<T, T> { enum { value = true }; };
	template<int N, int M> struct Equals { enum { value = (N == M) }; };
	template<int N, int M> struct And { enum { value = ((N) && (M)) }; };
	enum { ShouldUseSSE = And< Equals<Dims, 3>::value , SameType<Real, float>::value >::value };
	enum { BoundCapacity = ShouldUseSSE ? 4 : Dims };
	
	struct Bound
	{
		Real coors[BoundCapacity];
		
		inline Bound(const Vector &v)
		{
			Config c;
			for (int d = 0; d < Dims; ++d) coors[d] = c.getAxis(v, d);
		}
		
		inline Bound() {}
		
		inline Real operator[](size_t i) const
		{
			return coors[i];
		}
		
		inline Real & operator [](size_t i)
		{
			return coors[i];
		}
		
		inline bool anyLessThan(const Bound &b) const
		{
			for (int d = 0; d < Dims; ++d) if (coors[d] < b[d]) return true;
			return false;
		}
		
		inline Bound minWith(const Bound &b) const
		{
			Bound result;
			for (int d = 0; d < Dims; ++d) result[d] = coors[d] < b[d] ? coors[d] : b[d];
			return result;
		}
		
		inline Bound maxWith(const Bound &b) const
		{
			Bound result;
			for (int d = 0; d < Dims; ++d) result[d] = coors[d] > b[d] ? coors[d] : b[d];
			return result;
		}
	};
	
	struct Box
	{
		SpatialObject *object;
		Bound lowerBound;
		Bound upperBound;
		Bound center;
		
		inline Box() {}
		inline Box(SpatialObject *object): object(object) {}
		
	};
	
	struct BVHNode
	{
		size_t leftIndex;
		size_t rightIndex;
		Bound lowerBound;
		Bound upperBound;
	};
	
	std::vector<Box> objects;
	std::vector<size_t> objectIndices;
	std::vector<BVHNode> bvhNodes;
	size_t bvhRootIndex;
	
	// Temp vars for intersection testing
	std::vector<int> collided;
	std::vector<size_t> intersected;
	Bound currentBoxLowerBound;
	Bound currentBoxUpperBound;
	
	template <class T> inline void swap(T &a, T &b) { T t = a; a = b; b = t; }
	
	template <int Axis> inline void quickSelect(size_t select, size_t begin, size_t end)
	{
		if (end - 1 == select) return;
		size_t pi = ((end - begin) >> 1) + begin;
		pi -= pi > 0;
		Real pv = objects[objectIndices[pi]].center[Axis];
		size_t front = begin;
		size_t back = end - 1;
		while (front < back) {
			if (objects[objectIndices[front]].center[Axis] < pv) ++front;
			else if (objects[objectIndices[back]].center[Axis] > pv) --back;
			else swap(objectIndices[front++], objectIndices[back--]);
		}
		front += (front == back && objects[objectIndices[front]].center[Axis] <= pv);
		if (select < front) quickSelect<Axis>(select, begin, front);
		else quickSelect<Axis>(select, front, end);
	}
	
	// Helps improve speed marginally... the code is mostly (~90%) cache bound.
	inline void findPossibleIntersectedSSE(size_t bvhNodeIndex)
	{
		BVHNode &node = bvhNodes[bvhNodeIndex];
		
		if (node.rightIndex) {
			
			__m128 boxLB = _mm_loadu_ps((float*)&currentBoxLowerBound);
			__m128 boxUB = _mm_loadu_ps((float*)&currentBoxUpperBound);
			
			Bound &bvhUpperBound0 = bvhNodes[node.leftIndex].upperBound;
			Bound &bvhLowerBound0 = bvhNodes[node.leftIndex].lowerBound;
			Bound &bvhUpperBound1 = bvhNodes[node.rightIndex].upperBound;
			Bound &bvhLowerBound1 = bvhNodes[node.rightIndex].lowerBound;
			
			__m128 bvhUB0 = _mm_loadu_ps((float*)&bvhUpperBound0);
			__m128 bvhLB0 = _mm_loadu_ps((float*)&bvhLowerBound0);
			__m128 bvhUB1 = _mm_loadu_ps((float*)&bvhUpperBound1);
			__m128 bvhLB1 = _mm_loadu_ps((float*)&bvhLowerBound1);
			
			__m128 cmp0 = _mm_and_ps(_mm_cmpge_ps(bvhUB0, boxLB), _mm_cmpge_ps(boxUB, bvhLB0));
			__m128 cmp1 = _mm_and_ps(_mm_cmpge_ps(bvhUB1, boxLB), _mm_cmpge_ps(boxUB, bvhLB1));
			
			__m128 cmp = _mm_and_ps(_mm_and_ps(_mm_shuffle_ps(cmp0, cmp1, (0<<4)|0),
											   _mm_shuffle_ps(cmp0, cmp1, (1<<4)|1)),
									_mm_shuffle_ps(cmp0, cmp1, (2<<4)|2));
			
			if (((int*)&cmp)[0]) findPossibleIntersectedSSE(node.leftIndex);
			if (((int*)&cmp)[2]) findPossibleIntersectedSSE(node.rightIndex);
			
		} else {
			
			intersected.push_back(node.leftIndex);
		}
		
	}
	
	inline void findPossibleIntersected(size_t bvhNodeIndex)
	{
		BVHNode &node = bvhNodes[bvhNodeIndex];
		if (!(node.upperBound.anyLessThan(currentBoxLowerBound) ||
			  currentBoxUpperBound.anyLessThan(node.lowerBound))) {
			if (node.rightIndex) {
				findPossibleIntersected(node.leftIndex);
				findPossibleIntersected(node.rightIndex);
			} else {
				intersected.push_back(node.leftIndex);
			}
		}
	}
	
	inline void findPossibleIntersected(size_t bvhNodeIndex, const Bound &boxLowerBound, const Bound &boxUpperBound)
	{
		intersected.clear();
		
		currentBoxLowerBound = boxLowerBound;
		currentBoxUpperBound = boxUpperBound;
		
		if (ShouldUseSSE) findPossibleIntersectedSSE(bvhNodeIndex);
		else findPossibleIntersected(bvhNodeIndex);
		
	}
	
	
	inline BVHNode makeBVH(size_t objectIndex)
	{
		BVHNode node;
		node.leftIndex = objectIndex;
		node.rightIndex = 0;
		node.lowerBound = objects[objectIndex].lowerBound;
		node.upperBound = objects[objectIndex].upperBound;
		return node;
	}
	
	inline BVHNode makeBVH(size_t leftIndex, size_t rightIndex)
	{
		BVHNode node;
		node.leftIndex = leftIndex;
		node.rightIndex = rightIndex;
		node.lowerBound = bvhNodes[leftIndex].lowerBound.minWith(bvhNodes[rightIndex].lowerBound);
		node.upperBound = bvhNodes[leftIndex].upperBound.maxWith(bvhNodes[rightIndex].upperBound);
		return node;
	}
	
	template <int Axis> inline size_t makeBVH(size_t begin, size_t end)
	{
		size_t len = end - begin;
		if (len == 1) {
			bvhNodes.push_back(makeBVH(objectIndices[begin]));
			return bvhNodes.size() - 1;
		}
		
		size_t mid = (begin + end) >> 1;
		quickSelect<Axis>(mid, begin, end);
		
		bvhNodes.push_back(makeBVH(makeBVH<(Axis + 1) % Dims>(begin, mid),
								   makeBVH<(Axis + 1) % Dims>(mid, end)));
		return bvhNodes.size() - 1;
	}
	
	inline void createBVH()
	{
		Config c;
		for (size_t i = 0; i < objects.size(); ++i) {
			objects[i].lowerBound = Bound(c.getLowerBound(*objects[i].object));
			objects[i].upperBound = Bound(c.getUpperBound(*objects[i].object));
			objects[i].center = Bound(c.getPosition(*objects[i].object));
		}
		bvhNodes.clear();
		bvhRootIndex = makeBVH<0>(0, objects.size());
	}
	
	
public:
	void addObjects(const std::vector<SpatialObject*> &spatialObjects)
	{
		size_t objectsPrevSize = objects.size();
		size_t objectsSize = objects.size() + spatialObjects.size();
		objects.resize(objectsSize);
		objectIndices.resize(objectsSize);
		for (size_t i = objectsPrevSize, si = 0; i < objectsSize; ++i, ++si) {
			objects[i] = Box(spatialObjects[si]);
			objectIndices[i] = i;
		}
	}
	
	void setObjects(const std::vector<SpatialObject*> &spatialObjects)
	{
		objects.resize(spatialObjects.size());
		objectIndices.resize(spatialObjects.size());
		for (size_t i = 0; i < objects.size(); ++i) {
			objects[i] = Box(spatialObjects[i]);
			objectIndices[i] = i;
		}
	}
	
	inline std::vector<SpatialObject*> getCollidedObjects()
	{
		Config c;
		
		createBVH();
		
		collided.resize(objects.size());
		for (size_t i = 0; i < objectIndices.size(); ++i) collided[i] = 0;
		
		for (size_t i = 0; i < objectIndices.size(); ++i) {
			findPossibleIntersected(bvhRootIndex, objects[objectIndices[i]].lowerBound, objects[objectIndices[i]].upperBound);
			for (size_t j = 0; j < intersected.size(); ++j) {
				if (objects[intersected[j]].object != objects[objectIndices[i]].object) {
					if (c.checkCollision(*objects[intersected[j]].object, *objects[objectIndices[i]].object)) {
						collided[intersected[j]] = 1;
						collided[objectIndices[i]] = 1;
						
					}
				}
			}
		}
		
		size_t numCollided = 0;
		std::vector<SpatialObject*> collidedObjects(collided.size());
		for (size_t i = 0; i < collided.size(); ++i) {
			if (collided[i]) {
				collidedObjects[numCollided++] = objects[i].object;
			}
		}
		
		collidedObjects.resize(numCollided);
		
		return collidedObjects;
	}
};


#endif