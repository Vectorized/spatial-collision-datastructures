#ifndef BVHSpatialStructure_h
#define BVHSpatialStructure_h

#include "ISpatialStructure.h"
#include "BVHCollisionDetector.h"

namespace SpatialTest
{
	class BVHSpatialStructure : public ISpatialStructure
	{
		/**
		 A class to configure the BVH Collision Detector
		 */
		struct Config
		{
			/**
			 Returns the lower bound [x, y, z] of the spatial object
			 */
			inline Vector3 getLowerBound(const ISpatialObject &o)
			{
				float r = o.VGetRadius();
				return o.VGetPosition() - Vector3(r, r, r);
			}
			
			/**
			 Returns the upper bound [x, y, z] of the spatial object
			 */
			inline Vector3 getUpperBound(const ISpatialObject &o)
			{
				float r = o.VGetRadius();
				return o.VGetPosition() + Vector3(r, r, r);
			}
			
			/**
			 Returns the vector's component for the axis
			 */
			inline float getAxis(const Vector3 &v, const int axis)
			{
				if (axis == 2) return v.z;
				if (axis == 1) return v.y;
				return v.x;
			}
			
			/**
			 Returns the position of a spatial object
			 */
			inline Vector3 getPosition(const ISpatialObject &o) { return o.VGetPosition(); }
			
			/**
			 Returns true if the two spatial objects collide
			 */
			inline bool checkCollision(const ISpatialObject &a, const ISpatialObject &b) { return a.VCheckCollision(&b); }
		};
		
		BVHCollisionDetector<float, 3, Vector3, ISpatialObject, Config> collisionDetector;
		
		
	public:
		
		void VAddObjects(const std::vector<ISpatialObject*>& refObjects)
		{
			collisionDetector.addObjects(refObjects);
		}
		
		
		void VUpdate()
		{
			std::vector<ISpatialObject*> collidedObjects = collisionDetector.getCollidedObjects();
			
			for (size_t i = 0; i < collidedObjects.size(); ++i) {
				collidedObjects[i]->VCollisionOn();
			}
		}
		
	protected:
		
		
	};
};


#endif