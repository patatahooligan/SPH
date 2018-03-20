#pragma once

#include "stdafx.h"

#include "nanoflann.hpp"

#include "constants.h"
#include "particle.h"

class ParticleAdaptor {
	// This class provides an interface for nanoflann

	private:
		const ParticleContainer &particles;

	public:
		ParticleAdaptor(const ParticleContainer &target_array) :
			particles(target_array)
			{};

		inline size_t kdtree_get_point_count() const { return num_of_particles; }

		// This function should never be called with an invalid dim
		// Because it is used in a heavy computational loop, we have to omit
		// any error checking and rely on proper calling
		#pragma warning(suppress: 4715)
		inline float kdtree_get_pt(const size_t idx, int dim) const {
			assert(dim >= 0 && dim <= 2);
			switch (dim) {
				case 0:
					return particles[idx].position.x;
				case 1:
					return particles[idx].position.y;
				case 2:
					return particles[idx].position.z;
			}
		}

		template <class BBOX>
		inline bool kdtree_get_bbox(BBOX &bb) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, ParticleAdaptor>,
	ParticleAdaptor,
	3, size_t>
	ParticleKDTree;