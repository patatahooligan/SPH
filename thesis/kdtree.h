#pragma once

#include "stdafx.h"

#include "nanoflann.hpp"

#include "constants.h"
#include "particle.h"

class ParticleAdaptor {
	// This class provides an interface for nanoflann

	private:
		const particlearray &particles;

	public:
		ParticleAdaptor() = delete;		// Must construct with reference to particles
		ParticleAdaptor(const particlearray &target_array) :
			particles(target_array)
			{};

		inline size_t kdtree_get_point_count() const { return num_of_particles; }

		inline float kdtree_get_pt(const size_t idx, int dim) const {
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