#pragma once

#include "stdafx.h"

#include "nanoflann.hpp"

#include "constants.h"
#include "particle.h"

class ParticleAdaptor {
	// This class provides an interface for nanoflann

	private:
		// Just to make the following declaration more readable
		typedef Particle particlearray[num_of_particles];
		particlearray &particles;

	public:
		ParticleAdaptor() = delete;		// Must construct with reference to particles
		ParticleAdaptor(particlearray &target_array):
			particles(target_array) {};

		inline size_t kdtree_get_point_count() const { return num_of_particles; }

		inline float kdtree_get_pt(const size_t idx, int dim) const {
			// TODO: after testing this, consider removing the invalid argument check
			switch (dim) {
				case 0:
					return particles[idx].position.x;
				case 1:
					return particles[idx].position.y;
				case 2:
					return particles[idx].position.z;
				default:
					throw std::invalid_argument(std::string("Invalid dimension request ") + std::to_string(dim));
			}
		}
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::metric_L2, ParticleAdaptor, 3, size_t> ParticleKDTree;