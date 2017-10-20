#pragma once

#include "stdafx.h"

#include "nanoflann.hpp"

#include "constants.h"
#include "particle.h"

class ParticleAdaptor {
	// This class provides an interface for nanoflann

	private:
		// Just to make the following declaration more readable
		using particlearray = std::array<Particle, num_of_particles>;
		const particlearray &particles;

	public:
		ParticleAdaptor() = delete;		// Must construct with reference to particles
		ParticleAdaptor(const particlearray &target_array):
			particles(target_array) {};

		inline size_t kdtree_get_point_count() const { return num_of_particles; }

		inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t size) const {
			float
				dx = p1[0] - particles[idx_p2].position.x,
				dy = p1[1] - particles[idx_p2].position.y,
				dz = p1[2] - particles[idx_p2].position.z;
			return dx*dx + dy*dy + dz*dz;
		}

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

		template <class BBOX>
		inline bool kdtree_get_bbox(BBOX &bb) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, ParticleAdaptor>,
	ParticleAdaptor,
	3, size_t>
	ParticleKDTree;