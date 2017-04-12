#pragma once

#include <vector>

#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "particle.h"
#include "Vec3f.h"
#include "octree.h"

class Octree;

class ParticleSystem {
	// Holds the particles and handles the physics simulation

	private:
		float simulation_time;
		float smoothing_kernel(const Vec3f &r, const float h);
		Vec3f smoothing_kernel_derivative(const Vec3f &r, const float h);

		float calculate_time_step();
		void update_derivatives();
		void integrate_step();
		void conflict_resolution();

		Octree search_tree;

	public:
		Particle  particles[num_of_particles];

		ParticleSystem() : simulation_time(0.0f) {}
		~ParticleSystem(){}

		void randomize_particles();

		void calculate_initial_conditions();

		void simulation_step();
};