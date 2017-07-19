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

		// Calculate a time step that is stable.
		float calculate_time_step();

		void update_derivatives();

		// Integrate forward using leap-frog
		void integrate_step();

		// Handle particle-wall collision
		void conflict_resolution();

		Octree search_tree;

	public:		
		Particle  particles[num_of_particles];

		ParticleSystem() : simulation_time(0.0f) {}
		~ParticleSystem() {}

		float current_time() {return simulation_time;}

		// Randomly insert particles in the bounding box defined by const float size.
		void randomize_particles();

		// Calculates the initial conditions of the system. Must be used after particles are created
		// and before the first simulation_step() call.
		void calculate_initial_conditions();

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};