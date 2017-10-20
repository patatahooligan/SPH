#pragma once

#include <vector>

#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "particle.h"
#include "Vec3f.h"
#include "kdtree.h"

class Octree;

class ParticleSystem {
	// Holds the particles and handles the physics simulation

	private:
		float simulation_time;
		ParticleAdaptor kd_tree_adaptor;
		ParticleKDTree kd_tree;

		float smoothing_kernel(const Vec3f &r, const float h);
		Vec3f smoothing_kernel_derivative(const Vec3f &r, const float h);

		// Calculate a time step that is stable.
		float calculate_time_step();

		void update_derivatives();

		// Integrate forward using leap-frog
		void integrate_step();

		// Handle particle-wall collision
		void conflict_resolution();


	public:		
		std::array<Particle, num_of_particles> particles;

		ParticleSystem() :
			simulation_time(0.0f),
			kd_tree_adaptor(particles),
			kd_tree(3, kd_tree_adaptor) {}
		~ParticleSystem() {}

		// Delete these to make sure ParticleSystem is only ever passed by reference.
		ParticleSystem(const ParticleSystem &other) = delete;
		ParticleSystem& operator=(const ParticleSystem &other) = delete;

		float current_time() {return simulation_time;}

		// Randomly insert particles in the bounding box defined by const float size.
		void randomize_particles();

		// Calculates the initial conditions of the system. Must be used after particles are created
		// and before the first simulation_step() call.
		void calculate_initial_conditions();

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};