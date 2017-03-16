#pragma once

#include <vector>

#include "boost\geometry\geometry.hpp"
#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "Vec3f.h"
#include "octree.h"


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
		class Particle {
			// A single particle of the SPH simulation
		public:
			Vec3f position, velocity, velocity_half, acceleration;
			float
				density, density_derivative,
				pressure,
				temperature, temperature_derivative,
				viscocity;
			boost::numeric::ublas::matrix<float> stress_tensor;

			Particle(void);
		};

		Particle  particles[num_of_particles];

		void randomize_particles();

		void calculate_initial_conditions();

		void simulation_step();
};