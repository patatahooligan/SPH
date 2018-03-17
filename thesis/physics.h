#pragma once

#include <vector>

#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "particle.h"
#include "Vec3f.h"
#include "kdtree.h"

using kernel_function_t = float(const Vec3f&, const float);
using kernel_function_derivative_t = Vec3f(const Vec3f&, const float);

float piecewise_smoothing_kernel(const Vec3f &r, const float h);
Vec3f piecewise_smoothing_kernel_derivative(const Vec3f &r, const float h);

class ParticleSystem {
	// Holds the particles and handles the physics simulation

	private:
		float simulation_time;
		particlearray particles;
		ParticleAdaptor kd_tree_adaptor;
		ParticleKDTree kd_tree;
		CaseDef case_def;

		kernel_function_t &smoothing_kernel;
		kernel_function_derivative_t &smoothing_kernel_derivative;

		// Calculate a time step that is stable.
		float calculate_time_step();

		void update_derivatives();

		// Integrate forward using leap-frog
		void integrate_step();

		// Handle particle-wall collision
		void conflict_resolution();

		Vec3f boundary_force(const Particle& p);


	public:
		ParticleSystem(
			const CaseDef &case_def,
			kernel_function_t& smoothing_kernel = piecewise_smoothing_kernel,
			kernel_function_derivative_t& smoothing_kernel_derivative = piecewise_smoothing_kernel_derivative
		) :
			simulation_time(0.0f),
			kd_tree_adaptor(particles),
			kd_tree(3, kd_tree_adaptor),
			case_def(case_def),
			smoothing_kernel(smoothing_kernel),
			smoothing_kernel_derivative(piecewise_smoothing_kernel_derivative)
		{}

		// Delete these to make sure ParticleSystem is only ever passed by reference.
		ParticleSystem(const ParticleSystem &other) = delete;
		ParticleSystem& operator=(const ParticleSystem &other) = delete;

		const particlearray& get_particlearray() const { return particles; }

		float current_time() const {return simulation_time;}

		// Randomly insert particles in the bounding box defined by const float size.
		void randomize_particles();

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};