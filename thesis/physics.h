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
		ParticleContainer
			fluid_particles, boundary_particles,
			prev_fluid_particles, prev_boundary_particles,
			next_fluid_particles, next_boundary_particles;
		std::unique_ptr<Vec3f[]> acceleration;
		std::unique_ptr<float[]> density_derivative;
		ParticleAdaptor kd_tree_adaptor;
		ParticleKDTree kd_tree;
		const CaseDef case_def;
		int verlet_step;

		kernel_function_t &smoothing_kernel;
		kernel_function_derivative_t &smoothing_kernel_derivative;

		// Generate particles for the geomtery specified in the case_def member
		void generate_particles();

		// Calculate a time step that is stable.
		float calculate_time_step();

		void update_derivatives();

		// Integrate forward using verlet
		void integrate_verlet(float dt);

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
			kd_tree_adaptor(fluid_particles),
			kd_tree(3, kd_tree_adaptor),
			case_def(case_def),
			smoothing_kernel(smoothing_kernel),
			smoothing_kernel_derivative(piecewise_smoothing_kernel_derivative),
			verlet_step(0)
		{
			generate_particles();

			prev_fluid_particles.resize(fluid_particles.size());
			prev_boundary_particles.resize(boundary_particles.size());
			next_fluid_particles.resize(fluid_particles.size());
			next_boundary_particles.resize(boundary_particles.size());
			acceleration = std::make_unique<Vec3f[]>(fluid_particles.size());
			density_derivative = std::make_unique<float[]>(fluid_particles.size());
		}

		// Delete these to make sure ParticleSystem is only ever passed by reference.
		ParticleSystem(const ParticleSystem &other) = delete;
		ParticleSystem& operator=(const ParticleSystem &other) = delete;

		const ParticleContainer& get_particlearray() const { return fluid_particles; }

		float current_time() const {return simulation_time;}

		auto num_of_particles() const { return fluid_particles.size(); }

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};