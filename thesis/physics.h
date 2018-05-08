#pragma once

#include <vector>

#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "particle.h"
#include "Vec3f.h"
#include "searchgrid.h"


class CubicSpline {
	private:
		const float h, a;

	public:
		CubicSpline(const float h);

		float operator()(const Vec3f &r) const;
		Vec3f gradient(const Vec3f &r) const;
};

class ParticleSystem {
	// Holds the particles and handles the physics simulation

	private:
		float simulation_time;
		int num_of_fluid_particles;
		ParticleContainer
			particles, prev_particles, next_particles;
		std::unique_ptr<Vec3f[]> acceleration;
		std::unique_ptr<float[]> density_derivative;
		const CaseDef case_def;
		SearchGrid search_grid_fluid, search_grid_boundary;
		CubicSpline cubic_spline;
		int verlet_step;

		// Generate particles for the geomtery specified in the case_def member
		void generate_particles();

		// Calculate a time step that is stable.
		float calculate_time_step() const;

		void compute_derivatives();

		// Integrate forward using verlet
		void integrate_verlet(float dt);

		SearchGrid::cell_indices_container get_all_neighbors(const Vec3f &position) const;

	public:
		ParticleSystem(const CaseDef &case_def) :
			simulation_time(0.0f),
			case_def(case_def),
			search_grid_fluid(case_def.particles.point_min, case_def.particles.point_max, case_def.h),
			search_grid_boundary(case_def.particles.point_min, case_def.particles.point_max, case_def.h),
			cubic_spline(CubicSpline(case_def.h)),
			verlet_step(0)
		{
			generate_particles();

			// We need to ensure that all particle arrays are the same size. We also
			// want the boundary positions to be set here so we don't have to copy
			// them every time. Simply copying the whole vector is fast enough.
			prev_particles = particles;
			next_particles = particles;

			// Note: these need different sizes because acceleration is not needed for boundaries
			acceleration = std::make_unique<Vec3f[]>(num_of_fluid_particles);
			density_derivative = std::make_unique<float[]>(particles.size());
		}

		// Delete these to make sure ParticleSystem is only ever passed by reference.
		ParticleSystem(const ParticleSystem &other) = delete;
		ParticleSystem& operator=(const ParticleSystem &other) = delete;

		const ParticleContainer& get_particlearray() const { return particles; }

		float current_time() const {return simulation_time;}

		auto get_num_of_particles() const { return particles.size(); }
		auto get_num_of_fluid_particles() const { return num_of_fluid_particles; }

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};