#pragma once

#include <vector>

#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "particle.h"
#include "Vec3f.h"
#include "searchgrid.h"
#include "massspringdamper.h"


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
		CaseDef case_def;
		ParticleContainer
			particles, prev_particles, next_particles;
		std::vector<Vec3f> acceleration;
		std::vector<float> density_derivative, pressure;
		int num_of_fluid_particles;
		CaseDef::Box bounding_box;
		SearchGrid search_grid_fluid, search_grid_boundary;
		CubicSpline cubic_spline;
		std::vector<MassSpringDamper> mass_spring_damper;
		float simulation_time = 0.0f;
		int verlet_step = 0;

		void allocate_memory_for_verlet_variables() {
			// Note: these need different sizes because acceleration is not needed for boundaries
			acceleration.resize(num_of_fluid_particles);
			density_derivative.resize(particles.size());
			pressure.resize(particles.size());
		}

		// Generate particles for the geomtery specified in the case_def member
		ParticleContainer generate_particles();

		void generate_mass_spring_damper();

		CaseDef::Box get_particle_axis_aligned_bounding_box() const;

		// Calculate a time step that is stable.
		float calculate_time_step() const;

		template <ParticleType TypeOfPi, ParticleType TypeOfNeighbors>
		void compute_derivatives(const int i);

		void compute_derivatives();

		// Integrate forward using verlet
		void integrate_verlet(float dt);

		void remove_out_of_bounds_particles();

		SearchGrid::cell_indices_container get_all_neighbors(const Vec3f &position) const;

	public:
		ParticleSystem(const CaseDef &case_def);

		// Delete these to make sure ParticleSystem is only ever passed by reference.
		ParticleSystem(ParticleSystem &&) = default;
		ParticleSystem& operator=(const ParticleSystem &other) = delete;

		const ParticleContainer& get_particlearray() const { return particles; }
		const ParticleContainer& get_previous_particlearray() const { return prev_particles; }

		ParticleConstIterator get_fluid_begin() const { return particles.begin(); }
		ParticleConstIterator get_fluid_end() const { return particles.begin() + num_of_fluid_particles; }
		ParticleConstIterator get_boundary_begin() const { return particles.begin() + num_of_fluid_particles; }
		ParticleConstIterator get_boundary_end() const { return particles.end(); }

		auto current_time() const {return simulation_time;}
		auto current_step() const { return verlet_step; }

		auto get_num_of_particles() const { return particles.size(); }
		auto get_num_of_fluid_particles() const { return num_of_fluid_particles; }

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};