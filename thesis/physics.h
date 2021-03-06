#pragma once

#include "../common/constants.h"
#include "../common/particle.h"
#include "../common/../common/vec3f.h"
#include "searchgrid.h"
#include "massspringdamper.h"


class CubicSpline {
	private:
		const float h, a;

	public:
		CubicSpline(float h);

		float operator()(float length) const;
		float operator()(const Vec3f &r) const;

		float gradient_coef(float length) const;
		Vec3f gradient(const Vec3f &r) const;
};

class CubicSplinePrecalculated {
	private:
		std::vector<float> values, gradient_values;
		float step;

	public:
		CubicSplinePrecalculated(float h, int resolution);

		float operator()(float length) const;
		float operator()(const Vec3f &r) const;
		Vec3f gradient(const Vec3f & r, float length) const;
		Vec3f gradient(const Vec3f &r) const;
};

class ParticleSystem {
	// Holds the particles and handles the physics simulation

	private:
		CaseDef case_def;
		ParticleContainer
			particles, prev_particles;
		std::vector<Vec3f> acceleration;
		std::vector<float> density_derivative, pressure;
		int num_of_fluid_particles;
		CaseDef::Box bounding_box;
		SearchGrid search_grid_fluid, search_grid_boundary;
		CubicSplinePrecalculated cubic_spline;
		float simulation_time = 0.0f;
		int verlet_step = 0;

		std::vector<MassSpringSystem> fluid_spring_systems, boundary_spring_systems;

		void allocate_memory_for_verlet_variables() {
			// Note: these need different sizes because acceleration is not needed for boundaries
			acceleration.resize(num_of_fluid_particles);
			density_derivative.resize(particles.size());
			pressure.resize(particles.size());
		}

		// Generate particles for the geomtery specified in the case_def member
		ParticleContainer generate_particles();

		void generate_mass_spring_damper();

		// Calculate a time step that is stable.
		float calculate_time_step() const;

		template <ParticleType TypeOfPj>
		void spring_forces();

		template <ParticleType TypeOfPi, ParticleType TypeOfNeighbors>
		void compute_derivatives(const int i);

		void compute_derivatives();

		// Integrate forward using verlet
		void integrate_verlet(float dt);

		SearchGrid::cell_indices_container get_fluid_neighbors(const Vec3f &position) const;
		SearchGrid::cell_indices_container get_boundary_neighbors(const Vec3f &position) const;
		SearchGrid::cell_indices_container get_all_neighbors(const Vec3f &position) const;

	public:
		using SpringIterator = std::vector<MassSpringDamper>::iterator;
		using SpringConstIterator = std::vector<MassSpringDamper>::const_iterator;

		struct State {
			ParticleContainer
				prev_fluid_particles, prev_boundary_particles,
				fluid_particles, boundary_particles;
			std::vector<MassSpringSystem> fluid_spring_systems, boundary_spring_systems;
			float simulation_time;
			int verlet_step;
		};

		ParticleSystem(const CaseDef &case_def);
		ParticleSystem(const CaseDef &case_def, State state);

		// Delete these to make sure ParticleSystem is only ever passed by reference.
		ParticleSystem(ParticleSystem &&) = default;
		ParticleSystem& operator=(const ParticleSystem &other) = delete;

		State get_current_state() const;

		const ParticleContainer& get_particlearray() const { return particles; }
		const ParticleContainer& get_previous_particlearray() const { return prev_particles; }

		ParticleConstIterator get_fluid_begin() const { return particles.cbegin(); }
		ParticleConstIterator get_fluid_end() const { return particles.cbegin() + num_of_fluid_particles; }
		ParticleConstIterator get_boundary_begin() const { return particles.cbegin() + num_of_fluid_particles; }
		ParticleConstIterator get_boundary_end() const { return particles.cend(); }

		//MassSpringConstIterator get_fluid_fluid_springs_begin() const { return mass_spring_damper.cbegin(); }
		//MassSpringConstIterator get_fluid_fluid_springs_end() const { return mass_spring_damper.cbegin() + num_of_fluid_fluid_springs; }

		//MassSpringConstIterator get_fluid_boundary_springs_begin() const { return mass_spring_damper.cbegin() + num_of_fluid_fluid_springs; }
		//MassSpringConstIterator get_fluid_boundary_springs_end() const { return mass_spring_damper.cend(); }

		auto current_time() const {return simulation_time;}
		auto current_step() const { return verlet_step; }

		auto get_num_of_particles() const { return particles.size(); }
		auto get_num_of_fluid_particles() const { return num_of_fluid_particles; }

		// Update all derivatives and integrate a single step forward.
		void simulation_step();
};