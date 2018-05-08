#include "stdafx.h"

#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <cmath>
#include <cfloat>
#include <vector>

#include "physics.h"
#include "constants.h"
#include "vec3f.h"

CubicSpline::CubicSpline(const float h):
	h(h), a(1.0f / (pi * std::pow(h, 3))) {}

float CubicSpline::operator()(const Vec3f &r) const {
	const float
		q = r.length() / h;

	if (q < 1.0f)
		return a * (1 - (3.0f / 2.0f) * (q*q) + (3.0f / 4.0f) * std::pow(q, 3));
	else if (q < 2.0f)
		return a * (std::pow(2 - q, 3) / 4.0f);
	else
		return 0.0f;
}

Vec3f CubicSpline::gradient(const Vec3f &r) const {
	const float
		q = r.length() / h;

	// It is easier and faster to handle these edge cases separately
	if (q == 0.0f || q > 2.0f)
		return { 0.0f, 0.0f, 0.0f };

	const float
		dq = 1 / (h * h * q),
		dw = [&]() {
		if (q < 1.0f)
			return (-3.0f * q + 9.0f / 4.0f * q * q) * a;
		else if (q < 2.0f)
			return (3.0f * a * std::pow(1.0f - q, 2)) / 4.0f;
		else
			return 0.0f;
	}();

	return dw * dq * r;
}


void ParticleSystem::generate_particles() {
	auto fillbox = [&](const CaseDef::Box &box) {
		const float& density = case_def.particles.density;
		int
			x_increments = int(box.size.x / density) + 1,
			y_increments = int(box.size.y / density) + 1,
			z_increments = int(box.size.z / density) + 1;

		size_t num_of_particles = x_increments + y_increments + z_increments;

		particles.reserve(particles.size() + num_of_particles);

		for (int x = 0; x < x_increments; ++x) {
			for (int y = 0; y < y_increments; ++y) {
				for (int z = 0; z < z_increments; ++z) {
					Particle p;
					p.position = box.origin + Vec3f{ x * density, y * density, z * density };
					p.density = case_def.rhop0;
					particles.emplace_back(p);
				}
			}
		}
	};

	for (const auto& box : case_def.fluid_boxes)
		fillbox(box);

	// Storing fluid and bound particles in the same container greatly simplifies the iterative
	// algorithms. For the parts that apply specifically to one type of particle, we need to keep
	// track of the number of fluid particles and iterate accordingly.
	num_of_fluid_particles = particles.size();

	for (const auto& box : case_def.bound_boxes)
		fillbox(box);
}

float ParticleSystem::calculate_time_step() const {
	const float
		&h = case_def.h,
		&c = case_def.speedsound;

	float
		acceleration2_max = 0.0f,
		dt_cv = std::numeric_limits<float>::max();
	
	#pragma omp parallel for
	for (int i = 0; i < num_of_fluid_particles; ++i) {
		const Particle& Pi = particles[i];

		const auto index_ranges = get_all_neighbors(Pi.position);

		float sigma = 0.0f;
		for (const auto index_pair : index_ranges) {
			for (int j = index_pair.first; j < index_pair.second; ++j) {
				const Particle& Pj = particles[j];
				const Vec3f
					v_ij = Pj.velocity - Pi.velocity,
					r_ij = Pj.position - Pi.position;

				const float
					r2 = r_ij.length_squared();

				sigma = std::max(
					sigma,
					std::abs((case_def.h * dot_product(v_ij, r_ij)) / (r2 + 0.01f * h * h)));
			}
		}

		acceleration2_max = std::max(acceleration2_max, acceleration[i].length_squared());

		dt_cv = std::min(
			dt_cv,
			h / (c + sigma));
	}

	const float dt_f = std::sqrt(h / std::sqrt(acceleration2_max * case_def.particles.mass));

	return case_def.cflnumber * std::min(dt_cv, dt_f);
}

void ParticleSystem::compute_derivatives() {
	#pragma omp parallel for
	for (int i = 0; size_t(i) < particles.size(); ++i) {
		const Particle& Pi = particles[i];

		// Get neighbors of Pi
		const auto index_ranges = get_all_neighbors(Pi.position);

		constexpr int gamma = 7;
		const float
			beta = (case_def.speedsound * case_def.speedsound * case_def.rhop0) / gamma,
			Pi_pressure = beta * (std::pow(Pi.density / case_def.rhop0, gamma) - 1);

		if (i < num_of_fluid_particles)
			acceleration[i] = case_def.gravity;
		density_derivative[i] = 0.0f;

		for (const auto index_pair : index_ranges) {
			for (int j = index_pair.first; j < index_pair.second; ++j) {
				const Particle& Pj = particles[j];

				const Vec3f
					r_ij = Pi.position - Pj.position,
					v_ij = Pi.velocity - Pj.velocity,
					kernel_derivative_rij = cubic_spline.gradient(r_ij);

				density_derivative[i] += case_def.particles.mass * dot_product(v_ij, kernel_derivative_rij);

				// If this is a boundary particle skip the acceleration part
				if (i >= num_of_fluid_particles)
					continue;

				const float
					&h = case_def.h,
					r2 = r_ij.length_squared(),                           // Squared distance
					Pj_pressure = beta * (std::pow(Pj.density / case_def.rhop0, gamma) - 1),
					vel_pos_dot_product = dot_product(v_ij, r_ij),
					pi_ij = [&]() {
					if (vel_pos_dot_product > 0.0f) {
						constexpr float a = 0.01f;
						const float
							rho_ij = (Pi.density + Pj.density) / 2.0f,    // Mean density
							mu = (h * vel_pos_dot_product) / (r2 + 0.01f * h * h);

						return -(a * case_def.speedsound * mu) / rho_ij;
					}
					else
						return 0.0f;
				}();

				const float
					pressure_sum = Pj_pressure + Pi_pressure,
					density_product = Pi.density * Pj.density;

				if (i < num_of_fluid_particles)
					acceleration[i] -=
					case_def.particles.mass * ((pressure_sum / density_product) + pi_ij) *
					cubic_spline.gradient(r_ij);
			}
		}
	}
}

void ParticleSystem::integrate_verlet(const float dt) {
	// How often to use the Verlet corrective step
	// TODO: consider making this variable
	constexpr int corrective_step_interval = 50;

	if (verlet_step % corrective_step_interval == 0) {
		// Corrective step
		// Fluid particles
		#pragma omp parallel for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto &Pi = particles[i];
			auto &Pi_next = next_particles[i];

			Pi_next.velocity = Pi.velocity + dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi.density + dt * density_derivative[i];
		}

		// Boundary particles
		#pragma omp parallel for
		for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i)
			next_particles[i].density = particles[i].density + dt * density_derivative[i];
	}
	else {
		// Predictor step
		#pragma omp parallel for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto
				&Pi = particles[i],
				&Pi_next = next_particles[i],
				&Pi_prev = prev_particles[i];

			Pi_next.velocity = Pi_prev.velocity + 2 * dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi_prev.density + 2 * dt * density_derivative[i];
		}
		#pragma omp parallel for
		for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i)
			next_particles[i].density = prev_particles[i].density + 2 * dt * density_derivative[i];
	}

	++verlet_step;
}

SearchGrid::cell_indices_container ParticleSystem::get_all_neighbors(const Vec3f &position) const {
	SearchGrid::cell_indices_container neighbors;
	neighbors.reserve(54);

	search_grid_fluid.get_neighbor_indices(position, neighbors);
	const auto num_of_fluid_pairs = neighbors.size();
	search_grid_boundary.get_neighbor_indices(position, neighbors);

	// Apply offset to boundary particles
	for (auto i = num_of_fluid_pairs; i < neighbors.size(); ++i) {
		neighbors[i].first += num_of_fluid_particles;
		neighbors[i].second += num_of_fluid_particles;
	}

	return neighbors;
}

void ParticleSystem::simulation_step() {
	// Sort fluid and boundary particles separately
	search_grid_fluid.sort_containers(
		std::array<SearchGrid::iter, 3>{particles.begin(), prev_particles.begin(), next_particles.begin()},
		particles.begin() + num_of_fluid_particles
	);

	search_grid_boundary.sort_containers(
		std::array<SearchGrid::iter, 3>{
			particles.begin() + num_of_fluid_particles,
			prev_particles.begin() + num_of_fluid_particles,
			next_particles.begin() + num_of_fluid_particles
		},
		particles.end()
	);

	const float time_step = calculate_time_step();
	simulation_time += time_step;
	compute_derivatives();
	integrate_verlet(time_step);

	// Move variables of next step to current and current to prev
	// Using swap instead of move assignment retains the size of next_particles
	prev_particles.swap(particles);
	particles.swap(next_particles);
}