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
	auto fillbox = [this](const CaseDef::Box &box, ParticleContainer& target_container) {
		// Because we might need to fill any combination of faces of a box, break it down
		struct TargetBox { Vec3f origin; int x_increments, y_increments, z_increments; };

		const float& density = case_def.particles.density;

		const std::vector<TargetBox> targets = [&]() {
			const int
				x_increments = int(box.size.x / density) + 1,
				y_increments = int(box.size.y / density) + 1,
				z_increments = int(box.size.z / density) + 1;
			constexpr int thickness = 1;

			std::vector<TargetBox> targets;

			if (box.fillmode.solid) {
				// If box is solid, skip all other checks and return just this
				targets.push_back({ box.origin, x_increments, y_increments, z_increments });
				return targets;
			}

			// Otherwise, define a box for every face needed
			const Vec3f
				right_origin = box.origin + Vec3f{ box.size.x, 0.0f, 0.0f },
				back_origin = box.origin + Vec3f{ 0.0f, box.size.y, 0.0f },
				top_origin = box.origin + Vec3f{ 0.0f, 0.0f, box.size.z };
			if (box.fillmode.left) targets.push_back({ box.origin, thickness, y_increments, z_increments });
			if (box.fillmode.right) targets.push_back({ right_origin, thickness, y_increments, z_increments });
			if (box.fillmode.front) targets.push_back({ box.origin, x_increments, thickness, z_increments });
			if (box.fillmode.back) targets.push_back({ back_origin, x_increments, thickness, z_increments });
			if (box.fillmode.bottom) targets.push_back({ box.origin, x_increments, y_increments, thickness });
			if (box.fillmode.top) targets.push_back({ top_origin, x_increments, y_increments, thickness });
			return targets;
		} ();

		for (const auto& target : targets) {
			for (int x = 0; x < target.x_increments; ++x) {
				for (int y = 0; y < target.y_increments; ++y) {
					for (int z = 0; z < target.z_increments; ++z) {
						Particle p;
						p.position = box.origin + Vec3f{ x * density, y * density, z * density };
						p.density = case_def.rhop0;
						target_container.emplace_back(p);
					}
				}
			}
		}
	};

	auto remove_particles = [this](const CaseDef::Box &box, ParticleContainer& target_container) {
		target_container.erase(std::remove_if(target_container.begin(), target_container.end(),
			[&box](const Particle& p) {
				return box.contains(p.position);
			}), target_container.end());
	};

	ParticleContainer fluid_particles, boundary_particles;

	for (const auto& box : case_def.particle_boxes) {
		using Type = CaseDef::Box::Type;
		switch (box.type) {
		case Type::Fluid:
			fillbox(box, fluid_particles);
			break;
		case Type::Boundary:
			fillbox(box, boundary_particles);
			break;
		case Type::Void:
			remove_particles(box, fluid_particles);
			remove_particles(box, boundary_particles);
			break;
		}
	}

	// Make note of the number of fluid particles, then concatenate the two vectors
	// [0, num_of_particles)   -> fluid
	// [num_of_particles, end) -> boundary
	num_of_fluid_particles = fluid_particles.size();
	particles = std::move(fluid_particles);
	particles.insert(particles.end(), boundary_particles.begin(), boundary_particles.end());
	particles.shrink_to_fit();
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
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; size_t(i) < particles.size(); ++i) {
			constexpr int gamma = 7;
			const float
				beta = (case_def.speedsound * case_def.speedsound * case_def.rhop0) / gamma;

			pressure[i] = beta * (std::pow(particles[i].density / case_def.rhop0, gamma) - 1);

		}

		#pragma omp for
		for (int i = 0; size_t(i) < particles.size(); ++i) {
			const Particle& Pi = particles[i];

			// Get neighbors of Pi
			const auto index_ranges = get_all_neighbors(Pi.position);

			if (i < num_of_fluid_particles)
				acceleration[i] = case_def.gravity;
			density_derivative[i] = 0.0f;

			for (const auto index_pair : index_ranges) {
				for (int j = index_pair.first; j < index_pair.second; ++j) {
					const Particle& Pj = particles[j];

					const Vec3f
						r_ij = Pi.position - Pj.position,
						v_ij = Pi.velocity - Pj.velocity,
						kernel_gradient_rij = cubic_spline.gradient(r_ij);

					density_derivative[i] += case_def.particles.mass * dot_product(v_ij, kernel_gradient_rij);

					// If this is a boundary particle skip the acceleration part
					if (i >= num_of_fluid_particles)
						continue;

					const float
						&h = case_def.h,
						r2 = r_ij.length_squared(),                           // Squared distance
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
						&Pi_pressure = pressure[i],
						&Pj_pressure = pressure[j],
						pressure_sum = Pj_pressure + Pi_pressure,
						density_product = Pi.density * Pj.density;

					if (i < num_of_fluid_particles)
						acceleration[i] -=
						case_def.particles.mass * ((pressure_sum / density_product) + pi_ij) *
						kernel_gradient_rij;
				}
			}
		}
	}
}

void ParticleSystem::integrate_verlet(const float dt) {
	// How often to use the Verlet corrective step
	// TODO: consider making this variable
	constexpr int corrective_step_interval = 50;

	#pragma omp parallel
	if (verlet_step % corrective_step_interval == 0) {
		// Corrective step
		// Fluid particles
		#pragma omp for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto &Pi = particles[i];
			auto &Pi_next = next_particles[i];

			Pi_next.velocity = Pi.velocity + dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi.density + dt * density_derivative[i];
		}

		// Boundary particles
		#pragma omp for
		for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i)
			next_particles[i].density = particles[i].density + dt * density_derivative[i];
	}
	else {
		// Predictor step
		#pragma omp for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto
				&Pi = particles[i],
				&Pi_next = next_particles[i],
				&Pi_prev = prev_particles[i];

			Pi_next.velocity = Pi_prev.velocity + 2 * dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi_prev.density + 2 * dt * density_derivative[i];
		}
		#pragma omp for
		for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i)
			next_particles[i].density = prev_particles[i].density + 2 * dt * density_derivative[i];
	}

	++verlet_step;
}

SearchGrid::static_cell_indices_container ParticleSystem::get_all_neighbors(const Vec3f &position) const {
	auto neighbors = search_grid_fluid.get_neighbor_indices(position);
	const auto boundary_start = neighbors.size();
	search_grid_boundary.get_neighbor_indices(position, neighbors);

	// Apply offset to boundary particles
	for (auto i = boundary_start; i < neighbors.size(); ++i) {
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