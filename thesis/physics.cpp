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

auto radius_search(const ParticleKDTree &kd_tree, const Vec3f &center, const float radius) {
	float position[3] = { center.x, center.y, center.z };
	std::vector<std::pair<size_t, float>> indices_dists;
	kd_tree.radiusSearch(
		position, pow(2 * radius, 2), indices_dists, { 32, 0.0f, false });
	return indices_dists;
}

float piecewise_smoothing_kernel(const Vec3f &r, const float h) {
	// Piecewise quintic smoothing kernel.
	// r is the vector argument, h is the smoothing length

	const float q2 = r.length_squared() / (h*h);

	// Figure out which part of the piecewise function we want
	if (q2 > 9) {
		// Out of range of the smoothing kernel
		return 0.0f;
	}
	else if (q2 > 4) {
		return pow(3 - sqrt(q2), 5) * kernel_constant;
	}
	else if (q2 > 1) {
		const float q = sqrt(q2);
		return pow(3 - q, 5) - 6 * pow(2 - q, 5) * kernel_constant;
	}
	else {
		const float q = sqrt(q2);
		return pow(3 - q, 5) - 6 * pow(2 - q, 5) + 15 * pow(1 - q, 5) * kernel_constant;
	}
}

Vec3f piecewise_smoothing_kernel_derivative(const Vec3f &r, const float h) {
	const float q2 = r.length_squared() / (h*h);
	if (q2 >= 9) {
		return Vec3f(0, 0, 0);
	}
	else {
		const float
			q = sqrt(q2),
			xh = q * (h*h),
			c1 = (5.0f * pow(3 - q, 4)) / xh;

		if (q >= 2) {
			return Vec3f(
				c1 * r.x,
				c1 * r.y,
				c1 * r.z);
		}
		else {
			const float c2 = (-30.0f * pow(2 - q, 4)) / xh;
			if (q >= 1) {
				return Vec3f(
					(c1 + c2) * r.x,
					(c1 + c2) * r.y,
					(c1 + c2) * r.z);
			}
			else {
				const float c3 = (75.0f * pow(1 - q, 4)) / xh;
				return Vec3f(
					(c1 + c2 + c3) * r.x,
					(c1 + c2 + c3) * r.y,
					(c1 + c2 + c3) * r.z);
			}
		}
	}
}

float cubic_spline(const Vec3f &r, const float h) {
	assert(h >= 0.0f);
	const float
		q = std::sqrt(r.length_squared()) / h,
		a = 1.0f / (pi * pow(h, 3));

	if (q < 1.0f)
		return a * (1 - (3.0f / 2.0f) * (q*q) + (3.0f / 4.0f) * pow(q, 3));
	else if (q < 2.0f)
		return a * (pow(2 - q, 3) / 4.0f);
	else
		return 0.0f;
}


void ParticleSystem::generate_particles() {
	auto fillbox = [](const CaseDef::Box &box, ParticleContainer &particles, const float density) {
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
					particles.emplace_back(p);
				}
			}
		}
	};

	for (const auto& box : case_def.fluid_boxes)
		fillbox(box, particles, case_def.particles.density);

	// Storing fluid and bound particles in the same container greatly simplifies the iterative
	// algorithms. For the parts that apply specifically to one type of particle, we need to keep
	// track of the number of fluid particles and iterate accordingly.
	num_of_fluid_particles = particles.size();

	for (const auto& box : case_def.bound_boxes)
		fillbox(box, particles, case_def.particles.density);
}

float ParticleSystem::calculate_time_step() const {
	const float
		&h = case_def.h,
		&c = case_def.speedsound;

	float
		acceleration2_max = 0.0f,
		dt_cv = std::numeric_limits<float>::max();

	for (int i = 0; i < num_of_fluid_particles; ++i) {
		const Particle& Pi = particles[i];
		const auto indices_dists = radius_search(kd_tree, Pi.position, case_def.h);

		float sigma = 0.0f;
		for (const auto &index_distance : indices_dists) {
			const Particle& Pj = particles[index_distance.first];
			const Vec3f
				v_ij = Pj.velocity - Pi.velocity,
				r_ij = Pj.position - Pi.position;

			const float
				r2 = r_ij.length_squared();

			sigma = std::max(
				sigma,
				std::abs((case_def.h * dot_product(v_ij, r_ij)) / (r2 + 0.01f * h * h)));
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
	for (int i = 0, size = particles.size(); i < size; ++i) {
		const Particle& Pi = particles[i];

		// Get neighbors of Pi
		float position[3] = { Pi.position.x, Pi.position.y, Pi.position.z };
		const auto indices_dists = radius_search(kd_tree, Pi.position, case_def.h);

		constexpr int gamma = 7;
		const float
			beta = (case_def.speedsound * case_def.speedsound * case_def.rhop0) / gamma,
			Pi_pressure = beta * (std::pow(Pi.density / case_def.rhop0, gamma) - 1);

		acceleration[i] = case_def.gravity;
		density_derivative[i] = 0.0f;

		for (const auto &index_distance: indices_dists) {
			const Particle& Pj = particles[index_distance.first];

			const Vec3f
				r_ij = Pi.position - Pj.position,
				v_ij = Pi.velocity - Pj.velocity,
				kernel_derivative_rij = smoothing_kernel_derivative(r_ij, case_def.h);

			density_derivative[i] += case_def.particles.mass * dot_product(v_ij, kernel_derivative_rij);

			// If this is a boundary particle skip the acceleration part
			if (index_distance.first >= size_t(num_of_fluid_particles))
				continue;

			const float
				&h = case_def.h,
				&r2 = index_distance.second,                   // Squared distance
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

			acceleration[i] -=
				case_def.particles.mass * ((pressure_sum / density_product) + pi_ij) *
				smoothing_kernel_derivative(r_ij, case_def.h);
		}
	}
}

void ParticleSystem::integrate_verlet(const float dt) {
	// How often to use the Verlet corrective step
	// TODO: consider making this variable
	constexpr int corrective_step_interval = 50;

	if (verlet_step % corrective_step_interval) {
		// Fluid particles
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto &Pi = particles[i];
			auto &Pi_next = next_particles[i];

			Pi_next.velocity = Pi.velocity + dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi.density + dt * density_derivative[i];
		}

		// Boundary particles
		for (int i = num_of_fluid_particles, size = particles.size(); i < size; ++i)
			next_particles[i].density = particles[i].density + dt * density_derivative[i];
	}
	else {
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto
				&Pi = particles[i],
				&Pi_next = next_particles[i],
				&Pi_prev = prev_particles[i];

			Pi_next.velocity = Pi_prev.velocity + 2 * dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi_prev.density + 2 * dt * density_derivative[i];
		}
		for (int i = num_of_fluid_particles, size = particles.size(); i < size; ++i)
			next_particles[i].density = prev_particles[i].density + 2 * dt * density_derivative[i];
	}

	++verlet_step;
}

void ParticleSystem::conflict_resolution() {
	// If particles have intersected the wall, bring them back out to the surface.
	// Note that this does not implement a boundary condition because it is intented
	// to be called when integration position, whereas boundary conditions should be
	// calculated when calculating forces & acceleration.

	#pragma omp parallel for
	for (int i = 0; i < num_of_fluid_particles; ++i) {
		particles[i].position.x = std::max(particles[i].position.x, 0.0f);
		particles[i].position.x = std::min(particles[i].position.x, size);

		particles[i].position.y = std::max(particles[i].position.y, 0.0f);
		particles[i].position.y = std::min(particles[i].position.y, size);

		particles[i].position.z = std::max(particles[i].position.z, 0.0f);
		particles[i].position.z = std::min(particles[i].position.z, size);
	}
}

Vec3f ParticleSystem::boundary_force(const Particle& p) {

	// The fictitious particle that enforces the boundary condition should be slightly offset
	// so as to not intersect with the particles.
	constexpr float wall_offset = smoothing_length / 100.0f;

	auto lennard_jones_force = [](float r) {
		constexpr float
			delta = 1.0f,
			r_0 = smoothing_length;

		return (delta * (pow(r_0 / r, 12) - pow(r_0 / r, 6))) / r;
	};

	Vec3f F{ 0.0f, 0.0f, 0.0f };

	if (p.position.x > smoothing_length)
		F.x += lennard_jones_force(p.position.x);
	else if (p.position.x > size - smoothing_length)
		F.x -= lennard_jones_force(size - p.position.x);

	if (p.position.y > smoothing_length)
		F.y += lennard_jones_force(p.position.y);
	else if (p.position.y > size - smoothing_length)
		F.y -= lennard_jones_force(size - p.position.y);

	if (p.position.z > smoothing_length)
		F.z += lennard_jones_force(p.position.z);
	else if (p.position.z > size - smoothing_length)
		F.z -= lennard_jones_force(size - p.position.z);

	return F;
}

void ParticleSystem::simulation_step() {
	kd_tree.buildIndex();

	const float time_step = calculate_time_step();
	simulation_time += time_step;
	compute_derivatives();
	integrate_verlet(time_step);

	// Move variables of next step to current and current to prev
	// Using swap instead of move assignment retains the size of next_particles
	prev_particles.swap(particles);
	particles.swap(next_particles);
}