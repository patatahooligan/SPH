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
		a = 1.0 / (pi * pow(h, 3));

	if (q < 1.0f)
		return a * (1 - (3.0f / 2.0f) * (q*q) + (3.0f / 4.0f) * pow(q, 3));
	else if (q < 2.0f)
		return a * (pow(2 - q, 3) / 4.0f);
	else
		return 0.0f;
}


void ParticleSystem::generate_particles() {
	auto fillbox = [](const CaseDef::Box &box, ParticleContainer &particles, float density) {
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
		fillbox(box, fluid_particles, case_def.particles.density);

	for (const auto& box : case_def.bound_boxes)
		fillbox(box, boundary_particles, case_def.particles.density);
}

float ParticleSystem::calculate_time_step() {
	float
		max_velocity_magnitude_square = 0.0f;
	for (const auto& particle : fluid_particles) {
		const auto v2 = particle.velocity_half.length_squared();
		if (v2 > max_velocity_magnitude_square) {
			max_velocity_magnitude_square = v2;
		}
	}
	const float
		t1 = smoothing_length / (sqrt(max_velocity_magnitude_square) + speed_of_sound * speed_of_sound),
		t2 = (smoothing_length * smoothing_length) / (6 * viscocity);

	if (t1 < t2) 
		return t1/10.0f;
	else
		return t2/10.0f;
}

void ParticleSystem::update_derivatives() {
	// Calculate and update all derivatives of particle quantities needed to integrate
	// TODO: consider the case for zero neighbors and figure out if it needs handling

	// Calculate density
	// This has to be done for every particle before acceleration can be calculated
	#pragma omp parallel for
	for (int i = 0; i < fluid_particles.size(); ++i) {

		// This reference is used to write cleaner equations.
		Particle& Pi = fluid_particles[i];

		// Get a pair of <index, distance> for neighbors of pi
		float position[3] = { Pi.position.x, Pi.position.y, Pi.position.z };
		std::vector<std::pair<size_t, float>> indices_dists;
		kd_tree.radiusSearch(
			position, pow(smoothing_length, 2), indices_dists, { 32, 0.0f, false });
		assert(indices_dists[i].size() > 0);

		Pi.density = 0.0f;
		for (auto indice_dist_pair : indices_dists) {
			// This calculates the sum of the equation without the normalizing constant
			auto distance_squared = indice_dist_pair.second;
			Pi.density += pow((pow(smoothing_length, 2) - distance_squared), 3);
		}
		// Multiply by the normalizing constant
		assert(Pi.density > 0.0f && isfinite(Pi.density));
		Pi.density *= (4 * particle_mass) / (pi * pow(smoothing_length, 8));
	}

	// Calculate acceleration
	#pragma omp parallel for
	for (int i = 0; i < fluid_particles.size(); ++i) {

		// Get a pair of <index, distance> for neighbors of pi
		auto& Pi = fluid_particles[i];
		float position[3] = { Pi.position.x, Pi.position.y, Pi.position.z };
		std::vector<std::pair<size_t, float>> indices_dists;
		kd_tree.radiusSearch(
			position, pow(smoothing_length, 2), indices_dists, { 32, 0.0f, false });
		assert(indices_dists[i].size() > 0);

		// Sum of interaction forces
		Vec3f sumF{ 0.0f, 0.0f, 0.0f };
		auto h4 = pow(smoothing_length, 4);

		for (auto indice_dist_pair : indices_dists) {
			// If pi is pj or they are at the exact same point, do not compute interaction forces
			auto& Pj = fluid_particles[indice_dist_pair.first];
			if (Pi.position == Pj.position) continue;

			// A bunch of shorthands to make the following equations readable
			auto
				rhoi = Pi.density, rhoj = Pj.density, rho0 = reference_density,  // densities
				q_ij = sqrt(indice_dist_pair.second) / smoothing_length;         // normalized distance

			auto
				r_ij = Pi.position - Pj.position,
				v_ij = Pi.velocity - Pj.velocity;

			sumF +=
				((particle_mass * (1 - q_ij)) / (pi * h4 * rhoj)) *
				(15 * bulk_modulus * (rhoi + rhoj - 2 * rho0) * ((1 - q_ij) / q_ij) * r_ij - 40 * viscocity * v_ij);
			assert(isfinite(sumF));
		}
		// Multiply by the normalizing constant
		sumF *= particle_mass / (pi * h4);
		
		sumF += boundary_force(Pi);

		Pi.acceleration = sumF / Pi.density + Vec3f{ 0.0, -gravity_constant, 0.0f };
	}
}

void ParticleSystem::compute_forces() {
	for (int i = 0; i < fluid_particles.size(); ++i) {
		const Particle& Pi = fluid_particles[i];

		// Get neighbors of Pi
		float position[3] = { Pi.position.x, Pi.position.y, Pi.position.z };
		std::vector<std::pair<size_t, float>> indices_dists;
		kd_tree.radiusSearch(
			position, pow(2 * case_def.h, 2), indices_dists, { 32, 0.0f, false });

		constexpr int gamma = 7;
		const float
			beta = (case_def.speedsound * case_def.speedsound * case_def.rhop0) / gamma,
			Pi_pressure = beta * (std::pow(Pi.density / case_def.rhop0, gamma) - 1);

		acceleration[i] = case_def.gravity;

		for (const auto &index_distance: indices_dists) {
			const Particle& Pj = fluid_particles[index_distance.first];

			Vec3f
				r_ij = Pi.position - Pj.position,
				v_ij = Pi.velocity - Pj.velocity;

			const float
				&h = case_def.h,
				&r2 = index_distance.second,                   // Squared distance
				Pj_pressure = beta * (std::pow(Pj.density / case_def.rhop0, gamma) - 1),
				vel_pos_dot_product = dot_product(v_ij, r_ij);

			float
				pi_ij = 0.0f;                                 // Viscosity factor

			if (vel_pos_dot_product > 0.0f) {
				constexpr float a = 0.01f;
				const float
					rho_ij = (Pi.density + Pj.density) / 2.0f,    // Mean density
					mu = (h * vel_pos_dot_product) / (r2 + 0.01f * h * h);

				pi_ij = -(a * case_def.speedsound * mu) / rho_ij;
			}

			const float
				pressure_sum = Pj_pressure + Pi_pressure,
				density_product = Pi.density * Pj.density;

			acceleration[i] -=
				case_def.particles.mass * ((pressure_sum / density_product) + pi_ij) *
				smoothing_kernel_derivative(r_ij, case_def.h);
		}
	}
}

void ParticleSystem::integrate_verlet(float dt) {
	// How often to use the Verlet corrective step
	// TODO: consider making this variable
	constexpr int corrective_step_interval = 50;

	// Using swap to move the current particles to prev maintains the vector size we need
	prev_fluid_particles.swap(fluid_particles);
	fluid_particles.swap(next_fluid_particles);

	prev_boundary_particles.swap(boundary_particles);
	boundary_particles.swap(next_boundary_particles);

	if (verlet_step % corrective_step_interval) {
		// TODO!
	}
	else {
		for (int i = 0; i < fluid_particles.size(); ++i) {
			auto
				&pi = fluid_particles[i],
				&pi_next = next_fluid_particles[i],
				&pi_prev = prev_fluid_particles[i];

			pi_next.velocity = pi_prev.velocity + 2 * dt * acceleration[i];
			pi_next.position = pi.position + dt * pi.velocity + 0.5f * dt * dt * acceleration[i];
			pi_next.density = pi_prev.density + 2 * dt * density_derivative[i];
		}
	}

	++verlet_step;
}

void ParticleSystem::conflict_resolution() {
	// If particles have intersected the wall, bring them back out to the surface.
	// Note that this does not implement a boundary condition because it is intented
	// to be called when integration position, whereas boundary conditions should be
	// calculated when calculating forces & acceleration.

	#pragma omp parallel for
	for (int i = 0; i < fluid_particles.size(); ++i) {
		fluid_particles[i].position.x = std::max(fluid_particles[i].position.x, 0.0f);
		fluid_particles[i].position.x = std::min(fluid_particles[i].position.x, size);

		fluid_particles[i].position.y = std::max(fluid_particles[i].position.y, 0.0f);
		fluid_particles[i].position.y = std::min(fluid_particles[i].position.y, size);

		fluid_particles[i].position.z = std::max(fluid_particles[i].position.z, 0.0f);
		fluid_particles[i].position.z = std::min(fluid_particles[i].position.z, size);
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

	update_derivatives();
	integrate_step();
}