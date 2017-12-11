#include "stdafx.h"

#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <cmath>
#include <cfloat>
#include <vector>

#include "boost\geometry\geometry.hpp"
#include "boost\numeric\ublas\matrix.hpp"

#include "physics.h"
#include "constants.h"
#include "vec3f.h"


namespace ublas = boost::numeric::ublas;

Vec3f operator*(const ublas::matrix<float> &m, const Vec3f &v) {
	// Calculate the product of a 3x3 matrix m and a 3-vector v

	if (m.size1()!=3 || m.size2() != 3) {
		throw std::invalid_argument("Matrix must be 3x3");
	}

	return Vec3f(
		m(0, 0) * v.x + m(0, 1) * v.y + m(0, 2) * v.z,
		m(1, 2) * v.x + m(1, 1) * v.y + m(1, 2) * v.z,
		m(2, 0) * v.x + m(2, 1) * v.y + m(2, 2) * v.z);
}

float trace(const ublas::matrix<float> &m) {
	const size_t size = m.size1();
	if (size != m.size2()) throw std::invalid_argument("Trace undefined for non-square matrix");
	float sum = 0;
	for (size_t i = 0; i < size; ++i) {
		sum += m(i, i);
	}
	return sum;
}

ublas::matrix<float> dyadic_product(Vec3f &v1, Vec3f &v2) {
	ublas::matrix<float> ret((size_t)3, (size_t)3);
	ret(0, 0) = v1.x * v2.x;
	ret(0, 1) = v1.x * v2.y;
	ret(0, 2) = v1.x * v2.z;
	ret(1, 0) = v1.y * v2.x;
	ret(1, 1) = v1.y * v2.y;
	ret(1, 2) = v1.y * v2.z;
	ret(2, 0) = v1.z * v2.x;
	ret(2, 1) = v1.z * v2.y;
	ret(2, 2) = v1.z * v2.z;
	return ret;
}

ublas::matrix<float> reverse(ublas::matrix<float> m) {

	const auto
		size1 = m.size1(),
		size2 = m.size2();
	ublas::matrix<float> ret(size2, size1);
	for (size_t i = 1; i < size2; ++i) {
		for (size_t j = 1; j < size1; ++j) {
			ret(i, j) = m(j, i);
		}
	}

	return ret;
}


float ParticleSystem::calculate_time_step(void) {
	float
		max_velocity_magnitude_square = 0.0f;
	for (const auto& particle : particles) {
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

	// Save the results kd-tree searches here to re-use them in the second loop
	std::vector<std::pair<size_t, float> > indices_dists[num_of_particles];

	// Calculate density
	// This has to be done for every particle before acceleration can be calculated
	#pragma omp parallel for
	for (int i = 0; i < num_of_particles; ++i) {

		// This reference is used to write cleaner equations.
		Particle& Pi = particles[i];

		// Get a pair of <index, distance> for neighbors of pi
		float position[3] = { Pi.position.x, Pi.position.y, Pi.position.z };
		kd_tree.radiusSearch(
			position, pow(smoothing_length, 2), indices_dists[i], { 32, 0.0f, false });
		assert(indices_dists[i].size() > 0);

		Pi.density = 0.0f;
		for (auto indice_dist_pair : indices_dists[i]) {
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
	for (int i = 0; i < num_of_particles; ++i) {
		// Sum of interaction forces
		Vec3f sumF{ 0.0f, 0.0f, 0.0f };
		auto h4 = pow(smoothing_length, 4);

		Particle& Pi = particles[i];
		for (auto indice_dist_pair : indices_dists[i]) {
			// If pi is pj or they are at the exact same point, do not compute interaction forces
			auto& Pj = particles[indice_dist_pair.first];
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

		Pi.acceleration = (sumF + Vec3f{ 0.0, -gravity_constant, 0.0f }) / Pi.density;
	}
}

void ParticleSystem::integrate_step() {
	// Integrate forward. Derivatives need to have already been updated through update_derivatives.

	// This function implements a leap-frog scheme where derivatives are calculated for midpoints between
	// time steps, with the exclusion of acceleration (because it is 2nd order derivative). However, we
	// use a crude approximation for velocity at integer multiples of time step because it is needed for
	// the acceleration in update_derivatives().
	const float time_step = calculate_time_step();
	simulation_time += time_step;

	Vec3f new_velocity_half[num_of_particles];

	#pragma omp parallel for
	for (int i = 0; i < num_of_particles; ++i) {
		// Integrate velocity
		new_velocity_half[i] = Vec3f(
			particles[i].velocity_half.x + particles[i].acceleration.x * time_step,
			particles[i].velocity_half.y + particles[i].acceleration.y * time_step,
			particles[i].velocity_half.z + particles[i].acceleration.z * time_step);

		assert(isfinite(new_velocity_half[i]));

		// Crude approximation for velocity at current time
		particles[i].velocity = (particles[i].velocity_half + new_velocity_half[i]) / 2;

		// Replace old velocity_half with new one
		particles[i].velocity_half = new_velocity_half[i];

		// Integrate position
		particles[i].position += particles[i].velocity_half * time_step;
		assert(isfinite(particles[i].position));
	}
}

void ParticleSystem::conflict_resolution() {
	const float damping_factor = 0.2f;

	#pragma omp parallel for
	for (int i = 0; i < num_of_particles; ++i) {
		if (particles[i].position.x > sizex) {
			particles[i].position.x = 2 * sizex - particles[i].position.x;
			if (particles[i].velocity.x > 0) {
				particles[i].velocity.x = (damping_factor-1.0f) * particles[i].velocity.x;
			}
		}
		else if (particles[i].position.x < 0) {
			particles[i].position.x = -particles[i].position.y;
			if (particles[i].velocity.x < 0) {
				particles[i].velocity.x = (damping_factor-1.0f) * particles[i].velocity.x;
			}
		}
		if (particles[i].position.y > sizey) {
			particles[i].position.y = 2 * sizey - particles[i].position.y;
			if (particles[i].velocity.y > 0) {
				particles[i].velocity.y = (damping_factor-1.0f) * particles[i].velocity.y;
			}
		}
		else if (particles[i].position.y < 0) {
			particles[i].position.y = -particles[i].position.y;
			if (particles[i].velocity.y < 0) {
				particles[i].velocity.y = (damping_factor-1.0f) * particles[i].velocity.y;
			}
		}
		if (particles[i].position.z > sizez) {
			particles[i].position.z = 2 * sizez - particles[i].position.z;
			if (particles[i].velocity.z > 0) {
				particles[i].velocity.z = (damping_factor-1.0f) * particles[i].velocity.z;
			}
		}
		else if (particles[i].position.z < 0) {
			particles[i].position.z = -particles[i].position.z;
			if (particles[i].velocity.z < 0) {
				particles[i].velocity.z = (damping_factor-1.0f) * particles[i].velocity.z;
			}
		}
	}
}


float ParticleSystem::smoothing_kernel(const Vec3f &r, const float h) {
	// Piecewise quintic smoothing kernel.
	// r is the vector argument, h is the smoothing length

	const float q2 = r.length_squared()/(h*h);

	// Figure out which part of the piecewise function we want
	if (q2 > 9) {
		// Out of range of the smoothing kernel
		return 0.0f;
	}
	else if (q2 > 4) {
		return pow(3-sqrt(q2), 5) * kernel_constant;
	}
	else if (q2 > 1) {
		const float q = sqrt(q2);
		return pow(3-q, 5) - 6 * pow(2-q, 5) * kernel_constant;
	}
	else {
		const float q = sqrt(q2);
		return pow(3-q, 5) - 6 * pow(2-q, 5) + 15 * pow(1-q, 5) * kernel_constant;
	}
}

Vec3f ParticleSystem::smoothing_kernel_derivative(const Vec3f &r, const float h) {
	const float q2 = r.length_squared()/(h*h);
	if (q2 >= 9) {
		return Vec3f(0, 0, 0);
	}
	else {
		const float
			q = sqrt(q2),
			xh = q * (h*h),
			c1 = (5.0f * pow(3-q, 4)) / xh;

		if (q >= 2){
			return Vec3f(
				c1 * r.x,
				c1 * r.y,
				c1 * r.z);
		}
		else {
			const float c2 = (-30.0f * pow(2-q, 4)) / xh;
			if (q >= 1) {
				return Vec3f(
					(c1+c2) * r.x,
					(c1+c2) * r.y,
					(c1+c2) * r.z);
			}
			else {
				const float c3 = (75.0f * pow(1-q, 4)) / xh;
				return Vec3f(
					(c1+c2+c3) * r.x,
					(c1+c2+c3) * r.y,
					(c1+c2+c3) * r.z);
			}
		}
	}
}

void ParticleSystem::randomize_particles() {
	// Coefficients to normalize rand() to [0, size].
	constexpr float
		normalizing_coefx = size / RAND_MAX,
		normalizing_coefy = size / RAND_MAX,
		normalizing_coefz = size / RAND_MAX;

	srand((unsigned int)time(NULL));
	for (auto &particle : particles) {
		particle.position.x = rand() * normalizing_coefx;
		particle.position.y = rand() * normalizing_coefy;
		particle.position.z = rand() * normalizing_coefz;
	}
}

void ParticleSystem::simulation_step() {
	kd_tree.buildIndex();

	update_derivatives();
	integrate_step();
	conflict_resolution();
}