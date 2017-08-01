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
	size_t size = m.size1();
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

	unsigned
		size1 = m.size1(),
		size2 = m.size2();
	ublas::matrix<float> ret(size2, size1);
	for (unsigned i = 1; i < size2; ++i) {
		for (unsigned j = 1; j < size1; ++j) {
			ret(i, j) = m(j, i);
		}
	}

	return ret;
}


float ParticleSystem::calculate_time_step(void) {
	float
		max_velocity_magnitude_square = 0.0,
		max_viscocity = 0.0f;
	for (size_t i = 0; i < num_of_particles; i++) {
		float v2 = particles[i].velocity_half.length_squared();
		if (v2 > max_velocity_magnitude_square) {
			max_velocity_magnitude_square = v2;
		}
		if (particles[i].viscocity > max_viscocity) {
			max_viscocity = particles[i].viscocity;
		}
	}
	float
		t1 = smoothing_length / (sqrt(max_velocity_magnitude_square) + speed_of_sound * speed_of_sound),
		t2 = (smoothing_length * smoothing_length) / (6 * max_viscocity);

	if (t1 < t2) 
		return t1/10.0f;
	else
		return t2/10.0f;
}

void ParticleSystem::update_derivatives() {
	// Calculate and update all derivatives of particle quantities needed to integrate
	// TODO: consider the case for zero neighbors and figure out if it needs handling

	#pragma omp parallel for
	for (int i = 0; i < num_of_particles; ++i) {
		// There are multiple quantities we need to sum for the following equations
		float
			sum_density = 0.0f,
			sum_tempererature_laplacian = 0.0f;

		Vec3f
			sum_pressure(0, 0, 0),
			sum_acceleration(0, 0, 0),
			sum_stress_derivative(0, 0, 0);

		ublas::matrix<float> velocity_tensor_derivative = ublas::zero_matrix<float>((size_t)3);

		// This reference is used to write cleaner equations.
		Particle& pi = particles[i];

		std::vector<std::pair<size_t, float>> indices_dists;
		indices_dists.reserve(30);
		float position[3] = {pi.position.x, pi.position.y, pi.position.z};
		kd_tree.radiusSearch(position, 9 * (smoothing_length*smoothing_length), indices_dists, nanoflann::SearchParams(32, 0.0f, false));

		for (auto indice_dist_pair : indices_dists) {
			const size_t j = indice_dist_pair.first;
			if (i == j) continue;		// Needs to be for every OTHER particle

			Particle& pj = particles[j];

			// Calculate their relative position and skip to next particle if it's over the smoothing kernel
			Vec3f relative_position = pi.position - pj.position;
			if (relative_position.length_squared() > smoothing_length * smoothing_length) continue;

			Vec3f
				dw = smoothing_kernel_derivative(relative_position, smoothing_length),
				relative_velocity = pi.velocity - pj.velocity;

			float vdw = relative_velocity.dot_product(dw);

			sum_density += (vdw * particle_mass) / pj.density;
			sum_pressure -=
				(particle_mass *
				(pi.pressure/pow(pi.density, 2) + pj.pressure/pow(pj.density, 2))) *
				dw;

			float vx = relative_velocity.dot_product(relative_position);
			if (vx < 0) {
				float mu = (vx*smoothing_length) / relative_position.length_squared() + 0.01f * smoothing_length * smoothing_length;
				sum_acceleration -=
					(particle_mass / pi.density) *
					(visc_b * pow(mu, 2) - visc_a * mu * speed_of_sound) /
					0.5f * (pi.density + pj.density) *
					dw;
			}

			sum_stress_derivative +=
				(particle_mass / pj.density) *
				(pi.stress_tensor + pj.stress_tensor) *
				smoothing_kernel_derivative(relative_position, smoothing_length);

			sum_tempererature_laplacian +=
				(4 * particle_mass * pi.density *
				(pi.temperature - pj.temperature) * vdw) /
				(pj.density * (pi.density + pj.density) *
				(relative_velocity.length_squared() + 0.01f * smoothing_length*smoothing_length));

			velocity_tensor_derivative -= (particle_mass * dyadic_product(relative_velocity, dw)) / pj.density;
		}
		
		pi.density_derivative = pi.density * sum_density;

		pi.acceleration =
			sum_pressure +
			sum_stress_derivative / pi.density -
			Vec3f(0.0f, gravity_constant, 0.0f) -
			sum_acceleration;

		pi.temperature_derivative = thermal_diffusion_constant * sum_tempererature_laplacian;

		ublas::matrix<float> deformation_tensor = velocity_tensor_derivative;
		deformation_tensor += reverse(velocity_tensor_derivative);
		float
			intensity_of_deformation = sqrt(pow(trace(deformation_tensor), 2)/2);

		pi.viscocity =
			(1-exp(-(jump_number + 1) * intensity_of_deformation)) *
			(1/sqrt(intensity_of_deformation) + 1/intensity_of_deformation);
	}
}

void ParticleSystem::integrate_step() {
	// Integrate forward. Derivatives need to have already been updated through update_derivatives.

	// This function implements a leap-frog scheme where derivatives are calculated for midpoints between
	// time steps, with the exclusion of acceleration (because it is 2nd order derivative). However, we
	// use a crude approximation for velocity at integer multiples of time step because it is needed for
	// the acceleration in update_derivatives().
	float time_step = calculate_time_step();
	simulation_time += time_step;

	Vec3f new_velocity_half[num_of_particles];

	#pragma omp parallel for
	for (int i = 0; i < num_of_particles; ++i) {
		// Initial approximation for new velocity_half
		new_velocity_half[i] = Vec3f(
			particles[i].velocity_half.x + particles[i].acceleration.x * time_step,
			particles[i].velocity_half.y + particles[i].acceleration.y * time_step,
			particles[i].velocity_half.z + particles[i].acceleration.z * time_step);
		
		particles[i].density += particles[i].density_derivative * time_step;

		particles[i].temperature += particles[i].temperature_derivative * time_step;

		particles[i].pressure = (speed_of_sound * speed_of_sound) * (particles[i].density - reference_density);
	}

	#pragma omp parallel for
	for (int i = 0; i < num_of_particles; ++i) {
		Vec3f velocity_correction_sum(0.0f, 0.0f, 0.0f);

		// XSPH velocity correction
		for (int j = 0; j < num_of_particles; ++j) {
			Vec3f
				relative_velocity = new_velocity_half[j] - new_velocity_half[i],
				relative_position = particles[i].position - particles[j].position;
			velocity_correction_sum += (2 * relative_velocity * smoothing_kernel(relative_position, smoothing_length)) / (particles[i].density + particles[j].density);
		}
		new_velocity_half[i] += particle_mass * xsph_coeff * velocity_correction_sum;

		// Crude approximation for velocity at current time + time step
		particles[i].velocity = (particles[i].velocity_half + new_velocity_half[i]) / 2;
		
		// Replace old velocity_half with new one
		particles[i].velocity_half = new_velocity_half[i];

		particles[i].position += particles[i].velocity_half * time_step;
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

	float q2 = r.length_squared()/(h*h);

	// Figure out which part of the piecewise function we want
	if (q2 > 9) {
		// Out of range of the smoothing kernel
		return 0.0f;
	}
	else if (q2 > 4) {
		return pow(3-sqrt(q2), 5) * kernel_constant;
	}
	else if (q2 > 1) {
		float q = sqrt(q2);
		return pow(3-q, 5) - 6 * pow(2-q, 5) * kernel_constant;
	}
	else {
		float q = sqrt(q2);
		return pow(3-q, 5) - 6 * pow(2-q, 5) + 15 * pow(1-q, 5) * kernel_constant;
	}
}

Vec3f ParticleSystem::smoothing_kernel_derivative(const Vec3f &r, const float h) {
	float q2 = r.length_squared()/(h*h);
	if (q2 >= 9) {
		return Vec3f(0, 0, 0);
	}
	else {
		float
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
			float c2 = (-30.0f * pow(2-q, 4)) / xh;
			if (q >= 1) {
				return Vec3f(
					(c1+c2) * r.x,
					(c1+c2) * r.y,
					(c1+c2) * r.z);
			}
			else {
				float c3 = (75.0f * pow(1-q, 4)) / xh;
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
	const float
		normalizing_coefx = size / RAND_MAX,
		normalizing_coefy = size / RAND_MAX,
		normalizing_coefz = size / RAND_MAX;

	srand((unsigned int)time(NULL));
	for (size_t i=0; i < num_of_particles; ++i) {
		particles[i].position.x = rand() * normalizing_coefx;
		particles[i].position.y = rand() * normalizing_coefy;
		particles[i].position.z = rand() * normalizing_coefz;
	}
}

void ParticleSystem::calculate_initial_conditions() {
	for (size_t i = 0; i < num_of_particles; ++i) {
		particles[i].density = 0.0f;
		for (size_t j = 0; j < num_of_particles; ++j) {
			particles[i].density += particle_mass * smoothing_kernel(particles[i].position - particles[j].position, smoothing_length);
		}
		particles[i].pressure = speed_of_sound * speed_of_sound * (particles[i].density - reference_density);
	}
}

void ParticleSystem::simulation_step() {
	kd_tree.buildIndex();

	update_derivatives();
	integrate_step();
	conflict_resolution();
}