#pragma once

#include <stdexcept>

#include "boost\numeric\ublas\matrix.hpp"

#include "constants.h"
#include "vec3f.h"

class Particle {
	// A single particle of the SPH simulation
public:
	Vec3f position, velocity, velocity_half, acceleration;
	float
		density, density_derivative,
		pressure,
		temperature, temperature_derivative,
		viscocity;
	boost::numeric::ublas::matrix<float> stress_tensor;

	Particle() :
	stress_tensor (boost::numeric::ublas::zero_matrix<float>((size_t)3)),
	velocity(0.0f, 0.0f, 0.0f),
	velocity_half(0.0f, 0.0f, 0.0f),
	acceleration(0.0f, 0.0f, 0.0f) {
	}

	// Do no permit copying of Particle as it does not currently serve any purpose so using
	// it might be accidental and lead to bugs.
	Particle(const Particle& other) = delete;
	Particle& operator=(const Particle& other) = delete;
};

// This is used across translation units. Defining it here enforces consistency and facilitates changes.
using particlearray = std::array<Particle, num_of_particles>;