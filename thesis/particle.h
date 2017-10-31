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
	stress_tensor (boost::numeric::ublas::zero_matrix<float>((size_t)3)) {}

	// Would have been generated implicitly anyway, but explicit declaration is preferred
	Particle(const Particle& other) = default;
	Particle(Particle&& other) = default;
	Particle& operator=(const Particle& other) = default;
	Particle& operator=(Particle&& other) = default;
};

// This is used across translation units. Defining it here enforces consistency and facilitates changes.
using particlearray = std::array<Particle, num_of_particles>;