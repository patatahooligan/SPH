#pragma once

#include <stdexcept>

#include "boost\numeric\ublas\matrix.hpp"

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

	Particle(const Particle& other) {
		// There should be no code that needs this, so throw an exception to catch
		// typos that attempt to create a copy rather than use a reference of particle.
		throw std::exception("Copy constructor of particle was called");
	}
};