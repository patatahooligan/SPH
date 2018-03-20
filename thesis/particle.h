#pragma once

#include <stdexcept>
#include <array>

#include "constants.h"
#include "vec3f.h"

class Particle {
	// A single particle of the SPH simulation
public:
	Vec3f position, velocity, velocity_half, acceleration;
	float density;

	// Would have been generated implicitly anyway, but explicit declaration is preferred
	Particle() = default;
	Particle(const Particle&) = default;
	Particle(Particle&& other) = default;
	Particle& operator=(const Particle& other) = default;
	Particle& operator=(Particle&& other) = default;
};

// This is used across translation units. Defining it here enforces consistency and facilitates changes.
using ParticleContainer = std::array<Particle, num_of_particles>;