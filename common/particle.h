#pragma once

#include "../common/constants.h"
#include "vec3f.h"

struct Particle {
	// A single particle of the SPH simulation
	Vec3f position, velocity;
	float density;
};

// This is used across translation units. Defining it here enforces consistency and facilitates changes.
using ParticleContainer = std::vector<Particle>;
using ParticleIterator = ParticleContainer::iterator;
using ParticleConstIterator = ParticleContainer::const_iterator;