#pragma once
#include "stdafx.h"

#include "vec3f.h"
#include "particle.h"

struct MassSpringDamper {
	inline static float damping_coef;
	std::pair<int, int> particle_indices;
	float resting_length;

	Vec3f compute_force(const Particle &target, const Particle &other, const float k) const {
		const Vec3f
			relative_position = target.position - other.position,
			relative_velocity = target.velocity - other.velocity;
		const float length = relative_position.length();
		const Vec3f
			spring_unit_vector = relative_position / length,
			displacement = (length - resting_length) * spring_unit_vector;
		return -k * displacement - damping_coef * relative_velocity;
	}
};

using MassSpringContainer = std::vector<MassSpringDamper>;
using MassSpringIterator = MassSpringContainer::iterator;
using MassSpringConstIterator = MassSpringContainer::const_iterator;

struct MassSpringSystem {
	MassSpringContainer springs;
	float initial_k, start_of_melting, duration_of_melting;
};
