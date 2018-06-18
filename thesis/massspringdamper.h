#pragma once
#include "stdafx.h"

#include "vec3f.h"
#include "particle.h"

struct MassSpringDamper {
	inline static float resting_length, k, damping_coef;
	std::pair<int, int> particle_indices;

	static Vec3f compute_force(const Particle &target, const Particle &other) {
		const Vec3f
			relative_position = target.position - other.position,
			relative_velocity = target.velocity - other.velocity;
		const float length = relative_position.length();
		const Vec3f displacement = (length - resting_length) * (relative_position / length);
		return -k * displacement; // -damping_coef * relative_velocity;
	}
};
