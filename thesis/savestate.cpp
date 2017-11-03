#include "stdafx.h"

#include "savestate.h"

std::ostream& operator<<(std::ostream& out, const Vec3f &vec) {
	out << vec.x << vec.y << vec.z;
	return out;
}

std::ostream& operator<<(std::ostream& out, const Particle& particle) {
	out <<
		particle.position <<
		particle.velocity <<
		particle.acceleration <<
		particle.temperature <<
		particle.temperature_derivative <<
		particle.density <<
		particle.density_derivative <<
		particle.temperature <<
		particle.temperature_derivative <<
		particle.viscocity;
	return out;
}

void SaveState::save(const particlearray &data) {
	for (auto& particle : data)
		output_file << particle;
}