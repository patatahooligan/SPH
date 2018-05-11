#include "stdafx.h"

#include "savestate.h"

SaveState& SaveState::save(const ParticleContainer &data, const Mode mode) {
	for (auto& particle : data)
		save_particle(particle, mode);

	++step;
	return *this;
}

void SaveState::save_particle(const Particle& particle, const Mode mode) {
	output_file
		.write(reinterpret_cast<const char*>(&particle.position), sizeof(particle.position));

	if (mode == Mode::Full) {
		output_file
			.write(reinterpret_cast<const char*>(&particle.velocity), sizeof(particle.position))
			.write(reinterpret_cast<const char*>(&particle.density), sizeof(particle.density));
	}
}