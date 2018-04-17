#include "stdafx.h"

#include "savestate.h"

void SaveState::save(const ParticleContainer &data) {
	for (auto& particle : data)
		save_particle(particle);
}

void SaveState::save_particle(const Particle& particle) {
	output_file
		.write(reinterpret_cast<const char*>(&particle.position), sizeof(particle.position));
}