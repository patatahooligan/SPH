#include "stdafx.h"

#include "savestate.h"

void SaveState::save(const particlearray &data) {
	for (auto& particle : data)
		save_particle(particle);
}

void SaveState::save_particle(const Particle& particle) {
	output_file
		.write(reinterpret_cast<const char*>(&particle.position), sizeof(particle.position))
		.write(reinterpret_cast<const char*>(&particle.velocity), sizeof(particle.velocity))
		.write(reinterpret_cast<const char*>(&particle.acceleration), sizeof(particle.acceleration))
		.write(reinterpret_cast<const char*>(&particle.density), sizeof(particle.density));
}