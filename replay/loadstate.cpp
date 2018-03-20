#include "stdafx.h"

#include <iomanip>

#include "loadstate.h"

void LoadState::load(ParticleContainer& target_array) {
	for (auto& particle : target_array) {
		load_particle(particle);
	}
}

void LoadState::load_particle(Particle & particle) {
	input_file
		.read(reinterpret_cast<char*>(&particle.position), sizeof(particle.position))
		.read(reinterpret_cast<char*>(&particle.velocity), sizeof(particle.velocity))
		.read(reinterpret_cast<char*>(&particle.acceleration), sizeof(particle.acceleration))
		.read(reinterpret_cast<char*>(&particle.density), sizeof(particle.density));
}