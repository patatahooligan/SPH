#include "stdafx.h"

#include <iomanip>

#include "loadstate.h"

void LoadState::load(ParticleContainer& target_array, const Mode mode) {
	target_array.resize(num_of_particles);
	for (auto& particle : target_array) {
		load_particle(particle, mode);
	}
}

void LoadState::load_particle(Particle & particle, const Mode mode) {
	input_file
		.read(reinterpret_cast<char*>(&particle.position), sizeof(particle.position));

	if (mode == Mode::Full) {
		input_file
			.read(reinterpret_cast<char*>(&particle.velocity), sizeof(particle.position))
			.read(reinterpret_cast<char*>(&particle.density), sizeof(particle.density));
	}
}