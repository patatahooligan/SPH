#include "stdafx.h"

#include "loadstate.h"

bool LoadState::load(ParticleContainer& target_array, const Mode mode) {
	target_array.resize(num_of_particles);
	for (auto& particle : target_array) {
		load_particle(particle, mode);
	}

	return !input_file.eof();
}

std::optional<ParticleContainer> LoadState::load(const Mode mode) {
	// Same functionality as the other load, to offer a choice on return method
	ParticleContainer target_array;
	load(target_array, mode);
	if (input_file.good())
		return target_array;
	else
		return std::nullopt;
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