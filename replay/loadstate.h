#pragma once

#include <fstream>
#include <vector>

#include "particle.h"

class LoadState {
	public:
		// Disallow construction without a file to open, and copy-construction
		LoadState() = delete;
		LoadState(LoadState&) = delete;
		LoadState(LoadState&&) = default;

		LoadState(std::string input_filename) :
			input_file(input_filename, std::ios_base::binary)
		{
			if (!input_file.is_open())
				throw std::runtime_error(std::string("Could not open file ") + input_filename);

			input_file.read(reinterpret_cast<char*>(num_of_particles_m), sizeof(num_of_particles_m));
		}

		void load(ParticleContainer& target_array);

		ParticleContainer load() {
			// Same functionality as the other load, to offer a choice on return method
			ParticleContainer target_array;
			load(target_array);
			return target_array;
		}

		auto num_of_particles() { return num_of_particles_m; }

	private:
		std::ifstream input_file;
		int num_of_particles_m;

		void load_particle(Particle& particle);
};
