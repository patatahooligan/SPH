#pragma once

#include <fstream>
#include <vector>

#include "particle.h"

class LoadState {
	public:
		enum class Mode {
			Full,
			Position
		};

		// Disallow construction without a file to open, and copy-construction
		LoadState() = delete;
		LoadState(LoadState&) = delete;
		LoadState(LoadState&&) = default;

		LoadState(std::string input_filename) :
			input_file(input_filename, std::ios_base::binary)
		{
			if (!input_file.is_open())
				throw std::runtime_error(std::string("Could not open file ") + input_filename);

			input_file.read(reinterpret_cast<char*>(&num_of_fluid_particles), sizeof(num_of_fluid_particles));
			input_file.read(reinterpret_cast<char*>(&num_of_particles), sizeof(num_of_particles));
		}

		void load(ParticleContainer& target_array, const Mode mode);

		ParticleContainer load(const Mode mode) {
			// Same functionality as the other load, to offer a choice on return method
			ParticleContainer target_array;
			load(target_array, mode);
			return target_array;
		}

		auto get_num_of_fluid_particles() { return num_of_fluid_particles; }
		auto get_num_of_particles() { return num_of_particles; }

	private:
		std::ifstream input_file;
		int num_of_fluid_particles, num_of_particles;

		void load_particle(Particle& particle, const Mode mode);
};
