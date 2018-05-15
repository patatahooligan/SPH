#pragma once

#include <fstream>

#include "particle.h"

class SaveState {
	public:
		enum class Mode {
			Full,
			Position
		};

		// Disallow construction without a file to open, and copy-construction
		SaveState() = delete;
		SaveState(SaveState&) = delete;
		SaveState(SaveState&&) = default;

		SaveState(std::string_view output_filename, int num_of_fluid_particles, int num_of_particles) :
			output_file(output_filename.data(), std::ofstream::binary | std::ofstream::trunc) {
			if (!output_file.is_open())
				throw std::runtime_error(std::string("Could not open file ") + output_filename.data());

			output_file.write(reinterpret_cast<const char*>(&num_of_fluid_particles), sizeof(num_of_fluid_particles));
			output_file.write(reinterpret_cast<const char*>(&num_of_particles), sizeof(num_of_particles));
		}

		SaveState& save(const ParticleContainer &data, const Mode mode);

		int get_step() const { return step; }

	private:
		std::ofstream output_file;
		int step = 0;

		void save_particle(const Particle& particle, const Mode mode);
};