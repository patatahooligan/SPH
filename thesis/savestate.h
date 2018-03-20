#pragma once

#include <fstream>

#include "particle.h"

class SaveState {
	public:
		// Disallow construction without a file to open, and copy-construction
		SaveState() = delete;
		SaveState(SaveState&) = delete;
		SaveState(SaveState&&) = default;

		SaveState(std::string output_filename) :
			output_file(output_filename, std::ofstream::binary | std::ofstream::trunc) {
			if (!output_file.is_open())
				throw std::runtime_error(std::string("Could not open file ") + output_filename);
		}

		void save(const ParticleContainer &data);

	private:
		std::ofstream output_file;

		void save_particle(const Particle& particle);
};