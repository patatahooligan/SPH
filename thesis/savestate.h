#pragma once

#include <fstream>

#include "particle.h"

std::ostream& operator<<(std::ostream& out, Particle& particle) {
	// TODO
	return out;
}

class SaveState {
	public:
		SaveState(std::string output_filename) :
			output_file(output_filename, std::ios_base::binary) {}

		void save(particlearray &data) {
			for (auto& particle : data)
				output_file << particle;
		}

	private:
		std::ofstream output_file;
};