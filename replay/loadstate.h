#pragma once

#include <fstream>
#include <vector>

#include "particle.h"

class LoadState {
	public:
		// Disallow construction without a file to open
		LoadState() = delete;

		LoadState(std::string input_filename);

		void load(particlearray& target_array);

		particlearray load() {
			// Same functionality as the other load, to offer a choice on return method
			particlearray target_array;
			load(target_array);
			return target_array;
		}

	private:
		std::ifstream input_file;

		void load_particle(Particle& particle);
};
