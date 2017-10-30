#pragma once

#include <fstream>

#include "particle.h"

class SaveState {
	public:
		// Disallow construction without a file to open
		SaveState() = delete;

		SaveState(std::string output_filename) :
			output_file(output_filename, std::ios_base::binary) {
			// Write a header that indicates the format and the version number
			// This way if the format changes, backwards compatibility can be achieved.
			if (!output_file.is_open)
				throw std::runtime_error(std::string("Could not open file ") + output_filename);

			output_file << std::string("SaveState") << unsigned int{ 1 };
		}

		void save(particlearray &data);

	private:
		std::ofstream output_file;
};