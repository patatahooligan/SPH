#pragma once

#include <fstream>

#include "particle.h"
#include "loadstate.h"

class SaveState {
	public:
		using Mode = LoadState::Mode;
		enum class WriteMode {
			Append,
			Overwrite,
			DontModify
		};

		// Disallow construction without a file to open, and copy-construction
		SaveState() = delete;
		SaveState(SaveState&) = delete;
		SaveState(SaveState&&) = default;

		SaveState(
			const std::string output_filename, const WriteMode write_mode,
			const int num_of_fluid_particles, const int num_of_particles)
		{
			bool append_to_file = false;
			int old_num_of_fluid_particles, old_num_of_particles;
			if (write_mode != WriteMode::Overwrite) {
				std::ifstream old_file(output_filename, std::ios::binary);
				if (old_file.good()) {
					if (write_mode == WriteMode::DontModify)
						throw std::runtime_error("File exists and no write mode specified");
					append_to_file = true;
				}

				if (append_to_file) {
					old_file.read(reinterpret_cast<char*>(&old_num_of_fluid_particles),
						sizeof(old_num_of_fluid_particles));
					old_file.read(reinterpret_cast<char*>(&old_num_of_particles), sizeof(old_num_of_particles));
				}
			}

			const auto open_mode = [&]() {
				if (write_mode == WriteMode::Overwrite)
					return std::ofstream::trunc;
				if (write_mode == WriteMode::Append)
					return std::ofstream::app;
			}();

			output_file.open(output_filename, std::ofstream::binary | open_mode);

			if (!output_file.is_open())
				throw std::runtime_error(std::string("Could not open file ") + output_filename.data());

			if (append_to_file) {
				if (old_num_of_fluid_particles != num_of_fluid_particles ||
					old_num_of_particles != num_of_particles)
					throw std::runtime_error(
						"Previous binary file has different number of particles than current simulation");
			}
			else {
				output_file.write(reinterpret_cast<const char*>(&num_of_fluid_particles), sizeof(num_of_fluid_particles));
				output_file.write(reinterpret_cast<const char*>(&num_of_particles), sizeof(num_of_particles));
			}
		}

		SaveState& save(const ParticleContainer &data, const Mode mode);

		int get_step() const { return step; }

	private:
		std::ofstream output_file;
		int step = 0;

		void save_particle(const Particle& particle, const Mode mode);
};