#pragma once

#include <fstream>

#include "particle.h"
#include "loadstate.h"

class SaveState {
	public:
		enum class WriteMode {
			Append,
			Overwrite,
			DontModify
		};

		int get_step() const { return step; }

		virtual SaveState& save(ParticleConstIterator begin, ParticleConstIterator end) = 0;
		virtual ~SaveState() {}

	protected:
		int step = 0;
};


class SaveBinary : SaveState {
	public:
		using Mode = LoadState::Mode;

		SaveBinary(
			std::string output_filename, WriteMode write_mode,
			int num_of_fluid_particles, int num_of_particles);

		SaveBinary(SaveBinary&&) = default;

		SaveBinary& save(ParticleConstIterator begin, ParticleConstIterator end) override;

	private:
		std::ofstream output_file;
		Mode mode = Mode::Position;

		void save_particle(const Particle& particle);
};

class SaveVTK {
	private:
		vtkSmartPointer<vtkPoints> points;
		vtkSmartPointer<vtkPolyData> polydata;
		vtkSmartPointer<vtkCellArray> vertices;
	public:
		SaveVTK(ParticleConstIterator begin, ParticleConstIterator end);

		void save_particles(std::string output_filename);
		void save_surface(std::string output_filename);
};