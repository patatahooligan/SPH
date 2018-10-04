#pragma once

#include "particle.h"

class SaveVTK {
	private:
		vtkSmartPointer<vtkPoints> points;
		vtkSmartPointer<vtkPolyData> polydata;
		vtkSmartPointer<vtkCellArray> vertices;
	public:
		SaveVTK(ParticleConstIterator begin, ParticleConstIterator end);

		void save_particles(std::string output_filename);
		void save_surface(std::string output_filename, float h);
};