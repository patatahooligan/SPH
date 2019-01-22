#pragma once

#include "particle.h"

//vtkSmartPointer<vtkPolyData> polydata_from_vtp(std::string filename);
//
//vtkSmartPointer<vtkPolyData> surface_from_polydata(vtkSmartPointer<vtkPolyData> polydata, float h);
//
//void save_surface_to_vtp(vtkSmartPointer<vtkPolyData> surface, std::string output_filename);

class SaveVTK {
	private:
		vtkSmartPointer<vtkPoints> points;
		vtkSmartPointer<vtkPolyData> polydata;
		vtkSmartPointer<vtkCellArray> vertices;
	public:
		SaveVTK(ParticleConstIterator begin, ParticleConstIterator end);

		void save_particles(std::string output_filename);/*
		void save_surface(std::string output_filename, float h);*/
};