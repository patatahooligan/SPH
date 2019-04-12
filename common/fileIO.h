#pragma once

#include <string>

#include "particle.h"

inline void save_particles(ParticleConstIterator begin, ParticleConstIterator end, std::string output_filename) {
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto polydata = vtkSmartPointer<vtkPolyData>::New();
	auto vertices = vtkSmartPointer<vtkCellArray>::New();

	for (auto p = begin; p != end; ++p) {
		const auto id = points->InsertNextPoint(p->position.x, p->position.y, p->position.z);
		vertices->InsertNextCell(1, &id);
	}

	polydata->SetPoints(points);
	polydata->SetVerts(vertices);

	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(polydata);

	output_filename += '.';
	output_filename += writer->GetDefaultFileExtension();
	writer->SetFileName(output_filename.data());
	writer->Write();
}