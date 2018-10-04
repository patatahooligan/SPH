#include "stdafx.h"

#include "fileIO.h"

SaveVTK::SaveVTK(ParticleConstIterator begin, ParticleConstIterator end):
	points(vtkSmartPointer<vtkPoints>::New()),
	polydata(vtkSmartPointer<vtkPolyData>::New()),
	vertices(vtkSmartPointer<vtkCellArray>::New())
{
	// We pre-build the points because they are needed for either save format
	for (auto p = begin; p != end; ++p) {
		const auto id = points->InsertNextPoint(p->position.x, p->position.y, p->position.z);
		vertices->InsertNextCell(1, &id);
	}
	polydata->SetPoints(points);
	polydata->SetVerts(vertices);
}

void SaveVTK::save_particles(std::string output_filename) {
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(polydata);

	output_filename += '.';
	output_filename += writer->GetDefaultFileExtension();
	writer->SetFileName(output_filename.data());
	writer->Write();
}

void SaveVTK::save_surface(std::string output_filename, const float h) {
	double bounds[6];
	polydata->GetBounds(bounds);

	auto voxel_modeller = vtkSmartPointer<vtkVoxelModeller>::New();
	voxel_modeller->SetSampleDimensions(
		(bounds[1] - bounds[0]) / h, (bounds[3] -bounds[2]) / h, (bounds[5] - bounds[4]) / h);
	voxel_modeller->SetModelBounds(bounds);
	voxel_modeller->SetScalarTypeToFloat();
	voxel_modeller->SetMaximumDistance(0.1);
	voxel_modeller->SetInputData(polydata);

	auto marching_cubes = vtkSmartPointer<vtkMarchingCubes>::New();
	marching_cubes->SetInputConnection(voxel_modeller->GetOutputPort());
	marching_cubes->ComputeNormalsOn();
	marching_cubes->ComputeGradientsOn();
	marching_cubes->SetValue(0, 0.5);

	marching_cubes->Update();
	vtkSmartPointer<vtkPolyData> cubes_output = marching_cubes->GetOutput();

	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(cubes_output);

	output_filename += '.';
	output_filename += writer->GetDefaultFileExtension();
	writer->SetFileName(output_filename.data());
	writer->Write();
}