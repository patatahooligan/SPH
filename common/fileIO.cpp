#include "stdafx.h"

#include "fileIO.h"

vtkSmartPointer<vtkPolyData> polydata_from_vtp(const std::string filename) {
	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	return reader->GetOutput();
}

vtkSmartPointer<vtkPolyData> surface_from_polydata(
	const vtkSmartPointer<vtkPolyData> polydata, const float h) {

	double bounds[6];
	polydata->GetBounds(bounds);
	bounds[0] -= h; bounds[1] += h;
	bounds[2] -= h; bounds[3] += h;
	bounds[4] -= h; bounds[5] += h;

	const float resolution = 0.75f * h;

	auto voxel_modeller = vtkSmartPointer<vtkVoxelModeller>::New();
	voxel_modeller->SetSampleDimensions(
		(bounds[1] - bounds[0]) / resolution, (bounds[3] - bounds[2]) / resolution, (bounds[5] - bounds[4]) / resolution);
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
	return marching_cubes->GetOutput();
}

void save_surface_to_vtp(const vtkSmartPointer<vtkPolyData> surface, std::string output_filename) {
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(surface);

	output_filename += '.';
	output_filename += writer->GetDefaultFileExtension();
	writer->SetFileName(output_filename.data());
	writer->Write();
}

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
	const auto surface = surface_from_polydata(polydata, h);
	save_surface_to_vtp(surface, output_filename);
}