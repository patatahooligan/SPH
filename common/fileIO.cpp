#include "stdafx.h"

#include "fileIO.h"

vtkSmartPointer<vtkPolyData> polydata_from_vtp(const std::string filename) {
	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	return reader->GetOutput();
}

vtkSmartPointer<vtkPolyData> surface_from_polydata(
	const vtkSmartPointer<vtkPolyData> polydata, const float resolution) {

	double bounds[6];
	polydata->GetBounds(bounds);

	double offset = 2 * resolution;
	bounds[0] -= offset; bounds[1] += offset;
	bounds[2] -= offset; bounds[3] += offset;
	bounds[4] -= offset; bounds[5] += offset;

	for (auto &bound : bounds)
		bound = resolution * int(bound / resolution);

	//const float resolution = 0.05f * h;
	//auto implicit_modeller = vtkSmartPointer<vtkImplicitModeller>::New();
	//implicit_modeller->SetSampleDimensions(
	//	(bounds[1] - bounds[0]) / resolution, (bounds[3] - bounds[2]) / resolution, (bounds[5] - bounds[4]) / resolution);
	//implicit_modeller->SetModelBounds(bounds);
	//implicit_modeller->SetMaximumDistance(0.02);
	//implicit_modeller->SetProcessModeToPerVoxel();
	//implicit_modeller->SetInputData(polydata);

	auto voxel_modeller = vtkSmartPointer<vtkVoxelModeller>::New();
	voxel_modeller->SetSampleDimensions(
		(bounds[1] - bounds[0]) / resolution, (bounds[3] - bounds[2]) / resolution, (bounds[5] - bounds[4]) / resolution);
	voxel_modeller->SetModelBounds(bounds);
	voxel_modeller->SetScalarTypeToFloat();
	voxel_modeller->SetMaximumDistance(0.1);
	voxel_modeller->SetInputData(polydata);

	auto flyingedges3D = vtkSmartPointer<vtkFlyingEdges3D>::New();
	flyingedges3D->SetInputConnection(voxel_modeller->GetOutputPort());
	flyingedges3D->ComputeNormalsOn();
	flyingedges3D->ComputeGradientsOn();
	flyingedges3D->SetValue(0, 0.5);

	flyingedges3D->Update();
	return flyingedges3D->GetOutput();
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