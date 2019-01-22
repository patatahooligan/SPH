#include "stdafx.h"

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

	auto implicit_modeller = vtkSmartPointer<vtkImplicitModeller>::New();
	implicit_modeller->SetSampleDimensions(
		(bounds[1] - bounds[0]) / resolution, (bounds[3] - bounds[2]) / resolution, (bounds[5] - bounds[4]) / resolution);
	implicit_modeller->SetModelBounds(bounds);
	implicit_modeller->SetMaximumDistance(0.014);
	implicit_modeller->SetProcessModeToPerVoxel();
	implicit_modeller->SetInputData(polydata);
	implicit_modeller->SetOutputScalarTypeToFloat();

	auto flyingedges3D = vtkSmartPointer<vtkFlyingEdges3D>::New();
	flyingedges3D->SetInputConnection(implicit_modeller->GetOutputPort());
	flyingedges3D->SetValue(0, 0.5f);

	auto smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
	smoother->SetInputConnection(flyingedges3D->GetOutputPort());
	smoother->SetNumberOfIterations(15);
	smoother->BoundarySmoothingOff();
	smoother->FeatureEdgeSmoothingOff();
	smoother->SetPassBand(.01);
	smoother->NonManifoldSmoothingOn();
	smoother->NormalizeCoordinatesOn();
	smoother->Update();

	return smoother->GetOutput();
}

void save_surface_to_vtp(const vtkSmartPointer<vtkPolyData> surface, std::string output_filename) {
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(surface);

	output_filename += '.';
	output_filename += writer->GetDefaultFileExtension();
	writer->SetFileName(output_filename.data());
	writer->Write();
}

int main(int argc, char* argv[]) {
	if (argc < 3)
		return EXIT_FAILURE;

	constexpr float resolution = 0.006;
	const std::string
		input_prefix = argv[1],
		output_prefix = argv[2];

	const int	begin_index = argc >= 4 ? std::stoi(argv[3]) : 0;
	int end_index = begin_index;
	while (std::filesystem::exists(input_prefix + std::to_string(end_index) + ".vtp"))
		++end_index;

	for (int i = begin_index; i < end_index; ++i) {
		const auto polydata = polydata_from_vtp(input_prefix + std::to_string(i) + ".vtp");
		const auto surface = surface_from_polydata(polydata, resolution);
		const auto output_filename = output_prefix + std::to_string(i);
		save_surface_to_vtp(surface, output_filename);

		std::cout << "Saved output to " << output_filename << '\n';
	}
}