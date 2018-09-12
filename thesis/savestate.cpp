#include "stdafx.h"

#include "savestate.h"

SaveBinary::SaveBinary(
		const std::string output_filename, const WriteMode write_mode,
		const int num_of_fluid_particles, const int num_of_particles) {
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

SaveBinary& SaveBinary::save(const ParticleConstIterator begin, const ParticleConstIterator end) {
	for (auto particle_it = begin; particle_it != end; ++particle_it)
		save_particle(*particle_it);

	++step;
	return *this;
}

void SaveBinary::save_particle(const Particle& particle) {
	output_file
		.write(reinterpret_cast<const char*>(&particle.position), sizeof(particle.position));

	if (mode == Mode::Full) {
		output_file
			.write(reinterpret_cast<const char*>(&particle.velocity), sizeof(particle.position))
			.write(reinterpret_cast<const char*>(&particle.density), sizeof(particle.density));
	}
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