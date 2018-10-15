// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "physics.h"
#include "fileIO.h"
#include "XMLReader.h"

struct RunOptions {
	std::string case_filename;
	float time, output_period;

	std::optional<int> max_run_time;
	std::optional<std::string> input;
	std::optional<std::string> particles_output_filename, surface_output_filename;
};

auto get_options(int argc, char **argv) {
	cxxopts::Options cxx_options("SPH", "Smoothed Particle Hydrodynamics");
	RunOptions run_options;

	cxx_options.add_options()
		("particles-output", "Prefix for group of vtk files to hold the result of the simulation",
			cxxopts::value<std::string>())
		("surface-output", "Prefix for group of vtk files to hold the reconstructed surface",
			cxxopts::value<std::string>())
		("t, time", "Time of simulation run in seconds",
			cxxopts::value<float>())
		("c, case", "XML file that defines the case",
			cxxopts::value<std::string>())
		("p, output-period", "interval in seconds between saved snapshots in binary output",
			cxxopts::value<float>())
		("max-run-time", "Maximum time in minutes the app can run before termination",
			cxxopts::value<float>())
		("snapshot", "Filename for snapshot of final step (default \"snapshot.bin\")",
			cxxopts::value<std::string>());
	
	auto const_argv = const_cast<const char**>(argv);       // Workaround for cxxopts bug!!!
	const auto result = cxx_options.parse(argc, const_argv);

	if (result.count("case") == 0)
		throw std::runtime_error(std::string("No option specified for case"));
	else
		run_options.case_filename = result["case"].as<std::string>();

	if (result.count("time") == 0)
		throw std::runtime_error(std::string("No options specified for time"));
	else
		run_options.time = result["time"].as<float>();

	if (result.count("output-period") == 0)
		run_options.output_period = 0;
	else
		run_options.output_period = result["output-period"].as<float>();

	if (result.count("particles-output") != 0)
		run_options.particles_output_filename = result["particles-output"].as<std::string>();

	if (result.count("surface-output") != 0)
		run_options.surface_output_filename = result["surface-output"].as<std::string>();

	if (!run_options.particles_output_filename && !run_options.surface_output_filename)
		throw std::runtime_error(std::string("No option specified for output"));

	if (result.count("max-run-time") == 1)
		run_options.max_run_time = int(result["max-run-time"].as<float>() * 60.0f);

	if (result.count("input") == 1)
		run_options.input = result["input"].as<std::string>();

	return run_options;
}

int main(int argc, char **argv) {
	const auto options = get_options(argc, argv);

	const auto case_def = get_case_from_XML(options.case_filename.c_str());

	ParticleSystem ps(case_def);

	{
		SaveVTK save_VTK(ps.get_boundary_begin(), ps.get_boundary_end());
		if (options.particles_output_filename)
			save_VTK.save_particles(*options.particles_output_filename + "-boundary");

		if (options.surface_output_filename)
			save_VTK.save_surface(*options.surface_output_filename + "-boundary", case_def.h);
	}
	
	bool user_exit = false;
	auto& now = std::chrono::steady_clock::now;
	const auto start_time = now();
	int output_step = 0;
	while (
		ps.current_time() < options.time &&
		!user_exit &&
		(!options.max_run_time || now() - start_time < std::chrono::seconds{ *(options.max_run_time) }))
	{
		while (output_step * options.output_period <= ps.current_time()) {
			SaveVTK save_VTK(ps.get_fluid_begin(), ps.get_fluid_end());
			if (options.particles_output_filename)
				save_VTK.save_particles(*options.particles_output_filename + "-fluid" + std::to_string(output_step));

			if (options.surface_output_filename)
				save_VTK.save_surface(
					*options.surface_output_filename + "-fluid" + std::to_string(output_step), case_def.h);

			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";

			++output_step;
		}

		ps.simulation_step();

		while (_kbhit()) {
			auto key = _getch();
			if (key == 27)
				user_exit = true;
		}
	}
	std::cout << '\n';

	if (user_exit)
		std::cout << "Cancelled by user\n\n";

	std::cout << "Time of simulation :" << ps.current_time() << " s\n"
		<<"    (of requested " << options.time << " s)\n\n";

	// Save one last step, mostly for debugging
	// Consider removing this later
	SaveVTK save_VTK(ps.get_fluid_begin(), ps.get_fluid_end());
	if (options.particles_output_filename)
		save_VTK.save_particles(*options.particles_output_filename + "-fluid" + std::to_string(output_step));

	if (options.surface_output_filename)
		save_VTK.save_surface(
			*options.surface_output_filename + "-fluid" + std::to_string(output_step), case_def.h);

	std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";

	std::cout << "Verlet steps integrated : " << ps.current_step() << "\n";
	std::cout << "Steps saved in output file : " << output_step << "\n\n";

	std::cout << "Saving two snapshots of particles\n";

	save_particles_to_xml(ps.get_fluid_begin(), ps.get_fluid_end(), "prev_fluid_particles.xml", "particle");
	save_particles_to_xml(ps.get_boundary_begin(), ps.get_boundary_end(), "prev_boundary_particles.xml", "particle");
	ps.simulation_step();
	save_particles_to_xml(ps.get_fluid_begin(), ps.get_fluid_end(), "fluid_particles.xml", "particle");
	save_particles_to_xml(ps.get_boundary_begin(), ps.get_boundary_end(), "boundary_particles.xml", "particle");

	if (ps.get_springs_begin() != ps.get_springs_end()) {
		std::cout << "Saving mass-spring system\n";

		save_springs_to_xml(ps.get_springs_begin(), ps.get_springs_end(), "springs.xml", "spring");
	}

	return 0;
}