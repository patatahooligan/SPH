// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "physics.h"
#include "savestate.h"
#include "XMLReader.h"
#include "loadstate.h"

struct RunOptions {
	std::string case_filename;
	float time, output_period;
	SaveState::WriteMode write_mode;

	std::optional<int> max_run_time;
	std::optional<std::string> input;
	std::optional<std::string> binary_output_filename, vtk_output_filename;
};

auto get_options(int argc, char **argv) {
	cxxopts::Options cxx_options("SPH", "Smoothed Particle Hydrodynamics");
	RunOptions run_options;

	cxx_options.add_options()
		("b, binary-output", "Binary file to hold the result of the simulation",
			cxxopts::value<std::string>())
		("vtk-output", "Prefix for group of vtk files to hold the result of the simulation",
			cxxopts::value<std::string>())
		("overwrite", "Overwrite the binary output file if it exists (ignored if --append)")
		("append", "Append to binary output file if it exists (disables --overwrite)")
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

	// Mandatory options
	if (result.count("append") > 0)
		run_options.write_mode = SaveState::WriteMode::Append;
	else if (result.count("overwrite") > 0)
		run_options.write_mode = SaveState::WriteMode::Overwrite;
	else
		run_options.write_mode = SaveState::WriteMode::DontModify;

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

	// Optional
	if (result.count("binary-output") != 0)
		run_options.binary_output_filename = result["binary-output"].as<std::string>();

	if (result.count("vtk-output") != 0)
		run_options.vtk_output_filename = result["vtk-output"].as<std::string>();

	if (!run_options.binary_output_filename && !run_options.vtk_output_filename)
		throw std::runtime_error(std::string("No option specified for output"));

	if (result.count("max-run-time") == 1)
		run_options.max_run_time = int(result["max-run-time"].as<float>() * 60.0f);

	if (result.count("input") == 1)
		run_options.input = result["input"].as<std::string>();

	return run_options;
}

int main(int argc, char **argv) {
	omp_set_num_threads(5);

	const auto options = get_options(argc, argv);

	const auto case_def = get_case_from_XML(options.case_filename);

	ParticleSystem ps(case_def);

	save_VTK(ps.get_boundary_begin(), ps.get_boundary_end(),
		*options.vtk_output_filename + "-boundary");

	std::optional<SaveBinary> save_binary;
	if (options.binary_output_filename) {
		save_binary.emplace(
			*options.binary_output_filename, options.write_mode,
			ps.get_num_of_fluid_particles(), ps.get_num_of_particles());
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
		ps.simulation_step();

		while (output_step * options.output_period < ps.current_time()) {
			if (save_binary)
				save_binary->save(ps.get_fluid_begin(), ps.get_fluid_end());

			if (options.vtk_output_filename) {
				save_VTK(ps.get_fluid_begin(), ps.get_fluid_end(),
					*options.vtk_output_filename + "-fluid-" + std::to_string(output_step));
			}

			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";

			++output_step;
		}

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

	std::cout << "Verlet steps integrated : " << ps.current_step() << "\n";
	std::cout << "Steps saved in output file : " << output_step << '\n\n';

	return 0;
}