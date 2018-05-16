// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "physics.h"
#include "savestate.h"
#include "XMLReader.h"
#include "loadstate.h"

struct RunOptions {
	std::string output_filename, case_filename, snapshot_filename;
	float time, output_period;

	std::optional<int> max_run_time;
	std::optional<std::string> input;
};

auto get_options(int argc, char **argv) {
	cxxopts::Options cxx_options("SPH", "Smoothed Particle Hydrodynamics");
	RunOptions run_options;

	cxx_options.add_options()
		("o, output", "Binary file to hold the result of the simulation",
			cxxopts::value<std::string>())
		("t, time", "Time of simulation run in seconds",
			cxxopts::value<float>())
		("c, case", "XML file that defines the case",
			cxxopts::value<std::string>())
		("p, output-period", "interval in seconds between saved snapshots in binary output",
			cxxopts::value<float>())
		("max-run-time", "Maximum time in minutes the app can run before termination",
			cxxopts::value<float>())
		("i, input", "Binary snapshot file of previous run",
			cxxopts::value<std::string>())
		("snapshot", "Filename for snapshot of final step (default \"snapshot.bin\")",
			cxxopts::value<std::string>());
	
	auto const_argv = const_cast<const char**>(argv);       // Workaround for cxxopts bug!!!
	const auto result = cxx_options.parse(argc, const_argv);

	// Mandatory options
	if (result.count("output") == 0)
		throw std::runtime_error(std::string("No option specified for output"));
	else
		run_options.output_filename = result["output"].as<std::string>();

	if (result.count("case") == 0)
		throw std::runtime_error(std::string("No option specified for case"));
	else
		run_options.case_filename = result["case"].as<std::string>();

	if (result.count("time") == 0)
		throw std::runtime_error(std::string("No options specified for time"));
	else
		run_options.time = result["time"].as<float>();

	if (result.count("output-period") == 0)
		throw std::runtime_error(std::string("No options specified for output-period"));
	else
		run_options.output_period = result["output-period"].as<float>();

	// Optional
	if (result.count("snapshot") == 0)
		run_options.snapshot_filename = "snapshot.bin";
	else
		run_options.snapshot_filename = result["snapshot"].as<std::string>();

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

	auto ps = [&]() {
		if (options.input) {
			LoadState load_state{ *options.input };
			ParticleContainer previous_particles, current_particles;

			load_state.load(previous_particles, LoadState::Mode::Full);
			load_state.load(current_particles, LoadState::Mode::Full);

			return ParticleSystem{
				case_def, previous_particles, current_particles,
				load_state.get_num_of_fluid_particles() };
		}
		else
			return ParticleSystem{ case_def };
	} ();

	SaveState save_state(options.output_filename, ps.get_num_of_fluid_particles(), ps.get_num_of_particles());

	bool user_exit = false;
	auto& now = std::chrono::steady_clock::now;
	const auto start_time = now();
	while (
		ps.current_time() < options.time &&
		!user_exit &&
		options.max_run_time && now() - start_time < std::chrono::seconds{ *(options.max_run_time) })
	{
		ps.simulation_step();

		if (save_state.get_step() * options.output_period < ps.current_time()) {
			save_state.save(ps.get_particlearray(), SaveState::Mode::Position);
			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";
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

	std::cout << "Steps saved in binary file : " << save_state.get_step() << '\n\n';

	std::cout << "Saving snapshot of particles at moment of termination\n"
			<< "Saved in order (previous step, current step)\n\n";

	SaveState(options.snapshot_filename, ps.get_num_of_fluid_particles(), ps.get_num_of_particles())
		.save(ps.get_previous_particlearray(), SaveState::Mode::Full)
		.save(ps.get_particlearray(), SaveState::Mode::Full);

	return 0;
}