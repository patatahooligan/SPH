// replay.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "loadstate.h"
#include "render.h"
#include "XMLReader.h"


struct RunOptions {
	std::string case_filename, binary_input_filename;
	float time_step;
};

auto get_options(int argc, char **argv) {
	cxxopts::Options cxx_options("SPH", "Smoothed Particle Hydrodynamics");
	RunOptions run_options;

	cxx_options.add_options()
		("i, input", "Binary file that holds the result of a simulation",
		cxxopts::value<std::string>())
		("c, case", "XML file that defines the case",
		cxxopts::value<std::string>())
		("t, time-step", "interval in seconds between saved snapshots in binary input",
		cxxopts::value<float>());

	auto const_argv = const_cast<const char**>(argv);       // Workaround for cxxopts bug!!!
	const auto result = cxx_options.parse(argc, const_argv);

	// Mandatory options
	if (result.count("input") == 0)
		throw std::runtime_error(std::string("No option specified for output"));
	else
		run_options.binary_input_filename = result["input"].as<std::string>();

	if (result.count("case") == 0)
		throw std::runtime_error(std::string("No option specified for case"));
	else
		run_options.case_filename = result["case"].as<std::string>();

	if (result.count("time-step") == 0)
		throw std::runtime_error(std::string("No options specified for time step"));
	else
		run_options.time_step = result["time-step"].as<float>();

	return run_options;
}

int main(int argc, char **argv) {
	const auto options = get_options(argc, argv);
	const auto case_def = get_case_from_XML(options.case_filename);

	Renderer renderer{
		argc, argv, LoadState{options.binary_input_filename}, options.time_step,
		case_def.particles.point_min, case_def.particles.point_max };

	// Enter main loop
	renderer.main_loop();

    return 0;
}