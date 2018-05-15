// replay.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "loadstate.h"
#include "render.h"
#include "XMLReader.h"

// glut is designed to take pointer to void() for its callback type. This means that we can't
// pass a local variable to our callback, nor use a lambda that captures a local variable.
// A local variable enforces proper storage and lifetime duration, while a global pointer
// can be used inside a capture-less lambda.
static LoadState *load_state_pointer = nullptr;
static ParticleContainer *particle_array_pointer = nullptr;

void render_func() {
	assert(particle_array_pointer);
	render_particles(particle_array_pointer->begin(), particle_array_pointer->end());
}

void idle_func() {
	assert(load_state_pointer);
	assert(particle_array_pointer);

	// Declarations to clean-up the following lines
	constexpr float ms_per_frame = 1000.0f / framerate;
	auto& now = std::chrono::high_resolution_clock::now;

	// On first execution get current time as time of last frame
	static auto time_of_last_frame = now();

	// Measure time elapsed in ms because it offers sufficient accuracy
	auto current_time = now();
	auto time_elapsed_in_ms =
		std::chrono::duration_cast<std::chrono::milliseconds>(current_time - time_of_last_frame);

	if (time_elapsed_in_ms.count() > ms_per_frame) {
		time_of_last_frame = current_time;
		load_state_pointer->load(*particle_array_pointer);
		glutPostRedisplay();
	}
}

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
		("ts, time-step", "interval in seconds between saved snapshots in binary input",
		cxxopts::value<float>());

	auto const_argv = const_cast<const char**>(argv);       // Workaround for cxxopts bug!!!
	const auto result = cxx_options.parse(argc, const_argv);

	// Mandatory options
	if (result.count("input") == 0)
		throw std::runtime_error(std::string("No option specified for output"));
	else
		run_options.binary_input_filename = result["output"].as<std::string>();

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

	render_init(&argc, argv, render_func, idle_func,
		case_def.particles.point_min, case_def.particles.point_max, -case_def.gravity);

	LoadState load_state(options.binary_input_filename);
	ParticleContainer particle_array(load_state.num_of_particles());
	load_state_pointer = &load_state;
	particle_array_pointer = &particle_array;

	// Enter main loop
	glutMainLoop();

    return 0;
}