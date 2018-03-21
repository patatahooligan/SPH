// replay.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "loadstate.h"
#include "render.h"

// glut is designed to take pointer to void() for its callback type. This means that we can't
// pass a local variable to our callback, nor use a lambda that captures a local variable.
// A local variable enforces proper storage and lifetime duration, while a global pointer
// can be used inside a capture-less lambda.
static LoadState *load_state_pointer = nullptr;
static ParticleContainer *particle_array_pointer = nullptr;

void render_func() {
	assert(particle_array_pointer);
	render_particles(*particle_array_pointer);
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

int main(int argc, char **argv) {
	if (argc != 2)
		throw std::runtime_error("Wrong number of arguments");

	render_init(&argc, argv, render_func, idle_func);

	// Because glut might parse some of the arguments, we have to parse our own
	// after render_init
	LoadState load_state(argv[1]);
	ParticleContainer particle_array(num_of_particles);
	load_state_pointer = &load_state;
	particle_array_pointer = &particle_array;

	// Enter main loop
	glutMainLoop();

    return 0;
}