// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "boost/optional.hpp"

#include "render.h"
#include "physics.h"
#include "savestate.h"
#include "XMLReader.h"

// glut is designed to take pointer to void() for its callback type. This means that we can't
// pass a local variable to our callback, nor use a lambda that captures a local variable.
// A local variable enforces proper storage and lifetime duration, while a global pointer
// can be used inside a capture-less lambda.
static ParticleSystem *ps_pointer = nullptr;
static SaveState *save_state_pointer = nullptr;

void render_func() {
	const auto particles_begin = ps_pointer->get_particlearray().begin();
	render_particles(particles_begin, particles_begin + ps_pointer->get_num_of_fluid_particles());
}

void idle_func() {
	static size_t current_frame = 0;
	ps_pointer->simulation_step();
	if (ps_pointer->current_time() >= current_frame / framerate) {
		glutPostRedisplay();
		if (save_state_pointer)
			save_state_pointer->save(ps_pointer->get_particlearray());
		++current_frame;
	}
}

int main(int argc, char **argv) {
	omp_set_num_threads(5);

	// TODO: full command line argument parsing
	// For now 1st argument is case XML, 2nd is output file
	const auto case_def = get_case_from_XML(argv[1]);
	ParticleSystem ps(case_def, cubic_spline, cubic_spline_gradient);
	ps_pointer = &ps;

	render_init(&argc, argv, render_func, idle_func,
		case_def.particles.point_min, case_def.particles.point_max, -case_def.gravity);

	boost::optional<SaveState> save_state;
	if (argc == 3) {
		save_state.emplace(argv[2], ps.get_num_of_particles());
		save_state_pointer = save_state.get_ptr();
	}

	std::cout << glGetString(GL_VERSION) << std::endl;

	// Enter main loop
	glutMainLoop();

    return 0;
}