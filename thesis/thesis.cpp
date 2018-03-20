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
	render_particles(ps_pointer->get_particlearray());
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

boost::optional<SaveState> parse_arguments(int argc, char **argv) {
	if (argc == 1)
		return {};
	else if (argc == 2)
		return SaveState{ argv[1] };
	else
		throw std::runtime_error("Too many arguments");
}

int main(int argc, char **argv) {
	omp_set_num_threads(5);

	render_init(&argc, argv, render_func, idle_func);

	// Because glut might parse some of the arguments, we have to parse our own
	// after render_init

	// TODO: full command line argument parsing
	// For now 1st argument is case XML, 2nd is output file

	ParticleSystem ps(get_case_from_XML(argv[1]));
	ps_pointer = &ps;

	boost::optional<SaveState> save_state;
	if (argc == 2) {
		save_state.emplace(argv[2]);
		save_state_pointer = save_state.get_ptr();
	}
		
	ps.randomize_particles();

	GLubyte const *s = glGetString(GL_VERSION);
	std::cout << s << std::endl;

	// Enter main loop
	glutMainLoop();

    return 0;
}