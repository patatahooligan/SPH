// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "GL\freeglut.h"

#include "render.h"
#include "physics.h"
#include <omp.h>  

// glut is designed to take pointer to void() for its callback type. This means that we can't
// pass a local variable to our callback, nor use a lambda that captures a local variable.
// A local variable enforces proper storage and lifetime duration, while a global pointer
// can be used inside a capture-less lambda.
static ParticleSystem *ps_pointer;

void render_func() {
	render_particles(ps_pointer->get_particlearray());
}

void idle_func() {
	static size_t current_frame = 0;
	ps_pointer->simulation_step();
	if (ps_pointer->current_time() >= current_frame / framerate) {
		glutPostRedisplay();
		++current_frame;
	}
}

int main(int argc, char **argv) {
	ParticleSystem ps;
	ps_pointer = &ps;
	omp_set_num_threads(5);

	render_init(&argc, argv, render_func, idle_func);
	
	ps.randomize_particles();
	ps.calculate_initial_conditions();

	GLubyte const *s = glGetString(GL_VERSION);
	std::cout << s << std::endl;

	// Enter main loop
	glutMainLoop();

    return 0;
}