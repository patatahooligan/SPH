// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "GL\freeglut.h"

#include "render.h"
#include "physics.h"
#include <omp.h>  


ParticleSystem ps;

int main(int argc, char **argv) {
	omp_set_num_threads(5);

	render_init(&argc, argv);
	
	ps.randomize_particles();
	ps.calculate_initial_conditions();

	GLubyte const *s = glGetString(GL_VERSION);
	std::cout << s << std::endl;

	// Enter main loop
	glutMainLoop();

    return 0;
}