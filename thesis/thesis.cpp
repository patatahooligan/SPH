// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "GL\freeglut.h"

#include "render.h"
#include "physics.h"


ParticleSystem ps;

int main(int argc, char **argv)
{
	unsigned int _clearfp();
	//unsigned int current_control;
	//unsigned int err_no = _controlfp_s(&current_control, ~(_EM_ZERODIVIDE|_EM_OVERFLOW), _MCW_EM);

	render_init(argc, argv);
	
	ps.randomize_particles();
	ps.calculate_initial_conditions();

	GLubyte const *s = glGetString(GL_VERSION);
	std::cout << s << std::endl;

	// Enter main loop
	glutMainLoop();

    return 0;
}