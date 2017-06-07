// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "GL\freeglut.h"

#include "render.h"
#include "physics.h"
#include "video.h"


ParticleSystem ps;
Video video_output; 

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
	video_output.video_finalize();

    return 0;
}