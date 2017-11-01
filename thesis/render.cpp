#include "stdafx.h"

// Disable unsafe parameter error for ublas::matrix
#define _SCL_SECURE_NO_WARNINGS

#include "physics.h"
#include "render.h"
#include "constants.h"


constexpr double particle_display_size = 1.0;


void render_init(int *argc, char **argv, GlutCallbackType *render_function, GlutCallbackType *idle_callback) {
	// Initialize glut and create the window
	glutInit(argc, argv);
	glutInitWindowPosition(-1, -1);
	glutInitWindowSize(output_width, output_height);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("SPH");

	// Register callbacks
	if (render_function)
		glutDisplayFunc(render_function);
	else
		std::clog << "WARNING: no render function registered\n";
	if (idle_callback)
		glutIdleFunc(idle_callback);
	else
		std::clog << "WARNING: no idle function registered\n";
	glutKeyboardFunc(keyboardfunc);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_CONTINUE_EXECUTION);

	// Define a perspective projection for the camera
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// If x,y are not symmetric, perspective is skewed so use [-size/2,size/2]
	// instead of [0, size]. Offsetz is important because a clipping plane that
	// is too near also causes weird perspective.
	glFrustum(-size/2, size/2, -size/2, size/2, offsetz, offsetz+size);
	
	// Reposition camera to move back to [0, size].
	glTranslated(-size/2, -size/2, -offsetz);
}

void render_particles(const particlearray &particles) {
	// Draws the particles of ps as small spheres.

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	for (auto &particle : particles) {
		glPushMatrix();
		glTranslatef(
			particle.position.x,
			particle.position.y,
			particle.position.z);
		glutSolidSphere(particle_display_size, 10, 10);
		glPopMatrix();
	}
}

void keyboardfunc(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		glutLeaveMainLoop();
		break;
	}
}