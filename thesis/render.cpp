#include "stdafx.h"

// Disable unsafe parameter error for ublas::matrix
#define _SCL_SECURE_NO_WARNINGS

#include "physics.h"
#include "render.h"
#include "constants.h"


constexpr double particle_display_size = 1.0;
constexpr double fov = 60.0;


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

	// Define the projection of the world onto the camera
	// Arguments are fov(degrees), aspect ratio, near clipping plane, far clipping plane.
	// To ensure that no polygon is clipped, the clipping planes extend outside of the walls
	// of the container by a bit more than the particle display size.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov, 1.0, offsetz - 1.1 * particle_display_size, offsetz + size + 1.1 * particle_display_size);

	// Define the camera position and orientation
	// The arguments are as follows
	// x, y, z of camera position
	// x, y, z of point the camera is lookit towards
	// x, y, z of camera's up vector
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(size / 2, size / 2, -offsetz, size / 2, size / 2, size / 2, 0.0, 1.0, 0.0);
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
	glutSwapBuffers();
}

void keyboardfunc(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		glutLeaveMainLoop();
		break;
	}
}