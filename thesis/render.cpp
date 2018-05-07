#include "stdafx.h"

// Disable unsafe parameter error for ublas::matrix
#define _SCL_SECURE_NO_WARNINGS

#include "physics.h"
#include "render.h"
#include "constants.h"


constexpr double particle_display_size = 0.01f;

void render_init(int *argc, char **argv,
	GlutCallbackType* const render_function, GlutCallbackType* const idle_callback,
	const Vec3f &point_min, const Vec3f &point_max, const Vec3f &up_vector)
{
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

	set_camera(point_min, point_max, up_vector);
}

void set_camera(const Vec3f &point_min, const Vec3f &point_max, const Vec3f &up_vector) {
	const Vec3f
		center = (point_max + point_min) / 2,
		size = point_max - point_min,
		camera_relative_position = -dot_product(size, { 1.0f, 0.0f, 0.0f }),
		camera_position = center - camera_relative_position;

	// Use an offset such that a particle on the edge of the simulation space is not
	// clipped off
	constexpr float
		fov = 90.0f,
		offset = (float(particle_display_size) / 2.0f) * (1.0f + std::numeric_limits<float>::round_error());

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov, 1.0f, camera_relative_position.length() / 2.0f - offset, 
		camera_relative_position.length() * 1.5f + offset);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(
		camera_position.x, camera_position.y, camera_position.z,
		center.x, center.y, center.z,
		up_vector.x, up_vector.y, up_vector.z
	);
}

void render_particles(
	const ParticleContainer::const_iterator begin, const ParticleContainer::const_iterator end)
{
	// Draws the particles of ps as small spheres.

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	for (auto particle_iterator = begin; particle_iterator != end; ++particle_iterator) {
		glPushMatrix();
		glTranslatef(
			particle_iterator->position.x,
			particle_iterator->position.y,
			particle_iterator->position.z);
		glutSolidSphere(particle_display_size, 10, 10);
		glPopMatrix();
	}
	glutSwapBuffers();
}

void keyboardfunc(const unsigned char key, const int x, const int y) {
	switch (key) {
	case 27:
		glutLeaveMainLoop();
		break;
	}
}