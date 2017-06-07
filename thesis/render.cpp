#include "stdafx.h"

// Disable unsafe parameter error for ublas::matrix
#define _SCL_SECURE_NO_WARNINGS

#include "physics.h"
#include "render.h"
#include "constants.h"
#include "video.h"

extern ParticleSystem ps;
extern Video video_output;


// Render parameters
void (* const render_function)(const ParticleSystem&) = render_particles;

const double particle_display_size = 1.0;


void render_init(int argc, char **argv) {
	// Initialize glut and create the window
	glutInit(&argc, argv);
	glutInitWindowPosition(-1, -1);
	glutInitWindowSize(output_width, output_height);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("SPH");

	// Register callback
	glutDisplayFunc(render);
	glutIdleFunc([](){ps.simulation_step();});
	glutKeyboardFunc(keyboardfunc);

	// Define a perspective projection for the camera
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// If x,y are not symmetric, perspective is skewed so use [-size/2,size/2]
	// instead of [0, size]. Offsetz is important because a clipping plane that
	// is too near also causes weird perspective.
	glFrustum(-size/2, size/2, -size/2, size/2, offsetz, offsetz+size);
	
	// Reposition camera to move back to [0, size].
	glTranslated(-size/2, -size/2, -offsetz);

	video_output.video_init();
}

void render_sphere(ParticleSystem ps) {
	// Only exists to test glut functionality

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(25.0f, 25.0f, -25.0f);
	glutSolidSphere(particle_display_size, 100, 100);
	glTranslatef(50.0f, 50.0f, -50.0f);
	glutSolidSphere(particle_display_size, 100, 100);
	glPopMatrix();

	glutSwapBuffers();
}

void render_particles(const ParticleSystem &ps) {
	// Draws the particles of ps as small spheres.

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	for (unsigned i = 0; i < num_of_particles; i++) {
		glPushMatrix();
		glTranslatef(
			ps.particles[i].position.x,
			ps.particles[i].position.y,
			-ps.particles[i].position.z);
		glutSolidSphere(particle_display_size, 10, 10);
		glPopMatrix();
	}
}

void render() {
	// Function to be used as glut display function. Calls a
	// function given by global compile-time constant render_function.
	// render_function takes a ParticleSystem as the only argument.

	render_function(ps);
	video_output.encode_frame(ps.current_time());
	glutSwapBuffers();
}

void keyboardfunc(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		glutLeaveMainLoop();
		break;
	}
}