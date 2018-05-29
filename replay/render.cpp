#include "stdafx.h"

// Disable unsafe parameter error for ublas::matrix
#define _SCL_SECURE_NO_WARNINGS

#include "particle.h"
#include "render.h"
#include "constants.h"


constexpr double particle_display_size = 0.01f;

Renderer::Renderer(int &argc, char **argv,
	LoadState load_state, const float time_step,
	const Vec3f &point_min, const Vec3f &point_max) :
	load_state(std::move(load_state)),
	time_step(time_step),
	point_min(point_min), point_max(point_max)
{
	// Initialize glut and create the window
	glutInit(&argc, argv);
	glutInitWindowPosition(-1, -1);
	glutInitWindowSize(output_width, output_height);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("SPH");

	// Pre-load file to memory
	auto snapshot = this->load_state.load(LoadState::Mode::Position);
	while (snapshot) {
		snapshots.emplace_back(std::move(*snapshot));
		snapshot = this->load_state.load(LoadState::Mode::Position);
	}

	// Register callbacks
	global_renderer_pointer = this;
	glutDisplayFunc([]() {
		global_renderer_pointer->render_fluid_particles();
		});

	glutIdleFunc([]() {
		global_renderer_pointer->idle_func();
		});

	glutKeyboardFunc([](const unsigned char key, const int x, const int y) {
		global_renderer_pointer->keyboard_func(key, x, y);
		});

	glutMouseFunc([](const int button, const int state, const int x, const int y) {
		global_renderer_pointer->mouse_func(button, state, x, y);
		});

	glutMotionFunc([](const int x, const int y) {
		global_renderer_pointer->mouse_motion_func(x, y);
		});

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

	update_camera_position();
}

void Renderer::update_camera_position() {
	// Use an offset such that a particle on the edge of the simulation space is not
	// clipped off

	// Note that this will break for extreme field of view values
	constexpr float
		fov_degrees = 90.0f,
		fov_radians = (fov_degrees * pi) / 180.0f,
		offset = (float(particle_display_size) / 2.0f) * (1.0f + std::numeric_limits<float>::round_error());

	constexpr Vec3f up_vector = Vec3f::z_unit();

	const Vec3f
		center = (point_max + point_min) / 2,
		size = point_max - point_min;

	const float
		max_dimension = Vec3f{ size.x, size.y, 0.0f }.length_squared(),
		camera_distance_from_center = (max_dimension + offset) / std::tan(fov_radians / 2),
		near_clipping_distance = std::max(camera_distance_from_center - max_dimension / 2.0f - offset, 0.0f),
		far_clipping_distance = camera_distance_from_center + max_dimension / 2.0f + offset;
	
	const Vec3f
		camera_relative_position =
			Vec3f{ -std::cos(camera_angle), std::sin(camera_angle), 0.0f } * camera_distance_from_center,
		camera_position = center + camera_distance_coef * camera_relative_position;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov_degrees, 1.0f, near_clipping_distance, far_clipping_distance);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(
		camera_position.x, camera_position.y, camera_position.z,
		center.x, center.y, center.z,
		up_vector.x, up_vector.y, up_vector.z
	);
}

void Renderer::render_fluid_particles() {
	// Draws the particles of ps as small spheres.

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	for (auto particle_iterator = snapshots[step].begin();
		particle_iterator!= snapshots[step].begin() + load_state.get_num_of_fluid_particles();
		++particle_iterator)
	{
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

void Renderer::idle_func() {
	if (!play) return;

	const auto elapsed_time_ms =
		std::chrono::duration_cast<std::chrono::milliseconds>
		(clock::now() - playback_starting_time).count();

	const auto new_step = int(elapsed_time_ms / (1000.0f * time_step));

	if (new_step != step) {
		if (new_step >= snapshots.size()) {
			step = snapshots.size() - 1;
			play = false;
		}
		else
			step = new_step;

		glutPostRedisplay();
	}
}

void Renderer::keyboard_func(const unsigned char key, const int x, const int y) {
	switch (key) {
	case ' ':
		if (play) {
			step = 0;
			play = false;
		}
		else {
			play = true;
			playback_starting_time = clock::now();
		}
		break;
	case 27:
		glutLeaveMainLoop();
		break;
	}
}

void Renderer::mouse_func(const int button, const int state, const int x, const int y) {
	constexpr float zoom_factor = 0.05f;

	// Constants that glut doesn't name
	enum class Mouse : int {
		ScrollUp = 3,
		ScrollDown = 4
	};

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		previous_mouse_x = x;
		previous_mouse_y = y;
	}
	else if (button == int(Mouse::ScrollUp)) {
		camera_distance_coef /= (1.0f + zoom_factor);
	}
	else if (button == int(Mouse::ScrollDown)) {
		camera_distance_coef *= (1.0f + zoom_factor);
	}
	update_camera_position();
}

void Renderer::mouse_motion_func(const int x, const int y) {
	constexpr float sensitivity = 0.01f;

	camera_angle += sensitivity * (x - previous_mouse_x);
	camera_angle = std::fmod(camera_angle, 2 * pi);
	previous_mouse_x = x;

	update_camera_position();
	glutPostRedisplay();
}