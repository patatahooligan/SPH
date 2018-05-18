#pragma once

#include "particle.h"
#include "loadstate.h"

using GlutCallbackType = void();

class Renderer {
	private:
		inline static Renderer* global_renderer_pointer = nullptr;
		LoadState load_state;
		const Vec3f point_min, point_max;
		const float time_step;
		float camera_angle = 0.0f;
		int previous_mouse_x, previous_mouse_y;
		ParticleContainer particles;

		void render_fluid_particles();

		void idle_func();

		void update_camera_position();

		void keyboard_func(unsigned char key, int x, int y);

		void mouse_func(int button, int state, int x, int y);

		void mouse_motion_func(int x, int y);

	public:
		Renderer(int &argc, char **argv,
			LoadState load_state, float time_step,
			const Vec3f &point_min, const Vec3f &point_max);

		void main_loop() {
			glutMainLoop();
		}
};