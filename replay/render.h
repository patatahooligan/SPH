#pragma once

#include "particle.h"

using GlutCallbackType = void();

void render_init(int *argc, char **argv,
	GlutCallbackType* const render_function, GlutCallbackType* const idle_callback,
	const Vec3f &point_min, const Vec3f &point_max, const Vec3f &up_vector);

void render_particles(ParticleContainer::const_iterator begin, ParticleContainer::const_iterator end);

void set_camera(const Vec3f &point_min, const Vec3f &point_max, const Vec3f &up_vector);

void keyboardfunc(unsigned char key, int x, int y);