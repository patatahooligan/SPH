#pragma once

#include "particle.h"

using GlutCallbackType = void();

void render_init(
	int *argc, char **argv, GlutCallbackType *render_function = nullptr,
	GlutCallbackType *idle_callback = nullptr);

void render_particles(const particlearray &particles);

void keyboardfunc(unsigned char key, int x, int y);