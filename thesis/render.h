#pragma once

#include "physics.h"

void render_init(int *argc, char **argv);

void render_particles(const ParticleSystem &ps);

void render();

void keyboardfunc(unsigned char key, int x, int y);