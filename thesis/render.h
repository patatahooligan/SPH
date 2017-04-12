#pragma once

#include "physics.h"

void render_init(int argc, char **argv);

void render_sphere(ParticleSystem ps);

void render_particles(const ParticleSystem &ps);

void render();