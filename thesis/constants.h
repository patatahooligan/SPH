#pragma once

#include <cmath>

// Constants needed in more than one translation unit
constexpr float pi = 3.14159265359f;

// Size of system
constexpr float
	size = 1.0f,
	sizex = 1.0f,
	sizey = 1.0f,
	sizez = 1.0f,
	offsetz = 1.0f;

// Display window (and video output) dimensions
constexpr int output_width = 800;
constexpr int output_height = 800;
constexpr float framerate = 30.0f;

// Simulation parameters
constexpr size_t num_of_particles = 500;
constexpr float smoothing_length = 0.2f;
constexpr float kernel_constant = 3.0f/359.0f * pi * smoothing_length * smoothing_length * smoothing_length;
constexpr float particle_mass = 0.5f / num_of_particles;
constexpr float visc_a = 1.0f; // Bulk viscocity
constexpr float visc_b = 1.0f; // von Neumann-Ritchmyer viscocity
constexpr float viscocity = 0.0007978f;
constexpr float speed_of_sound = 140.0f;
constexpr float gravity_constant = 10.0f;
constexpr float reference_density = 1.0f;
constexpr float bulk_modulus = speed_of_sound / reference_density;