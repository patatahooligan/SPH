#pragma once

#include <cmath>

// Constants needed in more than one translation unit
constexpr float pi = 3.14159265359f;

// Size of system
constexpr float
	size = 20.0f,
	sizex = 20.0f,
	sizey = 20.0f,
	sizez = 20.0f,
	offsetz = 20.0f;

// Display window (and video output) dimensions
constexpr int output_width = 800;
constexpr int output_height = 800;

// Simulation parameters
constexpr size_t num_of_particles = 100;
constexpr float smoothing_length = 1.0f;
constexpr float kernel_constant = 3.0f/359.0f * pi * smoothing_length * smoothing_length * smoothing_length;
constexpr float particle_mass = 1.0f;
constexpr float visc_a = 1.0f;		// Bulk viscocity
constexpr float visc_b = 1.0f;		// von Neumann-Ritchmyer viscocity
constexpr float speed_of_sound = 140.0f;
constexpr float gravity_constant = 10.0f;
constexpr float thermal_diffusion_constant = 0.5910f;		// Thermal conductivity of water
constexpr float reference_density = 1.0f;
constexpr float jump_number = 10.0f;
constexpr float xsph_coeff = 0.2f;