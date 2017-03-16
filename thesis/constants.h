#pragma once

// Constants needed in more than one translation unit
const float pi = 3.14159265359f;

// Size of system
const float
	size = 20.0f,
	sizex = 20.0f,
	sizey = 20.0f,
	sizez = 20.0f,
	offsetz = 20.0f;

const unsigned int num_of_particles = 100;
const float smoothing_length = 1.0f;
const float kernel_constant = 3.0f/359.0f * pi * pow(smoothing_length, 3);
const float particle_mass = 1.0f;
const float visc_a = 1.0f;		// Bulk viscocity
const float visc_b = 1.0f;		// von Neumann-Ritchmyer viscocity
const float speed_of_sound = 140.0f;
const float gravity_constant = 10.0f;
const float thermal_diffusion_constant = 0.5910f;		// Thermal conductivity of water
const float reference_density = 1.0f;
const float jump_number = 10.0f;
const float xsph_coeff = 0.2f;