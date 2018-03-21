#pragma once

#include <cmath>

#include "vec3f.h"

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
constexpr float smoothing_length = 0.1f;
constexpr float kernel_constant = 3.0f/359.0f * pi * smoothing_length * smoothing_length * smoothing_length;
constexpr float particle_mass = 0.001f;
constexpr float visc_a = 1.0f; // Bulk viscocity
constexpr float visc_b = 1.0f; // von Neumann-Ritchmyer viscocity
constexpr float viscocity = 0.0007978f;
constexpr float speed_of_sound = 140.0f;
constexpr float gravity_constant = 10.0f;
constexpr float reference_density = 1.0f;
constexpr float bulk_modulus = speed_of_sound / reference_density;

struct CaseDef {
	Vec3f gravity;      // Gravitational acceleration
	float rhop0;        // reference density of the fluid
	float hswl;         // Maximum still water level to calculate speedofsound using coefsound
	float gamma;        // Polytropic constant for water used in the state equation
	float speedsystem;  // Maximum system speed (by default the dam-break propagation is used)
	float coefsound;    // Coefficient to multiply speedsystem
	float speedsound;   // Speed of sound to use in the simulation (by default speedofsound=coefsound*speedsystem)
	float coefh;        // Coefficient to calculate the smoothing length (h=coefh*sqrt(3*dp^2) in 3D)
	float cflnumber;    // Coefficient to multiply dt

	struct {
		float density;
		Vec3f point_min, point_max;
	} particles;

	struct Box {
		Vec3f origin, size;
	};

	std::vector<Box> fluid_boxes, bound_boxes;
};