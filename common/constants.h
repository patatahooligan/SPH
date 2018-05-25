#pragma once

#include <cmath>
#include <vector>

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
constexpr float framerate = 60.0f;

// Constants for the simulation, but intended to be read from a case file at run-time
// as they may vary from case to case.
struct CaseDef {
	Vec3f gravity;      // Gravitational acceleration
	float rhop0;        // reference density of the fluid
	float hswl;         // Maximum still water level to calculate speedofsound using coefsound
	float gamma;        // Polytropic constant for water used in the state equation
	float speedsystem;  // Maximum system speed (by default the dam-break propagation is used)
	float coefsound;    // Coefficient to multiply speedsystem
	float speedsound;   // Speed of sound to use in the simulation (by default speedofsound=coefsound*speedsystem)
	float h;            // Smoothing length
	float cflnumber;    // Coefficient to multiply dt
	float tensile_coef; // Coefficient to use in tensile correction method

	struct {
		float density, mass;
		Vec3f point_min, point_max;
	} particles;

	struct Box {
		enum class Type{
			Void,
			Fluid,
			Boundary
		} type = Type::Void;
		struct {
			bool
				solid = false,
				bottom = false,
				top = false,
				right = false,
				left = false,
				front = false,
				back = false;
		} fillmode;
		Vec3f origin, size;
	};

	std::vector<Box> particle_boxes;
};