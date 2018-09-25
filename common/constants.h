#pragma once

#include "vec3f.h"

// Constants needed in more than one translation unit
constexpr float pi = 3.14159265359f;

// Display window (and video output) dimensions
constexpr int output_width = 800;
constexpr int output_height = 800;
constexpr float framerate = 60.0f;

enum class ParticleType {
	Void,
	Fluid,
	Boundary
};

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
	float alpha;        // Coefficient in artificial viscosity

	float friction_coef;

	struct {
		bool on = false;
		float stiffness, start_of_stiffness_change, rate_of_stiffness_change, damping, max_length;
	} spring;

	struct {
		float density, mass;
		Vec3f point_min, point_max;
	} particles;

	struct Box {
		Vec3f origin, size;

		bool contains(const Vec3f &position) const {
			return
				position.x >= origin.x && position.x <= origin.x + size.x &&
				position.y >= origin.y && position.y <= origin.y + size.y &&
				position.z >= origin.z && position.z <= origin.z + size.z;
		};
	};

	struct CaseDefBox : Box {
		using Type = ParticleType;
		Type type = Type::Void;
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
	};

	std::vector<CaseDefBox> particle_boxes;

	struct PolyDataModel {
		vtkSmartPointer<vtkPolyData> poly_data;
		Vec3f offset, rotation;
		float scale;
	};
	std::vector<PolyDataModel> poly_data_models;
};