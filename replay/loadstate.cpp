#include "stdafx.h"

#include <iomanip>

#include "loadstate.h"

std::istream& operator>>(std::istream &in, Vec3f vec) {
	in >> vec.x >> vec.y >> vec.z;
	return in;
}

std::istream& operator>>(std::istream &in, Particle &particle) {
	in >>
		particle.position >>
		particle.velocity >>
		particle.acceleration >>
		particle.temperature >>
		particle.temperature_derivative >>
		particle.density >>
		particle.density_derivative >>
		particle.temperature >>
		particle.temperature_derivative >>
		particle.viscocity;
	return in;
}

LoadState::LoadState(std::string input_filename) :
	input_file(input_filename, std::ios_base::binary)
{
	if (!input_file.is_open())
		throw std::runtime_error("Cannot open file");

	// Verify the format of the file
	unsigned int version;
	std::string header;
	const std::string target("SaveState");

	input_file >> std::setw(target.length()) >> header;
	input_file >> version;
	if (header.compare(target) || version != 1)
		throw std::runtime_error("Invalid input file format");
}

void LoadState::load(particlearray& target_array) {
	for (auto& particle : target_array) {
		input_file >> particle;
	}
}