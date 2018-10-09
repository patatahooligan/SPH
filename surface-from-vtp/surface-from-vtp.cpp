#include "stdafx.h"

#include "fileIO.h"

int main(int argc, char* argv[]) {
	if (argc != 3)
		return EXIT_FAILURE;

	constexpr float h = 0.0346410161513775f;

	int i = 0;
	const std::string
		input_prefix = argv[1],
		output_prefix = argv[2];
	std::string input_filename = input_prefix + std::to_string(i) + ".vtp";
	while (std::filesystem::exists(input_filename)) {
		const auto polydata = polydata_from_vtp(input_filename);
		const auto surface = surface_from_polydata(polydata, h);
		save_surface_to_vtp(surface, output_prefix + std::to_string(i));

		input_filename = input_prefix + std::to_string(++i) + ".vtp";
	}
}