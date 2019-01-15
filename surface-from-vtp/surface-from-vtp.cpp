#include "stdafx.h"

#include "fileIO.h"

int main(int argc, char* argv[]) {
	if (argc < 3)
		return EXIT_FAILURE;

	constexpr float resolution = 0.006;
	const std::string
		input_prefix = argv[1],
		output_prefix = argv[2];

	const int	begin_index = argc >= 4 ? std::stoi(argv[3]) : 0;
	int end_index = begin_index;
	while (std::filesystem::exists(input_prefix + std::to_string(end_index) + ".vtp"))
		++end_index;

	for (int i = begin_index; i < end_index; ++i) {
		const auto polydata = polydata_from_vtp(input_prefix + std::to_string(i) + ".vtp");
		const auto surface = surface_from_polydata(polydata, resolution);
		const auto output_filename = output_prefix + std::to_string(i);
		save_surface_to_vtp(surface, output_filename);

		std::cout << "Saved output to " << output_filename << '\n';
	}
}