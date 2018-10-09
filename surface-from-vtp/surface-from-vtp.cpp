#include "stdafx.h"

#include "fileIO.h"

int main(int argc, char* argv[]) {
	if (argc != 3)
		return EXIT_FAILURE;

	constexpr float h = 0.0346410161513775f;

	const auto polydata = polydata_from_vtp(argv[1]);
	const auto surface = surface_from_polydata(polydata, h);
	save_surface_to_vtp(surface, argv[2]);
}