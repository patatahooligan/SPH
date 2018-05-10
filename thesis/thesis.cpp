// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "physics.h"
#include "savestate.h"
#include "XMLReader.h"


int main(int argc, char **argv) {
	omp_set_num_threads(5);

	// TODO: full command line argument parsing
	// For now 1st argument is case XML, 2nd is output file
	const auto case_def = get_case_from_XML(argv[1]);
	ParticleSystem ps(case_def);

	SaveState save_state(argv[2], ps.get_num_of_particles());

	constexpr float target_time = 2.0f;
	while (ps.current_time() < target_time) {
		ps.simulation_step();
		if (save_state.get_step() < ps.current_time() * framerate) {
			save_state.save(ps.get_particlearray());
			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";
		}
	}

    return 0;
}