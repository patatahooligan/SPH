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
	bool user_exit = false;
	while (ps.current_time() < target_time && !user_exit) {
		ps.simulation_step();
		if (save_state.get_step() < ps.current_time() * framerate) {
			save_state.save(ps.get_particlearray(), SaveState::Mode::Position);
			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";
		}

		while (_kbhit()) {
			auto key = _getch();
			if (key == 27)
				user_exit = true;
		}
	}
	std::cout << '\n';

	if (user_exit)
		std::cout << "Cancelled by user\n\n";

	std::cout << "Time of simulation :" << ps.current_time() << " s\n"
		<<"    (of requested " << target_time << " s)\n\n";

	std::cout << "Steps saved in binary file : " << save_state.get_step() << '\n\n';

	std::cout << "Saving snapshot of particles at moment of termination\n"
			<< "Saved in order (previous step, current step)\n\n";

	SaveState("snapshot.bin", ps.get_num_of_particles())
		.save(ps.get_previous_particlearray(), SaveState::Mode::Full)
		.save(ps.get_particlearray(), SaveState::Mode::Full);

	return 0;
}