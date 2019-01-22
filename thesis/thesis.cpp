// thesis.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "physics.h"
#include "fileIO.h"
#include "XMLReader.h"

struct RunOptions {
	std::string case_filename;
	float time, output_period;

	std::optional<int> max_run_time;
	std::optional<std::string> input;
	std::optional<std::string> particles_output_filename, initial_state_filename;
};

auto get_options(int argc, char **argv) {
	cxxopts::Options cxx_options("SPH", "Smoothed Particle Hydrodynamics");
	RunOptions run_options;

	cxx_options.add_options()
		("particles-output", "Prefix for group of vtk files to hold the result of the simulation",
			cxxopts::value<std::string>())
		("t, time", "Time of simulation run in seconds",
			cxxopts::value<float>())
		("c, case", "XML file that defines the case",
			cxxopts::value<std::string>())
		("p, output-period", "interval in seconds between saved snapshots in binary output",
			cxxopts::value<float>())
		("i, initial-state", "A saved state from older simulation to be used as initial state of this one",
			cxxopts::value<std::string>())
		("max-run-time", "Maximum time in minutes the app can run before termination",
			cxxopts::value<float>())
		("snapshot", "Filename for snapshot of final step (default \"snapshot.bin\")",
			cxxopts::value<std::string>());
	
	auto const_argv = const_cast<const char**>(argv);       // Workaround for cxxopts bug!!!
	const auto result = cxx_options.parse(argc, const_argv);

	if (result.count("case") == 0)
		throw std::runtime_error(std::string("No option specified for case"));
	else
		run_options.case_filename = result["case"].as<std::string>();

	if (result.count("time") == 0)
		throw std::runtime_error(std::string("No options specified for time"));
	else
		run_options.time = result["time"].as<float>();

	if (result.count("output-period") == 0)
		run_options.output_period = 0;
	else
		run_options.output_period = result["output-period"].as<float>();

	if (result.count("initial-state") != 0)
		run_options.initial_state_filename = result["initial-state"].as<std::string>();

	if (result.count("particles-output") != 0)
		run_options.particles_output_filename = result["particles-output"].as<std::string>();

	if (!run_options.particles_output_filename)
		throw std::runtime_error(std::string("No option specified for output"));

	if (result.count("max-run-time") == 1)
		run_options.max_run_time = int(result["max-run-time"].as<float>() * 60.0f);

	if (result.count("input") == 1)
		run_options.input = result["input"].as<std::string>();

	return run_options;
}

int main(int argc, char **argv) {
	const auto options = get_options(argc, argv);

	const auto case_def = get_case_from_XML(options.case_filename.c_str());

	int output_step = 0;
	ParticleSystem ps = [&]() {
		if (options.initial_state_filename) {
			ParticleSystem::State initial_state;
			tinyxml2::XMLDocument initial_state_xml;
			initial_state_xml.LoadFile(options.initial_state_filename->c_str());
			initial_state.prev_fluid_particles = load_particles_from_xml(initial_state_xml, "prev_fluid_particle");
			initial_state.prev_boundary_particles = load_particles_from_xml(initial_state_xml, "prev_boundary_particle");
			initial_state.fluid_particles = load_particles_from_xml(initial_state_xml, "fluid_particle");
			initial_state.boundary_particles = load_particles_from_xml(initial_state_xml, "boundary_particle");
			initial_state.fluid_spring_systems = load_multiple_spring_systems_from_xml(initial_state_xml, "fluid_spring_system");
			initial_state.boundary_spring_systems = load_multiple_spring_systems_from_xml(initial_state_xml, "boundary_spring_system");
			
			if (const auto simulation_time_elem = initial_state_xml.FirstChildElement("time"))
				initial_state.simulation_time = simulation_time_elem->FloatAttribute("value");
			
			if (const auto verlet_step_elem = initial_state_xml.FirstChildElement("verlet-step"))
				initial_state.verlet_step = verlet_step_elem->FloatAttribute("value");

			if (const auto output_step_elem = initial_state_xml.FirstChildElement("output-step"))
				output_step = output_step_elem->IntAttribute("value");

			return ParticleSystem(case_def, std::move(initial_state));
		}
		else
			return ParticleSystem(case_def);
	}();

	// Save initial state only if this is a clean run
	if (output_step == 0)	{
		SaveVTK save_VTK(ps.get_boundary_begin(), ps.get_boundary_end());
		if (options.particles_output_filename)
			save_VTK.save_particles(*options.particles_output_filename + "-boundary");
	}
	
	bool user_exit = false;
	auto& now = std::chrono::steady_clock::now;
	const auto start_time = now();
	while (
		ps.current_time() < options.time &&
		!user_exit &&
		(!options.max_run_time || now() - start_time < std::chrono::seconds{ *(options.max_run_time) }))
	{
		while (output_step * options.output_period <= ps.current_time()) {
			SaveVTK save_VTK(ps.get_fluid_begin(), ps.get_fluid_end());
			if (options.particles_output_filename)
				save_VTK.save_particles(*options.particles_output_filename + "-fluid" + std::to_string(output_step));

			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";

			++output_step;
		}

		ps.simulation_step();

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
		<<"    (of requested " << options.time << " s)\n\n";

	std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";

	std::cout << "Verlet steps integrated : " << ps.current_step() << "\n";
	std::cout << "Steps saved in output file : " << output_step << "\n\n";

	std::cout << "Saving two snapshots of particles\n";

	tinyxml2::XMLDocument final_state_output;
	const auto final_state = ps.get_current_state();

	save_particles_to_xml(final_state_output, final_state.prev_fluid_particles.begin(), final_state.prev_fluid_particles.end(),
	                      "prev_fluid_particle");
	save_particles_to_xml(final_state_output, final_state.fluid_particles.begin(), final_state.fluid_particles.end(),
	                      "fluid_particle");
	save_particles_to_xml(final_state_output, final_state.prev_boundary_particles.begin(), final_state.prev_boundary_particles.end(),
	                      "prev_boundary_particle");
	save_particles_to_xml(final_state_output, final_state.boundary_particles.begin(), final_state.boundary_particles.end(),
	                      "boundary_particle");
	save_multiple_spring_systems_to_xml(final_state_output, final_state.fluid_spring_systems, "fluid_spring_system");
	save_multiple_spring_systems_to_xml(final_state_output, final_state.boundary_spring_systems, "boundary_spring_system");

	append_element_to_xml(final_state_output, "time", ps.current_time());
	append_element_to_xml(final_state_output, "verlet-step", ps.current_step());
	append_element_to_xml(final_state_output, "output-step", output_step);

	final_state_output.SaveFile("final_state.xml");

	return 0;
}