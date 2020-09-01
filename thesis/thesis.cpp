// thesis.cpp : Defines the entry point for the console application.

#include "../common/stdafx.h"

#include "boost/program_options.hpp"

#include "physics.h"
#include "../common/fileIO.h"
#include "XMLReader.h"


std::atomic<bool> user_exit = false;

extern "C" void interrupt_handler(int ) {
	// Notify the main loop to exit
	user_exit = true;
}

struct ProgramOptions {
	std::string case_filename;
	float time, output_period;

	std::optional<int> max_run_time;
	std::optional<std::string> input;
	std::optional<std::string> particles_output_filename, initial_state_filename;

	bool exit = false;
};

auto get_options(int argc, char **argv) {
	namespace options = boost::program_options;

	options::options_description description(
		"Smooth particle hydrodynamics solver. Available options"
	);

	description.add_options()
		("help,h",
			"Produce help message")
		("particles-output,o",
			options::value<std::string>()->default_value("output"),
			"Prefix for group of vtk files to hold the result of the simulation")
		("time,t",
			options::value<float>()->required(),
			"Time of simulation run in seconds")
		("case,c",
			options::value<std::string>()->required(),
			"XML file that defines the case")
		("output-period,p",
			options::value<float>()->required(),
			"interval in seconds between saved snapshots in binary output")
		("initial-state,i",
			options::value<std::string>(),
			"A saved state from older simulation to be used as initial state "
			"of this one")
		("max-run-time",
			options::value<float>(),
			"Maximum time in minutes the app can run before termination")
		;

	options::variables_map variables_map;
	options::store(options::parse_command_line(argc, argv, description), variables_map);
	options::store(
		options::command_line_parser(argc, argv).
			options(description).
			run(),
		variables_map);

	ProgramOptions program_options;

	if (variables_map.count("help")) {
		std::cout << description << "\n";
		program_options.exit = true;
		return program_options;
	}

	options::notify(variables_map);
	if (variables_map.count("particles-output")) {
		program_options.particles_output_filename =
			variables_map["particles-output"].as<std::string>();
	}
	if (variables_map.count("time")) {
		program_options.time =
			variables_map["time"].as<float>();
	}
	if (variables_map.count("case")) {
		program_options.case_filename =
			variables_map["case"].as<std::string>();
	}
	if (variables_map.count("output-period")) {
		program_options.output_period =
			variables_map["output-period"].as<float>();
	}
	if (variables_map.count("initial-state")) {
		program_options.initial_state_filename =
			variables_map["initial-state"].as<std::string>();
	}
	if (variables_map.count("max-run-time")) {
		program_options.max_run_time =
			variables_map["max-run-time"].as<float>();
	}

	return program_options;
}

int main(int argc, char **argv) {
	const auto options = get_options(argc, argv);

	if (options.exit)
		return 0;

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
		if (options.particles_output_filename)
			save_particles(ps.get_boundary_begin(), ps.get_boundary_end(), *options.particles_output_filename + "-boundary");
	}

	auto& now = std::chrono::steady_clock::now;
	const auto start_time = now();
	while (
		ps.current_time() < options.time &&
		!user_exit &&
		(!options.max_run_time || now() - start_time < std::chrono::seconds{ *(options.max_run_time) }))
	{
		while (output_step * options.output_period <= ps.current_time()) {
			if (options.particles_output_filename)
				save_particles(ps.get_fluid_begin(), ps.get_fluid_end(), *options.particles_output_filename + "-fluid" + std::to_string(output_step));

			std::cout << "Snapshot saved at time " << ps.current_time() << "s\n";

			++output_step;
		}

		ps.simulation_step();
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
