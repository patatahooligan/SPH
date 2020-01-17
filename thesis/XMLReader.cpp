#include "../common/stdafx.h"

#include "../common/vec3f.h"
#include "../common/constants.h"
#include "XMLReader.h"
#include "../common/particle.h"
#include "massspringdamper.h"

using namespace tinyxml2;

Vec3f get_vec3f_from_element(const XMLElement& elem);
void set_element_to_vec3f(XMLElement& elem, const Vec3f &vec);

void get_constants_from_XML(XMLHandle& XML_root, CaseDef &case_def);
void get_geometry_from_XML(XMLHandle& XML_root, CaseDef &case_def);

Vec3f get_vec3f_from_element(const XMLElement& elem) {
	return { elem.FloatAttribute("x"), elem.FloatAttribute("y"), elem.FloatAttribute("z") };
}

void set_element_to_vec3f(XMLElement& elem, const Vec3f &vec) {
	elem.SetAttribute("x", vec.x);
	elem.SetAttribute("y", vec.y);
	elem.SetAttribute("z", vec.z);
}

CaseDef get_case_from_XML(const char*  xml_filename) {
	CaseDef case_def;
	XMLDocument case_XML;

	// TODO: Consider how to implement error handling
	//       For now just throw the return value for debugging

	auto error = case_XML.LoadFile(xml_filename);
	if (error != XML_SUCCESS)
		throw error;

	XMLHandle root(case_XML.RootElement());
	assert(root.ToNode());

	get_geometry_from_XML(root, case_def);

	get_constants_from_XML(root, case_def);

	case_def.particles.mass = case_def.rhop0 * std::pow(case_def.particles.density, 3);

	return case_def;
}

void get_constants_from_XML(XMLHandle& XML_root, CaseDef &case_def) {
	auto constants = XML_root.FirstChildElement("casedef").FirstChildElement("constantsdef");
	assert(constants.ToNode());

	if (auto alpha = constants.FirstChildElement("alpha").ToElement())
		case_def.alpha = alpha->FloatAttribute("value");

	if (auto gravity = constants.FirstChildElement("gravity").ToElement())
		case_def.gravity = get_vec3f_from_element(*gravity);

	if (auto rhop0 = constants.FirstChildElement("rhop0").ToElement())
		case_def.rhop0 = rhop0->FloatAttribute("value");

	if (auto spring = constants.FirstChildElement("spring").ToElement()) {
		case_def.spring.on = true;
		case_def.spring.damping = spring->FloatAttribute("damping");
		case_def.spring.max_length = spring->FloatAttribute("maxlength", case_def.particles.density);

		auto spring_system_element = spring->FirstChildElement("system");
		while (spring_system_element) {
			CaseDef::MassSpringSystem spring_system;

			spring_system.initial_k = spring_system_element->FloatAttribute("stiffness");
			spring_system.start_of_melting = spring_system_element->FloatAttribute("start_of_melting");
			spring_system.duration_of_melting = spring_system_element->FloatAttribute("duration_of_melting");
			spring_system.region.origin = get_vec3f_from_element(*(spring_system_element->FirstChildElement("pointmin")));
			const auto point_max = get_vec3f_from_element(*(spring_system_element->FirstChildElement("pointmax")));
			spring_system.region.size = point_max - spring_system.region.origin;

			case_def.spring.mass_spring_systems.push_back(spring_system);

			spring_system_element = spring_system_element->NextSiblingElement("system");
		}
	}

	if (auto ply = constants.FirstChildElement("ply").ToElement()) {
		do {
			const auto filename = ply->Attribute("filename");
			auto PLY_reader = vtkSmartPointer<vtkPLYReader>::New();
			PLY_reader->SetFileName(filename);
			PLY_reader->Update();
			CaseDef::PolyDataModel poly_data_model;
			poly_data_model.poly_data = PLY_reader->GetOutput();
			poly_data_model.offset = get_vec3f_from_element(*ply);
			poly_data_model.rotation = {
				ply->FloatAttribute("rotx"), ply->FloatAttribute("roty"), ply->FloatAttribute("rotz") };
			poly_data_model.scale = ply->FloatAttribute("scale", 1.0f);
			case_def.poly_data_models.emplace_back(poly_data_model);
			ply = ply->NextSiblingElement("ply");
		} while (ply);
	}

	if (auto hswl = constants.FirstChildElement("hswl").ToElement()) {
		if (hswl->BoolAttribute("auto")) {
			case_def.hswl = 0.0f;
			const Vec3f down_unit_vector = case_def.gravity.unit_vector();
			for (const auto &box : case_def.particle_boxes)
				if (box.type == CaseDef::CaseDefBox::Type::Fluid)
					case_def.hswl = std::max(
						std::abs(dot_product(box.size, down_unit_vector)),
						case_def.hswl);
			for (const auto &model : case_def.poly_data_models) {
				double bounds[6];
				model.poly_data->GetBounds(bounds);
				const Vec3f bounding_box = Vec3f{
					float(bounds[1] - bounds[0]),
					float(bounds[3] - bounds[2]),
					float(bounds[5] - bounds[4]) } * model.scale;
				case_def.hswl = std::max(
					std::abs(dot_product(bounding_box, down_unit_vector)),
					case_def.hswl);
			}
		}
		else
			case_def.hswl = hswl->FloatAttribute("value");
	}

	if (auto gamma = constants.FirstChildElement("gamma").ToElement())
		case_def.gamma = gamma->FloatAttribute("value");

	if (auto speedsystem = constants.FirstChildElement("speedsystem").ToElement()) {
		if (speedsystem->BoolAttribute("auto")) {
			case_def.speedsystem = std::sqrt(case_def.gravity.length() * case_def.hswl);
		}
		else
			case_def.speedsystem = speedsystem->FloatAttribute("value");
	}

	if (auto coefsound = constants.FirstChildElement("coefsound").ToElement())
		case_def.coefsound = coefsound->FloatAttribute("value");

	if (auto speedsound = constants.FirstChildElement("speedsound").ToElement()) {
		if (speedsound->BoolAttribute("auto")) {
			case_def.speedsound = case_def.coefsound * case_def.speedsystem;
		}
		else
			case_def.speedsound = speedsound->FloatAttribute("value)");
	}

	if (auto coefh = constants.FirstChildElement("coefh").ToElement())
		case_def.h = coefh->FloatAttribute("value") * float(std::sqrt(3)) * case_def.particles.density;

	const float
		deltap = 1.0f / 1.5f,
		w_deltap =
			(1.0f - 1.5f * deltap * deltap + 0.75f * deltap * deltap* deltap) /
			(pi * std::pow(case_def.h, 3));

	case_def.tensile_coef = 1.0f / w_deltap;

	if (auto cflnumber = constants.FirstChildElement("cflnumber").ToElement())
		case_def.cflnumber = cflnumber->FloatAttribute("value");

	if (auto friction = constants.FirstChildElement("friction").ToElement())
		case_def.friction_coef = friction->FloatAttribute("coefficient");
	else
		case_def.friction_coef = 0.0f;
}

void get_geometry_from_XML(XMLHandle& XML_root, CaseDef &case_def) {
	auto geometry = XML_root.FirstChildElement("casedef").FirstChildElement("geometry");
	assert(geometry.ToNode());

	auto definition = geometry.FirstChildElement("definition").ToElement();
	assert(definition);

	case_def.particles.density = definition->FloatAttribute("dp");

	// Minimum and maximum points of the domain, not the fluid
	auto pointmin = definition->FirstChildElement("pointmin");
	assert(pointmin);
	case_def.particles.point_min = get_vec3f_from_element(*pointmin);

	auto pointmax = definition->FirstChildElement("pointmax");
	assert(pointmax);
	case_def.particles.point_max = get_vec3f_from_element(*pointmax);

	// Execute the commands
	auto command_list = geometry.FirstChildElement("commands").FirstChildElement("mainlist");

	// The type of a box (fluid, boundary, void) is given in a separate command from the
	// geometry of the box, so we have to store its type
	using Type = CaseDef::CaseDefBox::Type;
	auto box_type = Type::Void;

	for (auto command = command_list.FirstChildElement().ToElement(); command;) {
		if (!std::strcmp(command->Name(), "setmkfluid"))
			box_type = Type::Fluid;
		else if (!std::strcmp(command->Name(), "setmkbound"))
			box_type = Type::Boundary;
		else if (!std::strcmp(command->Name(), "setmkvoid"))
			box_type = Type::Void;

		else if (!std::strcmp(command->Name(), "drawbox")) {			
			auto
				point = command->FirstChildElement("point"),
				size = command->FirstChildElement("size"),
				boxfill = command->FirstChildElement("boxfill");

			if (!point) throw std::runtime_error("No point given for drawbox");
			if (!size) throw std::runtime_error("No size given for drawbox");
			if (!boxfill) throw std::runtime_error("No boxfill mode given for drawbox");

			CaseDef::CaseDefBox box;
			box.origin = get_vec3f_from_element(*point);
			box.size = get_vec3f_from_element(*size);
			box.type = box_type;
			std::string_view boxfill_text(boxfill->GetText());
			if (boxfill_text.find("solid") != std::string_view::npos) box.fillmode.solid = true;
			if (boxfill_text.find("bottom") != std::string_view::npos) box.fillmode.bottom = true;
			if (boxfill_text.find("top") != std::string_view::npos) box.fillmode.top = true;
			if (boxfill_text.find("front") != std::string_view::npos) box.fillmode.front = true;
			if (boxfill_text.find("back") != std::string_view::npos) box.fillmode.back = true;
			if (boxfill_text.find("right") != std::string_view::npos) box.fillmode.right = true;
			if (boxfill_text.find("left") != std::string_view::npos) box.fillmode.left = true;

			case_def.particle_boxes.emplace_back(box);
		}

		command = command->NextSiblingElement();
	}
}


void save_particles_to_xml(XMLDocument &document, const ParticleConstIterator begin,
                           const ParticleConstIterator end, const char* element_name) {

	for (auto particle = begin; particle != end; ++particle) {
		auto particle_element = document.NewElement(element_name);

		auto position_element = document.NewElement("position");
		set_element_to_vec3f(*position_element, particle->position);
		particle_element->InsertEndChild(position_element);

		auto velocity_element = document.NewElement("velocity");
		set_element_to_vec3f(*velocity_element, particle->velocity);
		particle_element->InsertEndChild(velocity_element);

		auto density_element = document.NewElement("density");
		density_element->SetAttribute("value", particle->density);
		particle_element->InsertEndChild(density_element);

		document.InsertEndChild(particle_element);
	}
}

void save_spring_system_to_xml(XMLDocument &document, const MassSpringSystem&spring_system, const char* element_name) {
	auto system_element = document.NewElement(element_name);

	system_element->SetAttribute("initial_k", spring_system.initial_k);
	system_element->SetAttribute("start_of_melting", spring_system.start_of_melting);
	system_element->SetAttribute("duration_of_melting", spring_system.duration_of_melting);

	for (const auto &spring : spring_system.springs) {
		auto spring_element = document.NewElement("spring");
		spring_element->SetAttribute("restinglength", spring.resting_length);
		spring_element->SetAttribute("first", spring.particle_indices.first);
		spring_element->SetAttribute("second", spring.particle_indices.second);

		system_element->InsertEndChild(spring_element);
	}

	document.InsertEndChild(system_element);
}

void save_multiple_spring_systems_to_xml(XMLDocument &document, const std::vector<MassSpringSystem> &spring_systems,
                                         const char* element_prefix) {
	for (int i = 0; i < spring_systems.size(); ++i) {
		const std::string element_name = element_prefix + std::to_string(i);
		save_spring_system_to_xml(document, spring_systems[i], element_name.c_str());
	}
}

ParticleContainer load_particles_from_xml(const XMLDocument &document, const char* element_name) {
	ParticleContainer particles;
	auto element = document.FirstChildElement(element_name);
	while (element) {
		Particle particle;
		particle.position = get_vec3f_from_element(*(element->FirstChildElement("position")));
		particle.velocity = get_vec3f_from_element(*(element->FirstChildElement("velocity")));
		particle.density = element->FirstChildElement("density")->FloatAttribute("value");
		particles.push_back(particle);

		element = element->NextSiblingElement(element_name);
	}
	return particles;
}

MassSpringSystem load_spring_system_from_element(const XMLElement &spring_system_element) {
	MassSpringSystem spring_system;

	spring_system.initial_k = spring_system_element.FloatAttribute("initial_k");
	spring_system.start_of_melting = spring_system_element.FloatAttribute("start_of_melting");
	spring_system.duration_of_melting = spring_system_element.FloatAttribute("duration_of_melting");

	auto spring_element = spring_system_element.FirstChildElement("spring");
	while (spring_element) {
		MassSpringDamper spring;
		spring.resting_length = spring_element->FloatAttribute("restinglength");
		spring.particle_indices.first = spring_element->IntAttribute("first");
		spring.particle_indices.second = spring_element->IntAttribute("second");
		spring_system.springs.push_back(spring);

		spring_element = spring_element->NextSiblingElement("spring");
	}
	return spring_system;
}

MassSpringSystem load_spring_system_from_xml(const XMLDocument &document, const char* element_name) {
	return load_spring_system_from_element(*(document.FirstChildElement("element_name")));
}

std::vector<MassSpringSystem> load_multiple_spring_systems_from_xml(const XMLDocument &document, const char* element_prefix) {
	std::vector<MassSpringSystem> spring_systems;

	int i = 0;
	while (true) {
		const std::string element_name = element_prefix + std::to_string(i);
		const auto spring_system_element = document.FirstChildElement(element_name.c_str());
		if (!spring_system_element)
			break;

		spring_systems.push_back(load_spring_system_from_element(*spring_system_element));

		++i;
	}

	return spring_systems;
}