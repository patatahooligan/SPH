#include "stdafx.h"

#include "vec3f.h"
#include "constants.h"
#include "XMLReader.h"

using namespace tinyxml2;

Vec3f get_vec3f_from_element(XMLElement& elem);

void get_constants_from_XML(XMLHandle& XML_root, CaseDef &case_def);
void get_geometry_from_XML(XMLHandle& XML_root, CaseDef &case_def);

Vec3f get_vec3f_from_element(XMLElement& elem) {
	return { elem.FloatAttribute("x"), elem.FloatAttribute("y"), elem.FloatAttribute("z") };
}

CaseDef get_case_from_XML(const std::string_view xml_filename) {
	CaseDef case_def;
	XMLDocument case_XML;

	// TODO: Consider how to implement error handling
	//       For now just throw the return value for debugging

	auto error = case_XML.LoadFile(xml_filename.data());
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

	if (auto gravity = constants.FirstChildElement("gravity").ToElement())
		case_def.gravity = get_vec3f_from_element(*gravity);

	if (auto rhop0 = constants.FirstChildElement("rhop0").ToElement())
		case_def.rhop0 = rhop0->FloatAttribute("value");

	if (auto hswl = constants.FirstChildElement("hswl").ToElement()) {
		if (hswl->BoolAttribute("auto")) {
			case_def.hswl = 0.0f;
			const Vec3f down_unit_vector = case_def.gravity.unit_vector();
			for (auto &box : case_def.particle_boxes)
				if (box.type == CaseDef::Box::Type::Fluid)
					case_def.hswl = std::max(
						std::abs(dot_product(box.size, down_unit_vector)),
						case_def.hswl);
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

	if (auto cflnumber = constants.FirstChildElement("cflnumber").ToElement())
		case_def.cflnumber = cflnumber->FloatAttribute("value");
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
	using Type = CaseDef::Box::Type;
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

			CaseDef::Box box;
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