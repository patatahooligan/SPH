#include "stdafx.h"

#include "XMLReader.h"

using namespace tinyxml2;

void get_constants_from_XML(XMLHandle& XML_root, CaseDef &case_def);
void get_geometry_from_XML(XMLHandle& XML_root, CaseDef &case_def);

CaseDef get_case_from_XML(const char * xml_filename) {
	CaseDef case_def;
	XMLDocument case_XML;

	// TODO: Consider how to implement error handling
	//       For now just throw the return value for debugging

	auto error = case_XML.LoadFile(xml_filename);
	if (error != XML_SUCCESS)
		throw error;

	XMLHandle root(case_XML.RootElement());
	assert(root.ToNode());

	get_constants_from_XML(root, case_def);

	get_geometry_from_XML(root, case_def);

	return case_def;
}

void get_constants_from_XML(XMLHandle& XML_root, CaseDef &case_def) {
	auto constants = XML_root.FirstChildElement("casedef").FirstChildElement("constantsdef");
	assert(constants.ToNode());

	if (auto gravity = constants.FirstChildElement("gravity").ToElement()) {
		case_def.gravity.x = gravity->FloatAttribute("x");
		case_def.gravity.y = gravity->FloatAttribute("y");
		case_def.gravity.z = gravity->FloatAttribute("z");
	}

	if (auto rhop0 = constants.FirstChildElement("rhop0").ToElement())
		case_def.rhop0 = rhop0->FloatAttribute("value");

	if (auto hswl = constants.FirstChildElement("hswl").ToElement()) {
		if (hswl->BoolAttribute("auto")) {
			// TODO
		}
		else
			case_def.hswl = hswl->FloatAttribute("value");
	}

	if (auto gamma = constants.FirstChildElement("gamma").ToElement())
		case_def.gamma = gamma->FloatAttribute("value");

	if (auto speedsystem = constants.FirstChildElement("speedsystem").ToElement()) {
		if (speedsystem->BoolAttribute("auto")) {
			// TODO
		}
		else
			case_def.speedsystem = speedsystem->FloatAttribute("value");
	}

	if (auto coefsound = constants.FirstChildElement("coefsound").ToElement())
		case_def.coefsound = coefsound->FloatAttribute("value)");

	if (auto speedsound = constants.FirstChildElement("speedsound").ToElement()) {
		if (speedsound->BoolAttribute("auto")) {
			case_def.speedsound = case_def.coefsound * case_def.speedsystem;
		}
		else
			case_def.speedsound = speedsound->FloatAttribute("value)");
	}

	if (auto coefh = constants.FirstChildElement("coefh").ToElement())
		case_def.coefh = coefh->FloatAttribute("value)");

	if (auto cflnumber = constants.FirstChildElement("cflnumber").ToElement())
		case_def.cflnumber = cflnumber->FloatAttribute("value)");
}

void get_geometry_from_XML(XMLHandle& XML_root, CaseDef &case_def) {
	auto geometry = XML_root.FirstChildElement("casedef").FirstChildElement("geometry");
	assert(geometry.ToNode());

	auto definition = geometry.FirstChildElement("definition").ToElement();
	assert(definition);

	case_def.particles.density = definition->FloatAttribute("dp");

	auto pointmin = definition->FirstChildElement("pointmin");
	assert(pointmin);
	case_def.particles.point_min.x = pointmin->FloatAttribute("x");
	case_def.particles.point_min.y = pointmin->FloatAttribute("y");
	case_def.particles.point_min.z = pointmin->FloatAttribute("z");

	auto pointmax = definition->FirstChildElement("pointmax");
	assert(pointmax);
	case_def.particles.point_max.x = pointmax->FloatAttribute("x");
	case_def.particles.point_max.y = pointmax->FloatAttribute("y");
	case_def.particles.point_max.z = pointmax->FloatAttribute("z");
}