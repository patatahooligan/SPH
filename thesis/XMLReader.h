#pragma once

#include "stdafx.h"

#include "constants.h"
#include "particle.h"
#include "massspringdamper.h"

CaseDef get_case_from_XML(const char* xml_filename);

void save_particles_to_xml(tinyxml2::XMLDocument &document, ParticleConstIterator begin, ParticleConstIterator end,
                           const char* element_name);

void save_spring_system_to_xml(tinyxml2::XMLDocument &document, const MassSpringSystem&spring_system, const char* element_name);

void save_multiple_spring_systems_to_xml(tinyxml2::XMLDocument &document, const std::vector<MassSpringSystem> &spring_systems,
                                         const char* element_prefix);

template<class DataT>
void append_element_to_xml(tinyxml2::XMLDocument &document, const char* element_name, const DataT &element_value) {
	auto element = document.NewElement(element_name);
	element->SetAttribute("value", element_value);
	document.InsertEndChild(element);
}

ParticleContainer load_particles_from_xml(const tinyxml2::XMLDocument &document, const char* element_name);

MassSpringSystem load_spring_system_from_xml(const tinyxml2::XMLDocument &document, const char* element_name);

std::vector<MassSpringSystem> load_multiple_spring_systems_from_xml(const tinyxml2::XMLDocument &document, const char* element_prefix);