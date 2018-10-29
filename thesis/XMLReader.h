#pragma once

#include "stdafx.h"

#include "constants.h"
#include "particle.h"
#include "massspringdamper.h"

CaseDef get_case_from_XML(const char* xml_filename);

void save_particles_to_xml(tinyxml2::XMLDocument &document, ParticleConstIterator begin, ParticleConstIterator end,
                           const char* element_name);

void save_springs_to_xml(tinyxml2::XMLDocument &document, MassSpringConstIterator begin, MassSpringConstIterator end,
                         const char* element_name);

template<class DataT>
void append_element_to_xml(tinyxml2::XMLDocument &document, const char* element_name, const DataT &element_value) {
	auto element = document.NewElement(element_name);
	element->SetAttribute("value", element_value);
	document.InsertEndChild(element);
}