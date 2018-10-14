#pragma once

#include "constants.h"
#include "particle.h"
#include "massspringdamper.h"

CaseDef get_case_from_XML(const char* xml_filename);

void save_particles_to_xml(ParticleContainer::const_iterator begin, ParticleContainer::const_iterator end,
                           const char* xml_filename, const char* element_name);

void save_springs_to_xml(MassSpringConstIterator begin, MassSpringConstIterator end,
	const char* xml_filename, const char* element_name);
