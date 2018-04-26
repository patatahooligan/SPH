#include "stdafx.h"
#include "CppUnitTest.h"

#include <random>
#include <array>
#include <iostream>

#include "..\thesis\vec3f.h"
#include "..\thesis\particle.h"
#include "..\thesis\searchgrid.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace gridsearchtest {		
	TEST_CLASS(UnitTest1) {
	public:
		
		TEST_METHOD(SquareBoxNeighborSearch) {
			constexpr float
				h = 1.0f,
				dimension_min = 0.0f,
				dimension_max = 10.0f;
			constexpr Vec3f
				point_min{ dimension_min, dimension_min, dimension_min },
				point_max{ dimension_max, dimension_max, dimension_max };
			constexpr size_t
				num_of_particles = 500;

			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dist(dimension_min, dimension_max);
			
			std::array<ParticleContainer, 3> particles;
			for (auto& container : particles)
				container.reserve(num_of_particles);

			for (int i = 0; i < num_of_particles; ++i) {
				// Create three identical particle arrays (only position matters)
				Particle p;
				p.position = { float(dist(gen)), float(dist(gen)), float(dist(gen)) };
				for (auto& particle_container : particles)
					particle_container.emplace_back(p);
			}

			SearchGrid grid{ point_min, point_max, h };

			std::array<SearchGrid::iter, 3> begin_it{
				particles[0].begin(),
				particles[1].begin(),
				particles[2].begin() };
			grid.sort_containers(begin_it, particles[0].end());

			size_t neighbors = 0;
			for (int i = 0; i < num_of_particles; ++i) {
				// Assert that the arrays have been sorted in parallel
				Assert::IsTrue(particles[0][i].position == particles[1][i].position);
				Assert::IsTrue(particles[0][i].position == particles[2][i].position);

				// Assert that every neighbor returned is in the proximity of the current particle
				auto neighbor_indices = grid.get_neighbor_indices(particles[0][i].position);
				for (const auto& index_pair : neighbor_indices) {
					for (int j = index_pair.first; j < index_pair.second; ++j) {
						++neighbors;
						Vec3f relative_position =
							particles[0][i].position - particles[0][j].position;
						Assert::IsTrue(std::abs(relative_position.x) <= 2 * h);
						Assert::IsTrue(std::abs(relative_position.y) <= 2 * h);
						Assert::IsTrue(std::abs(relative_position.z) <= 2 * h);
					}
				}
			}

			std::cout << "Number of neighbors found: " << neighbors << '\n';
		}
	};
}