#include "stdafx.h"
#include "CppUnitTest.h"

#include <random>
#include <array>
#include <iostream>

#include "../common/vec3f.h"
#include "../common/particle.h"
#include "../thesis/searchgrid.h"

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
				num_of_particles = 10000;

			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dist(dimension_min, dimension_max);
			
			ParticleContainer particles, prev_particles;

			for (int i = 0; i < num_of_particles; ++i) {
				// Create two identical particle arrays (only position matters)
				Particle p;
				p.position = { float(dist(gen)), float(dist(gen)), float(dist(gen)) };
				particles.emplace_back(p);
				prev_particles.emplace_back(p);
			}

			// Keep a copy to test whether the output is a permutation of the input
			const auto unsorted = particles;

			SearchGrid grid{ point_min, point_max, h };

			grid.sort_containers(particles.begin(), particles.end(), prev_particles.begin());

			Assert::IsTrue(std::is_permutation(unsorted.begin(), unsorted.end(), particles.begin(),
				[](const Particle &lhs, const Particle &rhs) {
					return lhs.position == rhs.position;
				}));

			size_t neighbors = 0;
			for (int i = 0; i < num_of_particles; ++i) {
				// Assert that the arrays have been sorted in parallel
				Assert::IsTrue(particles[i].position == prev_particles[i].position);

				// Assert that every neighbor returned is in the proximity of the current particle
				const auto neighbor_indices = grid.get_neighbor_indices(particles[i].position);
				for (const auto& index_pair : neighbor_indices) {
					for (int j = index_pair.first; j < index_pair.second; ++j) {
						++neighbors;
						Vec3f relative_position =
							particles[i].position - particles[j].position;
						Assert::IsTrue(std::abs(relative_position.x) <= 2 * h);
						Assert::IsTrue(std::abs(relative_position.y) <= 2 * h);
						Assert::IsTrue(std::abs(relative_position.z) <= 2 * h);
					}
				}
			}
		}
	};
}