#include "stdafx.h"

#include <random>

#include "CppUnitTest.h"
#include "particle.h"
#include "../thesis/savestate.h"
#include "loadstate.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace savestateloadstatetest
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			constexpr int
				num_of_particles = 500,
				num_of_fluid_particles = num_of_particles/2,
				num_of_snapshots = 20;
			const float rand_float_min = -10.0f, rand_float_max = 10.0f;

			std::vector<ParticleContainer> to_save;

			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dist(rand_float_min, rand_float_max);
			for (int i = 0; i < num_of_snapshots; ++i) {
				ParticleContainer snapshot;
				for (int j = 0; j < num_of_particles; ++j) {
					Particle p;
					p.position = { float(dist(gen)), float(dist(gen)), float(dist(gen)) };
					p.velocity = { float(dist(gen)), float(dist(gen)), float(dist(gen)) };
					p.density = float(dist(gen));

					snapshot.emplace_back(p);
				}
				to_save.emplace_back(std::move(snapshot));
			}

			// Save vector to file
			SaveState::Mode modes[] = { SaveState::Mode::Full, SaveState::Mode::Position };
			for (const auto mode : modes) {
				const auto filename = "savestate-loadstate-test.bin";
				{
					SaveState save_state{ filename, num_of_fluid_particles, num_of_particles };
					for (const auto& snapshot : to_save)
						save_state.save(snapshot, mode);
				}
				LoadState load_state{ filename };
				Assert::AreEqual(load_state.get_num_of_fluid_particles(), num_of_fluid_particles);
				Assert::AreEqual(load_state.get_num_of_particles(), num_of_particles);
				for (int i = 0; i < num_of_snapshots; ++i) {
					auto snapshot = load_state.load(mode);
					Assert::AreEqual(snapshot.size(), size_t(num_of_particles));

					for (int j = 0; j < num_of_particles; ++j) {
						Assert::IsTrue(snapshot[j].position == to_save[i][j].position);
						if (mode == SaveState::Mode::Full) {
							Assert::IsTrue(snapshot[j].velocity == to_save[i][j].velocity);
							Assert::IsTrue(snapshot[j].density == to_save[i][j].density);
						}
					}
				}
			}
		}

	};
}