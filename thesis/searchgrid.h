#pragma once

#include "stdafx.h"

#include "particle.h"
#include "vec3f.h"

struct ParticleProxy {
	int index, cell;
	ParticleProxy(size_t index, size_t cell) :
		index(index),
		cell(cell) {}

	bool operator<(const ParticleProxy& other) const {
		return cell < other.cell;
	}
};

class SearchGrid {
	private:
		std::vector<ParticleProxy> proxies;
		const Vec3f point_min, point_max, size;
		const float h;
		const std::array<int, 3> grid_cells;

		std::array<int, 3> determine_number_of_cells(
			const Vec3f &size, const float h)
		{
			const Vec3f size = point_max - point_min;
			assert(size.x >= 0.0f, size.y >= 0.0f, z >= 0.0f);
			
			return {
				int(size.x / h) + 1,
				int(size.y / h) + 1,
				int(size.z / h) + 1
			};
		}

		int determine_cell(const Vec3f &position) {
			std::array<int, 3> cell;
			const Vec3f relative_position = position - size;

			for (size_t i = 0; i < 3; ++i) {
				// To robustly handle (slightly) out of bounds particles
				// group them to the nearest cell
				cell[i] = int(position[i] / h);
				cell[i] = std::max(cell[i], 0);
				cell[i] = std::min(cell[i], grid_cells[i] - 1);
			}

			return
				cell[0] +
				cell[1] * grid_cells[0] +
				cell[2] * grid_cells[0] * grid_cells[1];
		}

	public:
		using iter = ParticleContainer::iterator;

		SearchGrid(const Vec3f &point_min, const Vec3f &point_max, const float h) :
			point_min(point_min), point_max(point_max), size (point_max - point_min), h(h),
			grid_cells(determine_number_of_cells(size, h)) {}

		void determine_order(iter begin, iter end) {
			int size = std::distance(begin, end);
			proxies.resize(std::distance(begin, end));

			#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				proxies[i] = ParticleProxy(i, determine_cell(begin[i].position));

			// Get the proxies in proper order
			std::sort(begin, end);
		}

		void sort_containers(std::array<iter, 3> begin_it) {
			// Sort the containers in place
			// Precondition: determine_order() has been called to sort proxies
			// proxies[i] holds the index of the item that needs to be in i-th place

			// IMPORTANT: this destroys the proxies vector including the cell information
			//            indices that point to the cells must have been created before this call

			for (int i = 0; size_t(i) < proxies.size(); ++i) {
				// Keep a copy of [i]
				std::array<Particle, begin_it.size()> current_vec;
				for (size_t j = 0; j < begin_it.size(); ++j)
					current_vec[j] = (begin_it[j])[i];

				// Start by filling [i]
				auto current = i;

				// While the element we have to place in [current] is not the
				// one we are holding (the one initially in [i])
				while (i != proxies[current].index) {
					// [next] is the element that should be placed in [current]
					// It is also the next position we will have to fill
					int next = proxies[current].index;

					// Move the elements that belong to it[current]
					for (auto& it : begin_it)
						it[current] = it[next];

					// Now that [current] is occupied by the correct value,
					// make proxies[current] point to it
					proxies[current].index = proxies[next].index;

					// We've left a hole in [next], so that's the next index
					// we should handle
					current = next;
				}

				// When i == current, we've iterated to the index that current_vec
				// should be placed in
				for (size_t j = 0; j < begin_it.size(); ++j)
					(begin_it[j])[current] = current_vec[j];

				proxies[current].index = current;
			}
		}
};
