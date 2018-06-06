#pragma once

#include "stdafx.h"

#include <algorithm>
#include <execution>
#include <cassert>

#include "particle.h"
#include "vec3f.h"

struct ParticleProxy {
	int index, cell;

	ParticleProxy() = default;
	ParticleProxy(int index, int cell) :
		index(index),
		cell(cell) {}

	bool operator<(const ParticleProxy& other) const {
		return cell < other.cell;
	}
};

class SearchGrid {
	public:
		using iter = ParticleContainer::iterator;
		using index_pair = std::pair<int, int>;
		using cell_indices_container = std::vector<index_pair>;

	private:
		std::vector<ParticleProxy> proxies;
		Vec3f point_min, point_max, size;
		float h;
		std::array<int, 3> grid_cells;
		cell_indices_container cell_indices;

		std::array<int, 3> determine_number_of_cells (
			const Vec3f &size, const float h) const
		{
			assert(size.x >= 0.0f && size.y >= 0.0f && size.z >= 0.0f);
			
			return {
				int(size.x / h) + 1,
				int(size.y / h) + 1,
				int(size.z / h) + 1
			};
		}

		std::array<int, 3> determine_cell(const Vec3f &position) const {
			std::array<int, 3> cell;
			const Vec3f relative_position = position - point_min;

			for (size_t i = 0; i < 3; ++i) {
				// To robustly handle (slightly) out of bounds particles
				// group them to the nearest cell
				cell[i] = int(position[i] / h);
				cell[i] = std::max(cell[i], 0);
				cell[i] = std::min(cell[i], grid_cells[i] - 1);
			}

			return cell;
		}

		int cell_coordinates_to_index(const std::array<int, 3> &cell) const {
			// Map cell's 3D position to a 1D index
			return
				cell[2] +
				cell[1] * grid_cells[2] +
				cell[0] * grid_cells[1] * grid_cells[2];
		}

		void determine_order(iter begin, iter end) {
			const int size = std::distance(begin, end);
			proxies.resize(size);

			#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				proxies[i] = ParticleProxy(
					i,
					cell_coordinates_to_index(determine_cell(begin[i].position)));

			// Sort only the proxies, not the actual Particles
			std::sort(proxies.begin(), proxies.end());
		}

		void determine_cell_indices() {
			const auto num_of_cells = grid_cells[0] * grid_cells[1] * grid_cells[2];
			cell_indices.resize(num_of_cells);

			int proxy = 0;
			for (int cell = 0; cell < num_of_cells; ++cell) {
				cell_indices[cell].first = proxy;
				while (proxy < proxies.size() && proxies[proxy].cell == cell) {
					++proxy;
				}
				cell_indices[cell].second = proxy;
			}
		}

	public:
		SearchGrid(const Vec3f &point_min, const Vec3f &point_max, const float h) :
			point_min(point_min), point_max(point_max), size (point_max - point_min), h(h),
			grid_cells(determine_number_of_cells(size, h)) {}

		void sort_containers(iter target_begin, iter target_end, iter parallel_begin) {
			// Sort the containers in place based on [beg_it[0], end)

			determine_order(target_begin, target_end);

			determine_cell_indices();

			for (int i = 0; size_t(i) < proxies.size(); ++i) {
				// Keep a copy of [i]
				std::array<Particle, 2> current_vec;
				const Particle
					target_temp = target_begin[i],
					parallel_temp = parallel_begin[i];

				// Start by filling [i]
				auto current = i;

				// While the element we have to place in [current] is not the
				// one we are holding (the one initially in [i])
				while (i != proxies[current].index) {
					// [next] is the element that should be placed in [current]
					// It is also the next position we will have to fill
					int next = proxies[current].index;

					// Move the elements that belong to it[current]
					target_begin[current] = target_begin[next];
					parallel_begin[current] = parallel_begin[next];

					// Now that [current] is occupied by the correct value,
					// make proxies[current] point to it
					proxies[current].index = current;

					// We've left a hole in [next], so that's the next index
					// we should handle
					current = next;
				}

				// When i == current, we've iterated to the index that the temp particles
				// should be placed in
				target_begin[current] = target_temp;
				parallel_begin[current] = parallel_temp;

				proxies[current].index = current;
			}
		}

		void get_neighbor_indices(const Vec3f &position, cell_indices_container &container) const {
			const auto target_cell = determine_cell(position);
			const auto
				min_x = std::max(target_cell[0] - 1, 0),
				min_y = std::max(target_cell[1] - 1, 0),
				min_z = std::max(target_cell[2] - 1, 0),

				max_x = std::min(target_cell[0] + 1, grid_cells[0] - 1),
				max_y = std::min(target_cell[1] + 1, grid_cells[1] - 1),
				max_z = std::min(target_cell[2] + 1, grid_cells[2] - 1);

			for (int curr_x = min_x; curr_x <= max_x; ++curr_x) {
				for (int curr_y = min_y; curr_y <= max_y; ++curr_y) {
					for (int curr_z = min_z; curr_z <= max_z; ++curr_z) {
						// If current cell is valid (potentially empty), add it to neighbors
						container.emplace_back(
							cell_indices[cell_coordinates_to_index({ curr_x, curr_y, curr_z })]);
					}
				}
			}
		}

		cell_indices_container get_neighbor_indices(const Vec3f &position) const {
			cell_indices_container neighbor_indices;

			get_neighbor_indices(position, neighbor_indices);
			return neighbor_indices;
		}
};
