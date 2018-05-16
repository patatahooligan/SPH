#pragma once

#include "stdafx.h"

#include <algorithm>
#include <execution>

#include "boost/container/static_vector.hpp"

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
		using static_cell_indices_container = boost::container::static_vector<index_pair, 54>;

	private:
		std::vector<ParticleProxy> proxies;
		const Vec3f point_min, point_max, size;
		const float h;
		const std::array<int, 3> grid_cells;
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
			const Vec3f relative_position = position - size;

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
				cell[0] +
				cell[1] * grid_cells[0] +
				cell[2] * grid_cells[0] * grid_cells[1];
		}

		void determine_order(iter begin, iter end) {
			int size = std::distance(begin, end);
			proxies.resize(std::distance(begin, end));

			#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				proxies[i] = ParticleProxy(
					i,
					cell_coordinates_to_index(determine_cell(begin[i].position)));

			// Sort only the proxies, not the actual Particles
			std::sort(std::execution::par, proxies.begin(), proxies.end());
		}

		void determine_cell_indices() {
			auto num_of_cells = grid_cells[0] * grid_cells[1] * grid_cells[2];
			cell_indices.resize(num_of_cells);

			size_t proxy = 0;
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

		SearchGrid(SearchGrid&&) = default;

		void sort_containers(std::array<iter, 3> begin_it, iter end) {
			// Sort the containers in place based on [beg_it[0], end)

			determine_order(begin_it[0], end);

			determine_cell_indices();

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
					proxies[current].index = current;

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

		void get_neighbor_indices(const Vec3f &position, static_cell_indices_container &container) const {
			const auto target_cell = determine_cell(position);
			const auto
				&x = target_cell[0],
				&y = target_cell[1],
				&z = target_cell[2];

			for (int curr_z = std::max(z - 1, 0); curr_z <= std::min(z + 1, grid_cells[2]); ++curr_z) {
				for (int curr_y = std::max(y - 1, 0); curr_y <= std::min(y + 1, grid_cells[1]); ++curr_y) {
					for (int curr_x = std::max(x - 1, 0); curr_x <= std::min(x + 1, grid_cells[0]); ++curr_x) {
						// If current cell is valid (potentially empty), add it to neighbors
						container.push_back(
							cell_indices[cell_coordinates_to_index({ curr_x, curr_y, curr_z })]);
					}
				}
			}
		}

		static_cell_indices_container get_neighbor_indices(const Vec3f &position) const {
			static_cell_indices_container neighbor_indices;

			get_neighbor_indices(position, neighbor_indices);
			return neighbor_indices;
		}
};
