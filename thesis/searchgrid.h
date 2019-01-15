#pragma once

#include "stdafx.h"

#include "particle.h"
#include "vec3f.h"
#include "massspringdamper.h"

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
		using cell_indices_container = boost::container::static_vector<index_pair, 54>;

	private:
		// We need a separate internal type because it uses vastly more memory.
		// In theory they can both be the same type but for optimization reasons they might not.
		using internal_cell_indices_container = std::vector<index_pair>;

		std::vector<ParticleProxy> proxies;
		Vec3f point_min, point_max, size;
		float cell_size;
		std::array<int, 3> grid_cells;
		internal_cell_indices_container cell_indices;

		std::array<int, 3> determine_number_of_cells () const
		{
			assert(size.x >= 0.0f && size.y >= 0.0f && size.z >= 0.0f);
			
			return {
				int(size.x / cell_size) + 1,
				int(size.y / cell_size) + 1,
				int(size.z / cell_size) + 1
			};
		}

		std::array<int, 3> determine_cell(const Vec3f &position) const {
			std::array<int, 3> cell;
			const Vec3f relative_position = position - point_min;

			for (size_t i = 0; i < 3; ++i) {
				// To robustly handle (slightly) out of bounds particles
				// group them to the nearest cell
				cell[i] = int(relative_position[i] / cell_size);
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

		void determine_order(const iter begin, const iter end) {
			const int container_size = std::distance(begin, end);
			proxies.resize(container_size);

			#pragma omp parallel for
			for (int i = 0; i < container_size; ++i)
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
		SearchGrid(const Vec3f &point_min, const Vec3f &point_max, const float cell_size) :
			point_min(point_min), point_max(point_max), size (point_max - point_min), cell_size(cell_size) {}

		void sort_containers(
			const iter target_begin, const iter target_end, iter parallel_begin,
			std::vector<MassSpringSystem>* fluid_spring_systems = nullptr, std::vector<MassSpringSystem>* boundary_spring_systems = nullptr)
		{
			size = point_max - point_min;
			grid_cells = determine_number_of_cells();

			determine_order(target_begin, target_end);

			if (fluid_spring_systems && boundary_spring_systems) {
				// Create an array of reverse proxies and use it to update the mass_spring_damper indices
				std::vector<int> reverse_proxies;
				reverse_proxies.resize(proxies.size());
				#pragma omp parallel for
				for (int i = 0; i < proxies.size(); ++i) {
					reverse_proxies[proxies[i].index] = i;
				}

				if (fluid_spring_systems) {
					for (auto &fluid_spring_system : *fluid_spring_systems) {
						#pragma omp parallel for
						for (int i = 0; i < fluid_spring_system.springs.size(); ++i) {
							auto &mass_spring_damper = fluid_spring_system.springs[i];
							auto
								&first = mass_spring_damper.particle_indices.first,
								&second = mass_spring_damper.particle_indices.second;
							first = reverse_proxies[first];
							second = reverse_proxies[second];
						}
					}
				}

				if (boundary_spring_systems) {
					for (auto &boundary_spring_system : *boundary_spring_systems) {
						#pragma omp parallel for
						for (int i = 0; i < boundary_spring_system.springs.size(); ++i) {
							auto &mass_spring_damper = boundary_spring_system.springs[i];
							auto
								&first = mass_spring_damper.particle_indices.first;
							first = reverse_proxies[first];
						}
					}
				}
			}

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

		void build_index(iter target_begin, iter target_end) {
			size = point_max - point_min;
			grid_cells = determine_number_of_cells();

			const int container_size = std::distance(target_begin, target_end);
			proxies.resize(container_size);

			#pragma omp parallel for
			for (int i = 0; i < container_size; ++i)
				proxies[i] = ParticleProxy(
					i,
					cell_coordinates_to_index(determine_cell(target_begin[i].position)));

			if (!std::is_sorted(proxies.begin(), proxies.end()))
				throw std::runtime_error("Can't build index on non-sorted container");

			determine_cell_indices();
		}

		void get_neighbor_indices(const Vec3f &position, cell_indices_container &container) const {
			const auto target_cell = determine_cell(position);

			for (int dim = 0; dim < 3; ++dim) {
				if (target_cell[dim] < -1 || target_cell[dim] > grid_cells[dim] + 1)
					return;
			}

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
						// If current cell is valid and non-empty, add it to neighbors
						const auto& current_cell_indices =
							cell_indices[cell_coordinates_to_index({ curr_x, curr_y, curr_z })];
						if (current_cell_indices.first != current_cell_indices.second)
							container.emplace_back(current_cell_indices);
					}
				}
			}
		}

		cell_indices_container get_neighbor_indices(const Vec3f &position) const {
			cell_indices_container neighbor_indices;

			get_neighbor_indices(position, neighbor_indices);
			return neighbor_indices;
		}

		void set_point_min(const Vec3f new_point) { point_min = new_point; }
		void set_point_max(const Vec3f new_point) { point_max = new_point; }
};
