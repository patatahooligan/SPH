#include "stdafx.h"

#include "physics.h"
#include "constants.h"
#include "vec3f.h"

CaseDef::Box get_particle_axis_aligned_bounding_box(ParticleConstIterator begin, ParticleConstIterator end) {
	constexpr auto
		float_max = std::numeric_limits<float>::max(),
		float_lowest = std::numeric_limits<float>::lowest();

	Vec3f
		point_min = { float_max, float_max, float_max },
		point_max = { float_lowest, float_lowest, float_lowest };

	point_min = std::transform_reduce(std::execution::seq, begin, end, point_min,
		[](const Vec3f position1, const Vec3f position2) -> Vec3f {
		return {
			std::min(position1.x, position2.x),
			std::min(position1.y, position2.y),
			std::min(position1.z, position2.z)
		};
	},
		[](const Particle &particle) ->Vec3f {
		return particle.position;
	});

	point_max = std::transform_reduce(std::execution::seq, begin, end, point_max,
		[](const Vec3f position1, const Vec3f position2) -> Vec3f {
		return {
			std::max(position1.x, position2.x),
			std::max(position1.y, position2.y),
			std::max(position1.z, position2.z)
		};
	},
		[](const Particle &particle) ->Vec3f {
		return particle.position;
	});

	return { point_min, point_max - point_min };
}

CubicSpline::CubicSpline(const float h):
	h(h), a(1.0f / (pi * std::pow(h, 3))) {}

float CubicSpline::operator()(const float length) const {
	const float q = length / h;

	if (q < 1.0f)
		return a * (1 - (3.0f / 2.0f) * (q*q) + (3.0f / 4.0f) * std::pow(q, 3));
	else
		return a * (std::pow(2 - q, 3) / 4.0f);
}

float CubicSpline::operator()(const Vec3f &r) const {
	return (*this)(r.length());
}

float CubicSpline::gradient_coef(float length) const {
	const float
		q = length / h;

	const float
		dq = 1 / (h * h * q),
		dw = a * [&]() {
		if (q < 1.0f)
			return (9.0f / 4.0f) * q * q - 3.0f * q;
		else
			return -(3.0f / 4.0f) * (2.0f - q) * (2.0f - q);
	}();

	return dw * dq;
}

Vec3f CubicSpline::gradient(const Vec3f &r) const {
	const float
		q = r.length() / h;

	// It is easier and faster to handle these edge cases separately
	if (q == 0.0f || q > 2.0f)
		return { 0.0f, 0.0f, 0.0f };

	const float
		dq = 1 / (h * h * q),
		dw = a * [&]() {
		if (q < 1.0f)
			return (9.0f / 4.0f) * q * q - 3.0f * q;
		else if (q < 2.0f)
			return -(3.0f / 4.0f) * (2.0f - q) * (2.0f - q);
		else
			return 0.0f;
	}();

	return dw * dq * r;
}


CubicSplinePrecalculated::CubicSplinePrecalculated(const float h, const int resolution = 32767):
	h(h), resolution(resolution), step((2.0f * h) / resolution)
{
	values.resize(resolution + 1);
	gradient_values.resize(resolution + 1);

	CubicSpline spline(h);

	for (int i = 0; i <= resolution; ++i) {
		values[i] = spline(i * step);
		gradient_values[i] = spline.gradient_coef(i * step);
	}
}

float CubicSplinePrecalculated::operator()(const float length) const {
	return values[length / step];
}

float CubicSplinePrecalculated::operator()(const Vec3f & r) const {
	return (*this)(r.length());
}

Vec3f CubicSplinePrecalculated::gradient(const Vec3f & r, const float length) const {
	return gradient_values[length / step] * r;
}

Vec3f CubicSplinePrecalculated::gradient(const Vec3f & r) const {
	return gradient(r, r.length());
}


class ParticleGenerator {
	private:
		CaseDef case_def;
		ParticleContainer fluid_particles, boundary_particles;

		void fillPLY(CaseDef::PolyDataModel model, ParticleContainer &target_container) {
			double bounds[6];
			model.poly_data->ComputeBounds();
			model.poly_data->GetBounds(bounds);
			const double
				&x_min = bounds[0], &x_max = bounds[1],
				&y_min = bounds[2], &y_max = bounds[3],
				&z_min = bounds[4], &z_max = bounds[5];

			const float&
				density = case_def.particles.density;

			auto select_enclosed_points = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
			select_enclosed_points->SetTolerance(std::numeric_limits<double>::min());
			select_enclosed_points->SetCheckSurface(true);
			select_enclosed_points->Initialize(model.poly_data);
			for (float x = x_min; x <= x_max; x += density / model.scale) {
				for (float y = y_min; y <= y_max; y += density / model.scale) {
					for (float z = z_min; z <= z_max; z += density / model.scale) {
						if (select_enclosed_points->IsInsideSurface(x, y, z)) {
							Particle p;
							const auto& rot = model.rotation;
							p.position =
								Vec3f{ x, y, z }.rotate_x(rot.x).rotate_y(rot.y).rotate_z(rot.z)
								* model.scale + model.offset;
							p.density = case_def.rhop0;
							target_container.emplace_back(p);
						}
					}
				}
			}
			select_enclosed_points->Complete();
		}

		void remove_particles(const CaseDef::Box &box, ParticleContainer &target_container) {
			target_container.erase(
				std::remove_if(std::execution::par_unseq,
					target_container.begin(), target_container.end(), [&box](const Particle& p) {
						return box.contains(p.position);
					}),
				target_container.end());
		}

		void fillbox(const CaseDef::CaseDefBox &box, ParticleContainer& target_container) {
			// Because we might need to fill any combination of faces of a box, break it down
			const float& density = case_def.particles.density;

			const std::vector<CaseDef::Box> targets = [&]() {
				std::vector<CaseDef::Box> targets;

				if (box.fillmode.solid) {
					// If box is solid, skip all other checks and return just this
					targets.push_back({ box.origin, box.size });
					return targets;
				}

				// Otherwise, define a box for every face needed
				const float thickness = density;

				const Vec3f
					right_origin = box.origin + Vec3f{ box.size.x - thickness, 0.0f, 0.0f },
					back_origin = box.origin + Vec3f{ 0.0f, box.size.y - thickness, 0.0f },
					top_origin = box.origin + Vec3f{ 0.0f, 0.0f, box.size.z - thickness },

					x_size = { thickness, box.size.y, box.size.z },
					y_size = { box.size.x, thickness, box.size.z },
					z_size = { box.size.x, box.size.y, thickness };

				if (box.fillmode.left) targets.push_back({ box.origin, x_size });
				if (box.fillmode.right) targets.push_back({ right_origin, x_size });
				if (box.fillmode.front) targets.push_back({ box.origin, y_size });
				if (box.fillmode.back) targets.push_back({ back_origin, y_size });
				if (box.fillmode.bottom) targets.push_back({ box.origin, z_size });
				if (box.fillmode.top) targets.push_back({ top_origin, z_size });
				return targets;
			} ();

			for (const auto& target : targets) {
				const Vec3f
					fluid_offset = 0.4f * Vec3f{ density, density, density },
					boundary_offset = fluid_offset;
				remove_particles(
					{ target.origin - fluid_offset, target.size + 2.0f * fluid_offset },
					fluid_particles);
				remove_particles(
					{ target.origin - boundary_offset, target.size + 2.0f * boundary_offset },
					boundary_particles);
				const auto target_end = target.origin + target.size;
				for (float x = target.origin.x; x <= target_end.x; x += density) {
					for (float y = target.origin.y; y <= target_end.y; y += density) {
						for (float z = target.origin.z; z <= target_end.z; z += density) {
							Particle p;
							p.position = { x, y, z };
							p.density = case_def.rhop0;
							target_container.emplace_back(p);
						}
					}
				}
			}
		}

	public:
		ParticleGenerator(CaseDef case_def) :
			case_def(case_def)
		{
			for (const auto& model : case_def.poly_data_models)
				fillPLY(model, fluid_particles);

			for (const auto& box : case_def.particle_boxes) {
				using Type = CaseDef::CaseDefBox::Type;
				switch (box.type) {
				case Type::Fluid:
					fillbox(box, fluid_particles);
					break;
				case Type::Boundary:
					fillbox(box, boundary_particles);
					break;
				case Type::Void:
					remove_particles({ box.origin, box.size }, fluid_particles);
					remove_particles({ box.origin, box.size }, boundary_particles);
					break;
				}
			}
		}

		ParticleContainer get_fluid_particles() const {
			return fluid_particles;
		}

		ParticleContainer get_boundary_particles() const {
			return boundary_particles;
		}
};


ParticleContainer ParticleSystem::generate_particles() {
	ParticleGenerator generator(case_def);
	ParticleContainer
		fluid_particles = generator.get_fluid_particles(),
		boundary_particles = generator.get_boundary_particles();
	
	// Make note of the number of fluid particles, then concatenate the two vectors
	// [0, num_of_fluid_particles)   -> fluid
	// [num_of_fluid_particles, end) -> boundary
	num_of_fluid_particles = fluid_particles.size();

	ParticleContainer generated_particles = std::move(fluid_particles);
	generated_particles.insert(generated_particles.end(), boundary_particles.begin(), boundary_particles.end());
	generated_particles.shrink_to_fit();

	return generated_particles;
}

void ParticleSystem::generate_mass_spring_damper() {
	const auto &density = case_def.particles.density;
	MassSpringDamper::damping_coef = case_def.spring.damping;

	for (const auto &mass_spring_system : case_def.spring.mass_spring_systems) {
		MassSpringSystem system;

		system.duration_of_melting = mass_spring_system.duration_of_melting;
		system.initial_k = mass_spring_system.initial_k;
		system.start_of_melting = mass_spring_system.start_of_melting;

		fluid_spring_systems.emplace_back(system);
		boundary_spring_systems.emplace_back(system);
	}

	for (int i = 0; i < num_of_fluid_particles; ++i) {
		const int system_index = [this](const Vec3f &position) {
			for (int i = 0; i < case_def.spring.mass_spring_systems.size(); ++i) {
				if (case_def.spring.mass_spring_systems[i].region.contains(position))
					return i;
			}
			return -1;
		}(particles[i].position);

		if (system_index == -1)
			continue;

		const auto generate_springs = [&](auto& system, const auto& neighbor_indices) {
			for (const auto& index_pair : neighbor_indices) {
				if (index_pair.second <= i)
					continue;
				for (int j = index_pair.first; j < index_pair.second; ++j) {
					if (j <= i)
						continue;
					const float distance = (particles[i].position - particles[j].position).length();
					if (distance <= case_def.spring.max_length) {
						MassSpringDamper msp;
						msp.particle_indices = { i, j };
						msp.resting_length = distance;
						system.springs.emplace_back(msp);
					}
				}
			}
		};

		generate_springs(fluid_spring_systems[system_index], get_fluid_neighbors(particles[i].position));
		generate_springs(boundary_spring_systems[system_index], get_boundary_neighbors(particles[i].position));
	}
}

float ParticleSystem::calculate_time_step() const {
	const float sigma_max = std::transform_reduce(
		std::execution::par_unseq, get_fluid_begin(), get_fluid_end(), std::numeric_limits<float>::lowest(),
		[](const float sigma1, const float sigma2) { return std::max(sigma1, sigma2); },
		[this](const Particle& Pi) {
			const float
				&h = case_def.h;
			const auto index_ranges = get_all_neighbors(Pi.position);

			float sigma = std::numeric_limits<float>::lowest();
			for (const auto index_pair : index_ranges) {
				for (int j = index_pair.first; j < index_pair.second; ++j) {
					const Particle& Pj = particles[j];
					const Vec3f
						v_ij = Pi.velocity - Pj.velocity,
						r_ij = Pi.position - Pj.position;

					const float
						r2 = r_ij.length_squared();

					sigma = std::max(
						sigma,
						std::abs((h * dot_product(v_ij, r_ij)) / (r2 + 0.01f * h * h)));
				}
			}
			return sigma;
		});

	const float
		&h = case_def.h, c = case_def.speedsound,
		dt_cv = h / (c + sigma_max),
		acceleration2_max = std::transform_reduce(
			std::execution::par_unseq, acceleration.begin(), acceleration.end(),
			std::numeric_limits<float>::lowest(),
			[](const float a1, const float a2) { return std::max(a1, a2); },
			[](const Vec3f& accel) {return accel.length_squared(); }),
		dt_f = std::sqrt(h / std::sqrt(acceleration2_max * case_def.particles.mass));

	return case_def.cflnumber * std::min(dt_cv, dt_f);
}

template <ParticleType TypeOfPj>
void ParticleSystem::spring_forces() {
	auto& spring_systems = [this]() -> std::vector<MassSpringSystem>& {
		if (TypeOfPj == ParticleType::Fluid)
			return fluid_spring_systems;
		else
			return boundary_spring_systems;
	}();

	// Remove springs for parts that have melted
	spring_systems.erase(std::remove_if(spring_systems.begin(), spring_systems.end(), [this](const auto &system) {
		return simulation_time >= system.start_of_melting + system.duration_of_melting;
	}), spring_systems.end());

	for (const auto& spring_system : spring_systems) {
		const float k = [&]() {
			if (simulation_time > spring_system.start_of_melting) {
				const float
					time_from_melting_end = (spring_system.start_of_melting + spring_system.duration_of_melting) - simulation_time,
					normalization_coef =
					spring_system.initial_k / (spring_system.duration_of_melting * spring_system.duration_of_melting);
				return normalization_coef * (time_from_melting_end * time_from_melting_end);
			}
			else
				return spring_system.initial_k;
		}();

		for (const auto &spring : spring_system.springs) {
			const auto
				i = spring.particle_indices.first,
				j = spring.particle_indices.second;

			const auto force = spring.compute_force(particles[i], particles[j], k);

			acceleration[i] += force / case_def.particles.mass;
			if (TypeOfPj == ParticleType::Fluid)
				acceleration[j] -= force / case_def.particles.mass;
		}
	}
}

template <ParticleType TypeOfPi, ParticleType TypeOfNeighbors>
void ParticleSystem::compute_derivatives(const int i) {
	const Particle& Pi = particles[i];

	const auto index_ranges = [&](){
		if constexpr (TypeOfNeighbors == ParticleType::Fluid)
			return get_fluid_neighbors(Pi.position);
		else
			return get_boundary_neighbors(Pi.position);
	} ();

	for (const auto index_pair : index_ranges) {
		for (int j = index_pair.first; j < index_pair.second; ++j) {
			const Particle& Pj = particles[j];

			const Vec3f
				v_ij = Pi.velocity - Pj.velocity,
				r_ij = Pi.position - Pj.position;
			const float
				&h = case_def.h,
				r2 = r_ij.length_squared();                           // Squared distance

			if ((TypeOfPi == TypeOfNeighbors && i == j)
				|| r2 == 0.0f || r2 > 4.0f * h * h || v_ij == Vec3f{0.0f, 0.0f, 0.0f})
				continue;

			const float r = std::sqrt(r2);

			const Vec3f
				kernel_gradient_rij = cubic_spline.gradient(r_ij, r);

			density_derivative[i] += case_def.particles.mass * dot_product(v_ij, kernel_gradient_rij);

			// If this is a boundary particle skip the acceleration part
			if constexpr (TypeOfPi == ParticleType::Boundary)
				continue;

			const float
				vel_pos_dot_product = dot_product(v_ij, r_ij),
				&Pi_pressure = pressure[i],
				&Pj_pressure = pressure[j],
				pressure_sum = Pj_pressure + Pi_pressure,
				density_product = Pi.density * Pj.density,
				pi_ij = [&]() {
					if (vel_pos_dot_product < 0.0f) {
						const float &a = case_def.alpha;
						const float
							rho_ij = (Pi.density + Pj.density) / 2.0f,    // Mean density
							mu = (h * vel_pos_dot_product) / (r2 + 0.01f * h * h);

						return -(a * case_def.speedsound * mu) / rho_ij;
					}
					else
						return 0.0f;
				}(),
				tensile_correction_term = [&]() {
					const float
						f_ij = cubic_spline(r) * case_def.tensile_coef,
						coef_i = Pi_pressure > 0 ? 0.01f : -0.2f,
						coef_j = Pj_pressure > 0 ? 0.01f : -0.2f,
						tensile_i = coef_i * (Pi_pressure / (Pi.density * Pi.density)),
						tensile_j = coef_j * (Pj_pressure / (Pj.density * Pj.density));

					return std::pow(f_ij, 4) * (tensile_i + tensile_j);
				} ();

			acceleration[i] -=
			case_def.particles.mass *
			((pressure_sum / density_product) + pi_ij + tensile_correction_term) *
			kernel_gradient_rij;
		}
	}
}

void ParticleSystem::compute_derivatives() {
	#pragma omp parallel for
	for (int i = 0; i < particles.size(); ++i) {
		constexpr int gamma = 7;
		const float
			beta = (case_def.speedsound * case_def.speedsound * case_def.rhop0) / gamma;

		pressure[i] = beta * (std::pow(particles[i].density / case_def.rhop0, gamma) - 1);
	}

	// Fluid-Fluid and Fluid-Boundary interactions
	#pragma omp parallel for
	for (int i = 0; i < num_of_fluid_particles; ++i) {
		acceleration[i] = case_def.gravity;
		density_derivative[i] = 0.0f;
		compute_derivatives<ParticleType::Fluid, ParticleType::Fluid>(i);
		compute_derivatives<ParticleType::Fluid, ParticleType::Boundary>(i);
	}

	// Boundary-Fluid and Boundary-Boundary interactions
	#pragma omp parallel for
	for (int i = num_of_fluid_particles; i < particles.size(); ++i) {
		density_derivative[i] = 0.0f;
		compute_derivatives<ParticleType::Boundary, ParticleType::Fluid>(i);
	}

	if (case_def.spring.on) {
		spring_forces<ParticleType::Fluid>();
		spring_forces<ParticleType::Boundary>();
	}
}

void ParticleSystem::integrate_verlet(const float dt) {
	// How often to use the Verlet corrective step
	// TODO: consider making this variable
	constexpr int corrective_step_interval = 40;

	// To conserve memory, calculate the next state inside prev_particles and then swap the two vectors

	#pragma omp parallel
	if (verlet_step % corrective_step_interval == 0) {
		// Corrective step
		// Fluid particles
		#pragma omp for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto &Pi = particles[i];
			auto &Pi_prev = prev_particles[i];

			Pi_prev.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_prev.velocity = Pi.velocity + dt * acceleration[i];
			Pi_prev.density = Pi.density + dt * density_derivative[i];
		}

		// Boundary particles
		#pragma omp for
		for (int i = num_of_fluid_particles; i < particles.size(); ++i)
			prev_particles[i].density = particles[i].density + dt * density_derivative[i];
	}
	else {
		// Predictor step
		#pragma omp for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto
				&Pi = particles[i],
				&Pi_prev = prev_particles[i];

			Pi_prev.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_prev.velocity = Pi_prev.velocity + 2 * dt * acceleration[i];
			Pi_prev.density = Pi_prev.density + 2 * dt * density_derivative[i];
		}
		#pragma omp for
		for (int i = num_of_fluid_particles; i < particles.size(); ++i)
			prev_particles[i].density = prev_particles[i].density + 2 * dt * density_derivative[i];
	}

	prev_particles.swap(particles);
	++verlet_step;
}

ParticleSystem::ParticleSystem(const CaseDef &case_def) :
	case_def(case_def),
	bounding_box({ case_def.particles.point_min, case_def.particles.point_max - case_def.particles.point_min }),
	search_grid_fluid(bounding_box.origin, bounding_box.origin + bounding_box.size, case_def.h * 2),
	search_grid_boundary(bounding_box.origin, bounding_box.origin + bounding_box.size, case_def.h * 2),
	cubic_spline(case_def.h)
{
	// We need to ensure that all particle arrays are the same size. We also
	// want the boundary positions to be set here so we don't have to copy
	// them every time. Simply copying the whole vector is fast enough.
	particles = generate_particles();
	prev_particles = particles;

	const auto
		fluid_bounding_box = get_particle_axis_aligned_bounding_box(get_fluid_begin(), get_fluid_end()),
		boundary_bounding_box =	get_particle_axis_aligned_bounding_box(get_boundary_begin(), get_boundary_end());

	search_grid_fluid.set_point_min(fluid_bounding_box.origin);
	search_grid_fluid.set_point_max(fluid_bounding_box.origin + fluid_bounding_box.size);

	search_grid_boundary.set_point_min(boundary_bounding_box.origin);
	search_grid_boundary.set_point_max(boundary_bounding_box.origin + boundary_bounding_box.size);

	search_grid_fluid.sort_containers(
		particles.begin(),
		particles.begin() + num_of_fluid_particles,
		prev_particles.begin(),
		nullptr
	);

	search_grid_boundary.sort_containers(
		particles.begin() + num_of_fluid_particles,
		particles.end(),
		prev_particles.begin() + num_of_fluid_particles,
		nullptr
	);

	if (case_def.spring.on)
		generate_mass_spring_damper();

	allocate_memory_for_verlet_variables();
}

ParticleSystem::ParticleSystem(const CaseDef & case_def, State state) :
	case_def(case_def),
	bounding_box({ case_def.particles.point_min, case_def.particles.point_max - case_def.particles.point_min }),
	search_grid_fluid(bounding_box.origin, bounding_box.origin + bounding_box.size, case_def.h * 2),
	search_grid_boundary(bounding_box.origin, bounding_box.origin + bounding_box.size, case_def.h * 2),
	cubic_spline(case_def.h)
{
	num_of_fluid_particles = state.fluid_particles.size();
	particles = std::move(state.fluid_particles);
	particles.insert(particles.end(), state.boundary_particles.begin(), state.boundary_particles.end());
	particles.shrink_to_fit();

	prev_particles = std::move(state.prev_fluid_particles);
	prev_particles.insert(prev_particles.end(), state.prev_boundary_particles.begin(), state.prev_boundary_particles.end());
	prev_particles.shrink_to_fit();

	const auto boundary_bounding_box = get_particle_axis_aligned_bounding_box(get_boundary_begin(), get_boundary_end());

	search_grid_boundary.set_point_min(boundary_bounding_box.origin);
	search_grid_boundary.set_point_max(boundary_bounding_box.origin + boundary_bounding_box.size);

	search_grid_boundary.build_index(particles.begin() + num_of_fluid_particles, particles.end());

	allocate_memory_for_verlet_variables();

	fluid_spring_systems = std::move(state.fluid_spring_systems);
	boundary_spring_systems = std::move(state.boundary_spring_systems);

	MassSpringDamper::damping_coef = case_def.spring.damping;

	verlet_step = state.verlet_step;
	simulation_time = state.simulation_time;
}

ParticleSystem::State ParticleSystem::get_current_state() const {
	State current_state;
	std::copy(particles.begin(), particles.begin() + num_of_fluid_particles, std::back_inserter(current_state.fluid_particles));
	std::copy(particles.begin() + num_of_fluid_particles, particles.end(), std::back_inserter(current_state.boundary_particles));

	std::copy(prev_particles.begin(), prev_particles.begin() + num_of_fluid_particles,
	          std::back_inserter(current_state.prev_fluid_particles));
	std::copy(prev_particles.begin() + num_of_fluid_particles, prev_particles.end(),
	          std::back_inserter(current_state.prev_boundary_particles));

	current_state.fluid_spring_systems = fluid_spring_systems;
	current_state.boundary_spring_systems = boundary_spring_systems;
	current_state.simulation_time = simulation_time;
	current_state.verlet_step = verlet_step;

	return current_state;
}

SearchGrid::cell_indices_container ParticleSystem::get_fluid_neighbors(const Vec3f & position) const
{
	return search_grid_fluid.get_neighbor_indices(position);
}

SearchGrid::cell_indices_container ParticleSystem::get_boundary_neighbors(const Vec3f & position) const
{
	auto neighbors = search_grid_boundary.get_neighbor_indices(position);

	std::for_each(neighbors.begin(), neighbors.end(), [this](auto &neighbor) {
		neighbor.first += num_of_fluid_particles;
		neighbor.second += num_of_fluid_particles;
	});

	return neighbors;
}

SearchGrid::cell_indices_container ParticleSystem::get_all_neighbors(const Vec3f &position) const {
	auto neighbors = get_fluid_neighbors(position);
	const auto boundary_start = neighbors.size();
	search_grid_boundary.get_neighbor_indices(position, neighbors);

	std::for_each(neighbors.begin() + boundary_start, neighbors.end(), [this](auto &neighbor) {
		neighbor.first += num_of_fluid_particles;
		neighbor.second += num_of_fluid_particles;
	});

	return neighbors;
}

void ParticleSystem::simulation_step() {
	// Readjust the search space to only include the bounding box of particles in the simulation
	// Make sure that the search space is not larger than the requested simulation space.
	auto
		fluid_aabb = get_particle_axis_aligned_bounding_box(get_fluid_begin(), get_fluid_end());

	if (!bounding_box.contains(fluid_aabb.origin))
		fluid_aabb.origin = bounding_box.origin;
	if (!bounding_box.contains(fluid_aabb.origin + fluid_aabb.size))
		fluid_aabb.size = (bounding_box.origin + bounding_box.size) - fluid_aabb.origin;
	
	search_grid_fluid.set_point_min(fluid_aabb.origin);
	search_grid_fluid.set_point_max(fluid_aabb.origin + fluid_aabb.size);

	// Sort fluid and boundary particles separately
	search_grid_fluid.sort_containers(
		particles.begin(), particles.begin() + num_of_fluid_particles,
		prev_particles.begin(),
		&fluid_spring_systems,
		&boundary_spring_systems
	);

	compute_derivatives();

	const float time_step = calculate_time_step();
	simulation_time += time_step;
	integrate_verlet(time_step);
}