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

Vec3f CubicSplinePrecalculated::gradient(const Vec3f & r) const {
	return gradient_values[r.length() / step] * r;
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
					fluid_offset = 0.4f * Vec3f{ case_def.h, case_def.h, case_def.h },
					boundary_offset = 0.4f * Vec3f{ density, density, density };
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
	// [0, num_of_particles)   -> fluid
	// [num_of_particles, end) -> boundary
	num_of_fluid_particles = fluid_particles.size();

	ParticleContainer generated_particles = std::move(fluid_particles);
	generated_particles.insert(generated_particles.end(), boundary_particles.begin(), boundary_particles.end());
	generated_particles.shrink_to_fit();

	return generated_particles;
}

void ParticleSystem::generate_mass_spring_damper() {
	const auto &density = case_def.particles.density;
	MassSpringDamper::k = case_def.spring.stiffness;
	MassSpringDamper::damping_coef = case_def.spring.damping;

	// Fluid - fluid springs
	for (int i = 0; i < num_of_fluid_particles; ++i) {
		const auto neighbor_indices = search_grid_fluid.get_neighbor_indices(particles[i].position);
		for (const auto& index_pair : neighbor_indices) {
			for (int j = index_pair.first; j < index_pair.second; ++j) {
				const float distance = (particles[i].position - particles[j].position).length();
				if (distance <= case_def.spring.max_length && distance > 0) {
					MassSpringDamper msp;
					msp.particle_indices = { i, j };
					msp.resting_length = distance;
					mass_spring_damper.emplace_back(msp);
				}
			}
		}
	}

	num_of_fluid_fluid_springs = mass_spring_damper.size();

	// Fluid - boundary springs
	// Create these separately so they are at the end of the vector
	// This simplifies handling them when calculating acceleration
	for (int i = 0; i < num_of_fluid_particles; ++i) {
		const auto neighbor_indices = search_grid_boundary.get_neighbor_indices(particles[i].position);
		for (const auto& index_pair : neighbor_indices) {
			for (int j = index_pair.first; j < index_pair.second; ++j) {
				const float distance = (particles[i].position - particles[j].position).length();
				if (distance <= case_def.spring.max_length && distance > 0) {
					MassSpringDamper msp;
					msp.particle_indices = { i, j };
					msp.resting_length = distance;
					mass_spring_damper.emplace_back(msp);
				}
			}
		}
	}

	mass_spring_damper.shrink_to_fit();
}

void ParticleSystem::generate_friction_boxes() {
	if (case_def.friction_coef == 0.0f)
		return;

	for (const auto& box : case_def.particle_boxes) {
		if (box.type != CaseDef::CaseDefBox::Type::Boundary)
			continue;

		const float
			wall_thickness = case_def.particles.density,
			friction_thickness = 3 * case_def.h;

		const Vec3f
			left_origin = box.origin + wall_thickness * Vec3f::x_unit(),
			right_origin = box.origin + box.size.x - (wall_thickness + friction_thickness) * Vec3f::x_unit(),
			front_origin = box.origin + wall_thickness * Vec3f::y_unit(),
			back_origin = box.origin + box.size.y - (wall_thickness + friction_thickness) * Vec3f::y_unit(),
			bottom_origin = box.origin + wall_thickness * Vec3f::z_unit(),

			x_size = { friction_thickness, box.size.y,         box.size.z },
			y_size = { box.size.x,         friction_thickness, box.size.z },
			z_size = { box.size.x,         box.size.y,         friction_thickness };

		using Plane = FrictionBox::Plane;
		if (box.fillmode.left) friction_boxes.push_back({ { left_origin, x_size }, Plane::YZ });
		if (box.fillmode.right) friction_boxes.push_back({{ right_origin, x_size }, Plane::YZ });
		if (box.fillmode.front) friction_boxes.push_back({ { front_origin, y_size }, Plane::XZ });
		if (box.fillmode.back) friction_boxes.push_back({{ back_origin, y_size }, Plane::XZ });
		if (box.fillmode.bottom) friction_boxes.push_back({ { bottom_origin, z_size }, Plane::XY });
	}

	friction_boxes.shrink_to_fit();
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

template <ParticleType TypeOfPi, ParticleType TypeOfNeighbors>
void ParticleSystem::compute_derivatives(const int i) {
	const Particle& Pi = particles[i];

	const auto index_ranges = [&](){
		if constexpr (TypeOfNeighbors == ParticleType::Fluid)
			return search_grid_fluid.get_neighbor_indices(Pi.position);
		else
			return search_grid_boundary.get_neighbor_indices(Pi.position);
	} ();

	for (const auto index_pair : index_ranges) {
		for (int j = index_pair.first; j < index_pair.second; ++j) {
			if constexpr (TypeOfNeighbors == ParticleType::Boundary)
				j += num_of_fluid_particles;
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

			const Vec3f
				kernel_gradient_rij = cubic_spline.gradient(r_ij);

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
						f_ij = cubic_spline(r_ij) * case_def.tensile_coef,
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

	// Friction forces
	if constexpr (TypeOfPi == ParticleType::Fluid && TypeOfNeighbors == ParticleType::Boundary) {
		for (const auto& friction_box : friction_boxes) {
			if (friction_box.box.contains(Pi.position)) {
				const auto tangent_velocity = [&]() {
					switch (friction_box.plane) {
					case FrictionBox::Plane::XY:
						return Vec3f{ Pi.velocity.x, Pi.velocity.y, 0.0f };
					case FrictionBox::Plane::XZ:
						return Vec3f{ Pi.velocity.x, 0.0f, Pi.velocity.z };
					case FrictionBox::Plane::YZ:
						return Vec3f{ 0.0f, Pi.velocity.y, Pi.velocity.z };
					}
				} ();

				acceleration[i] -=
					(case_def.particles.mass * case_def.friction_coef * tangent_velocity) / particles[i].density;
			}
		}
	}
}

void ParticleSystem::compute_derivatives() {
	#pragma omp parallel for
	for (int i = 0; size_t(i) < particles.size(); ++i) {
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
	for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i) {
		density_derivative[i] = 0.0f;
		compute_derivatives<ParticleType::Boundary, ParticleType::Fluid>(i);
	}

	if (!case_def.spring.on)
		return;

	// Mass-Spring system forces
	if (simulation_time > case_def.spring.start_of_stiffness_change) {
		const float duration_of_change = simulation_time - case_def.spring.start_of_stiffness_change;
		MassSpringDamper::k =
			case_def.spring.stiffness + duration_of_change * case_def.spring.rate_of_stiffness_change;
		if (MassSpringDamper::k <= 0.0f) {
			mass_spring_damper.clear();
			case_def.spring.on = false;
		}
	}

	for (int k = 0; k < num_of_fluid_fluid_springs; ++k) {
		const auto& spring = mass_spring_damper[k];
		const auto
			i = spring.particle_indices.first,
			j = spring.particle_indices.second;

		const auto force = spring.compute_force(particles[i], particles[j]);

		acceleration[i] += force / case_def.particles.mass;
		acceleration[j] -= force / case_def.particles.mass;
	}

	for (int k = num_of_fluid_fluid_springs; k < mass_spring_damper.size(); ++k) {
		const auto& spring = mass_spring_damper[k];
		const auto
			i = spring.particle_indices.first,
			j = spring.particle_indices.second;

		const auto force = spring.compute_force(particles[i], particles[j]);

		acceleration[i] += force / case_def.particles.mass;
	}
}

void ParticleSystem::integrate_verlet(const float dt) {
	// How often to use the Verlet corrective step
	// TODO: consider making this variable
	constexpr int corrective_step_interval = 40;

	#pragma omp parallel
	if (verlet_step % corrective_step_interval == 0) {
		// Corrective step
		// Fluid particles
		#pragma omp for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto &Pi = particles[i];
			auto &Pi_next = next_particles[i];

			Pi_next.velocity = Pi.velocity + dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi.density + dt * density_derivative[i];
		}

		// Boundary particles
		#pragma omp for
		for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i)
			next_particles[i].density = particles[i].density + dt * density_derivative[i];
	}
	else {
		// Predictor step
		#pragma omp for
		for (int i = 0; i < num_of_fluid_particles; ++i) {
			auto
				&Pi = particles[i],
				&Pi_next = next_particles[i],
				&Pi_prev = prev_particles[i];

			Pi_next.velocity = Pi_prev.velocity + 2 * dt * acceleration[i];
			Pi_next.position = Pi.position + dt * Pi.velocity + 0.5f * dt * dt * acceleration[i];
			Pi_next.density = Pi_prev.density + 2 * dt * density_derivative[i];
		}
		#pragma omp for
		for (int i = num_of_fluid_particles; size_t(i) < particles.size(); ++i)
			next_particles[i].density = prev_particles[i].density + 2 * dt * density_derivative[i];
	}

	++verlet_step;
}

void ParticleSystem::remove_out_of_bounds_particles() {
	// Figure out the particles in current time step that are out of bounds and remove them
	// Remove the same particles and do any necessary re-ordering to prev particles to keep
	// the two containers consistent. Resize next particles for size constistency but don't care
	// about the data because it's supposed to be garbage at the point this function needs to run.
	
	const auto indices_to_remove = [this]() {
		std::vector<size_t> indices_to_remove;
		for (size_t i = 0; i < num_of_fluid_particles; ++i) {
			if (!bounding_box.contains(particles[i].position))
				indices_to_remove.emplace_back(i);
		}
		return indices_to_remove;
	} ();


	// Simply remove all springs that point to particles about to be removed
	auto to_be_removed = [&indices_to_remove](const size_t i) {
		return std::binary_search(indices_to_remove.begin(), indices_to_remove.end(), i);
	};
	mass_spring_damper.erase(
		std::remove_if(std::execution::par_unseq, mass_spring_damper.begin(), mass_spring_damper.end(),
			[to_be_removed](const MassSpringDamper& spring) {
				return
					to_be_removed(spring.particle_indices.first) ||
					to_be_removed(spring.particle_indices.second);
			}),
		mass_spring_damper.end()
	);

	std::for_each(std::execution::par_unseq,
		mass_spring_damper.begin(), mass_spring_damper.end(), [&](auto &spring) {
		auto adjust_index = [&](int &index) {
			auto distance_from_end = num_of_fluid_particles - index;
			while (distance_from_end <= indices_to_remove.size()) {
				index = indices_to_remove[indices_to_remove.size() - distance_from_end];
				distance_from_end = num_of_fluid_particles - index;
			}
		};

		adjust_index(spring.particle_indices.first);
		adjust_index(spring.particle_indices.second);
	});

	// To preserve the contiguous storage of the container with the minimum complexity cost,
	// we will move the out-of-bounds particles to the end of the range and remove them from there
	// This not only minimizes the number of moves required, it also invalidates the minimum amount of
	// indices to particles.

	size_t next_slot_to_swap = 1;
	for (auto i = indices_to_remove.crbegin(); i != indices_to_remove.crend(); ++i) {
		auto particle_temp = std::move(particles[*i]);
		auto prev_particle_temp = std::move(prev_particles[*i]);

		particles[*i] = std::move(particles[num_of_fluid_particles - next_slot_to_swap]);
		prev_particles[*i] = std::move(prev_particles[num_of_fluid_particles - next_slot_to_swap]);

		particles[num_of_fluid_particles - next_slot_to_swap] =
			std::move(particles[particles.size() - next_slot_to_swap]);
		prev_particles[num_of_fluid_particles - next_slot_to_swap] =
			std::move(prev_particles[prev_particles.size() - next_slot_to_swap]);

		particles[particles.size() - next_slot_to_swap] = std::move(particle_temp);
		prev_particles[prev_particles.size() - next_slot_to_swap] = std::move(prev_particle_temp);

		++next_slot_to_swap;
	}

	// Readjust fluid and total sizes
	const size_t
		number_of_removed_particles = indices_to_remove.size(),
		total_number_of_particles = particles.size() - number_of_removed_particles;
	num_of_fluid_particles -= number_of_removed_particles;

	// A simple resize works because all out-of-bounds particles are at the end
	particles.resize(total_number_of_particles);
	prev_particles.resize(total_number_of_particles);
	next_particles.resize(total_number_of_particles);
}

ParticleSystem::ParticleSystem(const CaseDef &case_def) :
	case_def(case_def),
	particles(generate_particles()),
	bounding_box({ case_def.particles.point_min, case_def.particles.point_max - case_def.particles.point_min }),
	search_grid_fluid(bounding_box.origin, bounding_box.origin + bounding_box.size, case_def.h),
	search_grid_boundary(bounding_box.origin, bounding_box.origin + bounding_box.size, case_def.h),
	cubic_spline(case_def.h)
{
	// We need to ensure that all particle arrays are the same size. We also
	// want the boundary positions to be set here so we don't have to copy
	// them every time. Simply copying the whole vector is fast enough.
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

	next_particles = particles;

	if (case_def.spring.on)
		generate_mass_spring_damper();

	generate_friction_boxes();

	allocate_memory_for_verlet_variables();
}

SearchGrid::cell_indices_container ParticleSystem::get_all_neighbors(const Vec3f &position) const {
	SearchGrid::cell_indices_container neighbors;
	neighbors.reserve(54);
	search_grid_fluid.get_neighbor_indices(position, neighbors);
	const auto boundary_start = neighbors.size();
	search_grid_boundary.get_neighbor_indices(position, neighbors);

	// Apply offset to boundary particles
	for (auto i = boundary_start; i < neighbors.size(); ++i) {
		neighbors[i].first += num_of_fluid_particles;
		neighbors[i].second += num_of_fluid_particles;
	}

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
		&mass_spring_damper
	);

	compute_derivatives();

	const float time_step = calculate_time_step();
	simulation_time += time_step;
	integrate_verlet(time_step);

	// Move variables of next step to current and current to prev
	// Using swap instead of move assignment retains the size of next_particles
	// next_particles practically holds garbage values after this
	prev_particles.swap(particles);
	particles.swap(next_particles);

	// Temporarily remove this because it doesn't function correctly with springs attached
	// to boundary particles
	//remove_out_of_bounds_particles();
}