#ifndef CGOGN_SIMULATION_SIMULATION_SOLVER_H
#define CGOGN_SIMULATION_SIMULATION_SOLVER_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class Simulation_solver
{
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MAP>::Vertex;

public:
	Simulation_constraint<MAP>* constraint_;
	std::shared_ptr<Attribute<Vec3>> speed_;
	std::shared_ptr<Attribute<Vec3>> forces_ext_;
	std::shared_ptr<Attribute<bool>> fixed_vertex;
	std::shared_ptr<Attribute<std::array<Vec3, 9>>> RK_coeff;
	Vec3 gravity_;

	Simulation_solver()
		: constraint_(nullptr), speed_(nullptr), forces_ext_(nullptr), fixed_vertex(nullptr), RK_coeff(nullptr),
		  gravity_(0, 0, 0)
	{
	}

	void init_solver(MAP& m, Simulation_constraint<MAP>* sc, const std::shared_ptr<Attribute<Vec3>>& speed = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& forces = nullptr)
	{
		if (speed != nullptr)
			speed_ = speed;
		if (speed_ == nullptr)
		{
			speed_ = get_attribute<Vec3, Vertex>(m, "simulation_solver_vitesse");
			if (speed_ == nullptr)
			{
				speed_ = add_attribute<Vec3, Vertex>(m, "simulation_solver_vitesse");
				parallel_foreach_cell(m, [&](Vertex v) -> bool {
					value<Vec3>(m, speed_.get(), v) = Vec3(0, 0, 0);
					return true;
				});
			}
		}
		if (forces != nullptr)
			forces_ext_ = forces;
		if (forces_ext_ == nullptr)
		{
			forces_ext_ = get_attribute<Vec3, Vertex>(m, "simulation_solver_forces_ext");
			if (forces_ext_ == nullptr)
			{
				forces_ext_ = add_attribute<Vec3, Vertex>(m, "simulation_solver_forces_ext");
				parallel_foreach_cell(m, [&](Vertex v) -> bool {
					value<Vec3>(m, forces_ext_.get(), v) = Vec3(0, 0, 0);
					return true;
				});
			}
		}

		RK_coeff = get_attribute<std::array<Vec3, 9>, Vertex>(m, "Solver_RK_coeff");
		if (RK_coeff == nullptr)
			RK_coeff = add_attribute<std::array<Vec3, 9>, Vertex>(m, "Solver_RK_coeff");
		constraint_ = sc;
	}

	void compute_time_step(MAP& m, Attribute<Vec3>* vertex_position, Attribute<double>* masse, double time_step)
	{
		if (!constraint_)
			return;
		/*constraint_->solve_constraint(m, vertex_position, forces_ext_.get(), time_step);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (fixed_vertex && value<bool>(m, fixed_vertex.get(), v))
				return true;
			// compute speed
			value<Vec3>(m, speed_.get(), v) =
				0.995 * value<Vec3>(m, speed_.get(), v) +
				time_step * value<Vec3>(m, forces_ext_.get(), v) / value<double>(m, masse, v);
			value<Vec3>(m, vertex_position, v) =
				value<Vec3>(m, vertex_position, v) + time_step * value<Vec3>(m, speed_, v);
			value<Vec3>(m, forces_ext_.get(), v) = Vec3(0, 0, 0);
			return true;
		});*/

		constraint_->solve_constraint(m, vertex_position, this->forces_ext_.get(), time_step);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
				return true;
			// k1
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[0] =
				(time_step * value<Vec3>(m, this->forces_ext_.get(), v) / value<double>(m, masse, v)) +
				time_step * gravity_;
			// j1
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[1] =
				time_step * (0.995 * value<Vec3>(m, this->speed_.get(), v) +
							 value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[0]);

			value<Vec3>(m, vertex_position, v) += value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[1] / 2.0f;

			value<Vec3>(m, this->forces_ext_.get(), v) = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[8];
			return true;
		});

		constraint_->solve_constraint(m, vertex_position, this->forces_ext_.get(), time_step);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
				return true;

			// k2
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[2] =
				(time_step * value<Vec3>(m, this->forces_ext_.get(), v) / value<double>(m, masse, v)) +
				time_step * gravity_;
			// j2
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[3] =
				time_step * (0.995 * value<Vec3>(m, this->speed_.get(), v) +
							 value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[2] / 2.0f);

			value<Vec3>(m, vertex_position, v) -= value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[1] / 2.0f;

			value<Vec3>(m, vertex_position, v) += value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[3] / 2.0f;

			value<Vec3>(m, this->forces_ext_.get(), v) =
				value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[8] + gravity_ * value<double>(m, masse, v);
			return true;
		});

		constraint_->solve_constraint(m, vertex_position, this->forces_ext_.get(), time_step);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
				return true;

			// k3
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[4] =
				(time_step * value<Vec3>(m, this->forces_ext_.get(), v) / value<double>(m, masse, v)) +
				time_step * gravity_;
			// j3
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[5] =
				time_step * (0.995 * value<Vec3>(m, this->speed_.get(), v) +
							 value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[4] / 2.0f);

			value<Vec3>(m, vertex_position, v) -= value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[3] / 2.0f;

			value<Vec3>(m, vertex_position, v) += value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[5];

			value<Vec3>(m, this->forces_ext_.get(), v) =
				value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[8] + gravity_ * value<double>(m, masse, v);
			return true;
		});

		constraint_->solve_constraint(m, vertex_position, this->forces_ext_.get(), time_step);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
				return true;

			// k4
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[6] =
				(time_step * value<Vec3>(m, this->forces_ext_.get(), v) / value<double>(m, masse, v)) +
				time_step * gravity_;
			// j4
			value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[7] =
				time_step * (0.995 * value<Vec3>(m, this->speed_.get(), v) +
							 value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[6]);

			value<Vec3>(m, vertex_position, v) -= value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[5];

			Vec3 k1 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[0];
			Vec3 k2 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[2];
			Vec3 k3 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[4];
			Vec3 k4 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[6];

			Vec3 j1 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[1];
			Vec3 j2 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[3];
			Vec3 j3 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[5];
			Vec3 j4 = value<std::array<Vec3, 9>>(m, this->RK_coeff.get(), v)[7];

			Vec3 diff_pos = 1. / 6. * (j1 + 2 * j2 + 2 * j3 + j4);

			Vec3 diff_speed = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);

			if (diff_pos.norm() < 1.0e-10)
				diff_pos = Vec3(0, 0, 0);

			if (diff_speed.norm() < 1.0e-10)
				diff_speed = Vec3(0, 0, 0);

			value<Vec3>(m, vertex_position, v) += diff_pos;

			value<Vec3>(m, this->speed_.get(), v) = 0.995 * value<Vec3>(m, this->speed_.get(), v) + diff_speed;
			value<Vec3>(m, forces_ext_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_SOLVER_H
