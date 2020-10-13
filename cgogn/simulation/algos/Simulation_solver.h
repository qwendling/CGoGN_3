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

	Simulation_solver() : constraint_(nullptr), speed_(nullptr), forces_ext_(nullptr), fixed_vertex(nullptr)
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
		constraint_ = sc;
	}

	void compute_time_step(MAP& m, Attribute<Vec3>* vertex_position, Attribute<double>* masse, double time_step)
	{
		if (!constraint_)
			return;
		constraint_->solve_constraint(m, vertex_position, forces_ext_.get(), time_step);
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
		});
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_SOLVER_H
