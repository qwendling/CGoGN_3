#ifndef CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
#define CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>
#include <cgogn/simulation/algos/Simulation_solver.h>
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_constraint.h>

namespace cgogn
{
namespace simulation
{
template <typename MR_MAP>
class Simulation_solver_multiresolution : public Simulation_solver<MR_MAP>
{
	using Self = Simulation_solver_multiresolution<MR_MAP>;
	using Inherit = Simulation_solver<MR_MAP>;
	template <typename T>
	using Attribute = typename mesh_traits<MR_MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MR_MAP>::Vertex;

public:
	Propagation_Constraint<MR_MAP>* pc_;
	std::shared_ptr<Attribute<std::array<Vertex, 3>>> parents_;
	std::shared_ptr<Attribute<Vec3>> relative_pos_;
	Simulation_solver_multiresolution()
		: Simulation_solver<MR_MAP>(), pc_(nullptr), parents_(nullptr), relative_pos_(nullptr)
	{
	}
	void init_solver(MR_MAP& m, Simulation_constraint<MR_MAP>* sc, Propagation_Constraint<MR_MAP>* pc = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& speed = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& forces = nullptr)
	{
		Inherit::init_solver(m, sc, speed, forces);
		pc_ = pc;
	}

	void reset_forces(MR_MAP& m)
	{
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			value<Vec3>(m, this->speed_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
	}

	void compute_time_step(MR_MAP& m_meca, MR_MAP& m_geom, Attribute<Vec3>* vertex_position, Attribute<double>* masse,
						   double time_step)
	{
		if (!this->constraint_)
			return;
		this->constraint_->solve_constraint(m_meca, vertex_position, this->forces_ext_.get(), time_step);
		parallel_foreach_cell(m_meca, [&](Vertex v) -> bool {
			// compute speed
			value<Vec3>(m_meca, this->speed_.get(), v) =
				0.995 * value<Vec3>(m_meca, this->speed_.get(), v) +
				time_step * value<Vec3>(m_meca, this->forces_ext_.get(), v) / value<double>(m_meca, masse, v);
			value<Vec3>(m_meca, vertex_position, v) =
				value<Vec3>(m_meca, vertex_position, v) + time_step * value<Vec3>(m_meca, this->speed_, v);
			return true;
		});
		if (!pc_)
			return;
		pc_->propagate(
			m_meca, m_geom, vertex_position, this->forces_ext_.get(), relative_pos_.get(), parents_.get(),
			[&](Vertex v) {
				value<Vec3>(m_geom, this->speed_.get(), v) =
					0.995 * value<Vec3>(m_geom, this->speed_.get(), v) +
					time_step * value<Vec3>(m_geom, this->forces_ext_.get(), v) / value<double>(m_geom, masse, v);
				value<Vec3>(m_geom, vertex_position, v) =
					value<Vec3>(m_geom, vertex_position, v) + time_step * value<Vec3>(m_geom, this->speed_, v);
			},
			time_step);
		parallel_foreach_cell(m_geom, [&](Vertex v) -> bool {
			value<Vec3>(m_geom, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
