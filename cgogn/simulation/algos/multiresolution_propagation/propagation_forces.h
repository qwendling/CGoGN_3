#ifndef CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_FORCES_H
#define CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_FORCES_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/simulation/algos/multiresolution_propagation/propagation_constraint.h>

namespace cgogn
{
namespace simulation
{
template <typename MR_MAP>
class Propagation_Forces : public Propagation_Constraint<MR_MAP>
{
	template <typename T>
	using Attribute = typename mesh_traits<MR_MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MR_MAP>::Vertex;

public:
	Propagation_Forces()
	{
	}
	void propagate(MR_MAP& m_meca, MR_MAP& m_geom, Attribute<Vec3>*, Attribute<Vec3>* result_forces, Attribute<Vec3>*,
				   Attribute<std::pair<Vertex, Vertex>>* parent, const std::function<void(Vertex)>& integration,
				   double) override
	{

		std::vector<std::vector<Vertex>> vect_vertex_per_resolution;
		vect_vertex_per_resolution.resize(m_geom.maximum_level_ + 1);
		foreach_cell(m_geom, [&](Vertex v) -> bool {
			if (!m_meca.dart_is_visible(v.dart))
				vect_vertex_per_resolution[m_geom.dart_level(v.dart)].push_back(v);
			return true;
		});
		CPH3 cph_(m_geom);
		for (int i = 0; i < vect_vertex_per_resolution.size(); i++)
		{
			cph_.current_level_ = i;
			for (auto& v : vect_vertex_per_resolution[i])
			{
				std::pair<Vertex, Vertex> p = value<std::pair<Vertex, Vertex>>(m_geom, parent, v);
				value<Vec3>(m_geom, result_forces, v) +=
					(value<Vec3>(m_geom, result_forces, p.first) + value<Vec3>(m_geom, result_forces, p.second)) / 2;
				integration(v);
			}
		}
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_FORCES_H
