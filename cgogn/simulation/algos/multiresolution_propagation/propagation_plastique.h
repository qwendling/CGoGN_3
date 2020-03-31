#ifndef CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_PLASTIQUE_H
#define CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_PLASTIQUE_H

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
class Propagation_Plastique : public Propagation_Constraint<MR_MAP>
{
	template <typename T>
	using Attribute = typename mesh_traits<MR_MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MR_MAP>::Vertex;

public:
	double alpha_;

	Propagation_Plastique() : alpha_(0.5f)
	{
	}
	Propagation_Plastique(double alpha) : alpha_(alpha)
	{
	}
	void propagate(MR_MAP& m_meca, MR_MAP& m_geom, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces,
				   Attribute<Vec3>* pos_relative, Attribute<std::pair<Vertex, Vertex>>* parent,
				   const std::function<void(Vertex)>& integration, double time_step) override
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
				if (is_boundary(cph_, v.dart))
					v.dart = phi3(cph_, v.dart);
				std::pair<Vertex, Vertex> p = value<std::pair<Vertex, Vertex>>(m_geom, parent, v);
				Vec3 C1 = (value<Vec3>(m_geom, pos, p.first) + value<Vec3>(m_geom, pos, p.second)) / 2;
				Vec3 C2 = geometry::centroid<Vec3>(cph_, CPH3::Volume(v.dart), pos);
				Vec3 A = value<Vec3>(m_geom, pos, p.first);
				Vec3 V1 = (A - C1).normalized();
				Vec3 V2 = (C2 - C1).normalized();
				Vec3 V3 = (V1.cross(V2));
				V3 = V3.normalized();
				Vec3 p_r = value<Vec3>(m_geom, pos_relative, v);
				Vec3 dest = C1 + V1 * p_r[0] + V2 * p_r[1] + V3 * p_r[2];
				Vec3 cur_pos = value<Vec3>(m_geom, pos, v);
				value<Vec3>(m_geom, result_forces, v) += alpha_ * (dest - cur_pos) / (time_step * time_step);
				integration(v);
				value<Vec3>(m_geom, pos, v) = dest;
			}
		}
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_PLASTIQUE_H
