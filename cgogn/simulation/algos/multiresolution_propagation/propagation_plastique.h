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
	void propagate(MR_MAP& m_meca, MR_MAP& m_geom, Attribute<Vec3>* pos, Attribute<Vec3>* speed, Attribute<Vec3>*,
				   Attribute<double>*, Attribute<Vec3>* pos_relative, Attribute<std::array<Vertex, 3>>* parent,
				   double time_step) const override
	{

		std::vector<std::vector<Vertex>> vect_vertex_per_resolution;
		vect_vertex_per_resolution.resize(m_geom.maximum_level_ + 1);
		std::unordered_set<uint32> list_vertices;

		std::function<void(Vertex)> add_vertex;
		add_vertex = [&](Vertex v) {
			if (m_meca.dart_is_visible(v.dart))
				return;
			if (list_vertices.count(index_of(m_meca, v)))
				return;
			vect_vertex_per_resolution[m_geom.dart_level(v.dart)].push_back(v);

			list_vertices.insert(index_of(m_meca, v));
			std::array<Vertex, 3> p = value<std::array<Vertex, 3>>(m_geom, parent, v);
			add_vertex(p[0]);
			add_vertex(p[1]);
			add_vertex(p[2]);
		};

		foreach_cell(m_geom, [&](Vertex v) -> bool {
			if (!m_meca.dart_is_visible(v.dart))
			{
				add_vertex(v);
			}
			return true;
		});
		std::clock_t start;
		double duration;

		start = std::clock();
		for (uint32 i = 1; i < vect_vertex_per_resolution.size(); i++)
		{
			for (auto& v : vect_vertex_per_resolution[i])
			{
				std::array<Vertex, 3>& p = value<std::array<Vertex, 3>>(m_geom, parent, v);
				Vec3 A = value<Vec3>(m_geom, pos, p[0]);
				Vec3 B = value<Vec3>(m_geom, pos, p[1]);
				Vec3 C = value<Vec3>(m_geom, pos, p[2]);
				Vec3 V1 = (B - A).normalized();
				Vec3 V2 = (C - A).normalized();
				Vec3 V3 = (V1.cross(V2));
				V3 = V3.normalized();
				Vec3 p_r = value<Vec3>(m_geom, pos_relative, v);
				Vec3 dest = A + V1 * p_r[0] + V2 * p_r[1] + V3 * p_r[2];
				Vec3 cur_pos = value<Vec3>(m_geom, pos, v);
				if (speed != nullptr)
					value<Vec3>(m_geom, speed, v) = (dest - cur_pos) / time_step;
				value<Vec3>(m_geom, pos, v) = dest;
			}
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "Real P time : " << duration << std::endl;
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_PLASTIQUE_H
