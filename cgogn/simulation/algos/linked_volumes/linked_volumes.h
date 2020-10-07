#ifndef CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
#define CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class Linked_volumes
{
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Face = typename mesh_traits<MAP>::Face;
	using Volume = typename mesh_traits<MAP>::Volume;

	MAP& m_;
	std::shared_ptr<Attribute<Vec3>> centroid_;
	std::shared_ptr<Attribute<double>> distance_plan_;

public:
	Linked_volumes(MAP& m) : m_(m)
	{
		centroid_ = get_attribute<Vec3, Volume>(m, "centroid");
		if (centroid_ == nullptr)
			centroid_ = add_attribute<Vec3, Volume>(m, "centroid");
		distance_plan_ = get_attribute<Vec3, Volume>(m, "distance_plan");
		if (distance_plan_ == nullptr)
			distance_plan_ = add_attribute<Vec3, Volume>(m, "distance_plan");
	}

	void compute_cut_plan(Vec3 dir_plan, double w, Attribute<Vec3> pos, bool compute_centroid = true)
	{
		// Ajout plan de coupe en parametre
		if (compute_centroid)
		{
			geometry::compute_centroid<Vec3, Volume>(m_, pos, centroid_.get());
		}
		parallel_foreach_cell(m_, [&](Volume v) -> bool {
			value<double>(m_, this->distance_plan_.get(), v) = dir_plan.dot(value<Vec3>(m_, this->centroid_.get(), v));
			return true;
		});
		std::vector<Face> face_vect;
		foreach_cell(m_, [&](Face f) -> bool {
			if (is_incident_to_boundary(m_, f))
			{
				return true;
			}
			double v1 = value<double>(m_, this->distance_plan_.get(), Volume(f.dart)) - w;
			double v2 = value<double>(m_, this->distance_plan_.get(), Volume(phi3(m_, f.dart))) - w;
			if (v1 * v2 < 0)
			{
				face_vect.push_back(f);
			}
			return true;
		});

		// unsew faces
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
