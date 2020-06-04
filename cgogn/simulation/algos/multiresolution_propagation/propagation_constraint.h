#ifndef CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_CONSTRAINT_H
#define CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_CONSTRAINT_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{
namespace simulation
{
template <typename MR_MAP>
class Propagation_Constraint
{
	template <typename T>
	using Attribute = typename mesh_traits<MR_MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MR_MAP>::Vertex;

public:
	virtual void propagate(MR_MAP& m_meca, MR_MAP& m_geom, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces,
						   Attribute<Vec3>* pos_relative, Attribute<std::array<Vertex, 3>>* parent,
						   const std::function<void(Vertex)>& integration, double time_step = 0.005f) = 0;
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_MULTIRESOLTION_PROPAGATION_PROPAGATION_CONSTRAINT_H
