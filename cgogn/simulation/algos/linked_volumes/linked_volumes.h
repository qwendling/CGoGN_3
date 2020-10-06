#ifndef CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
#define CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H

#include <cgogn/core/functions/attributes.h>
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

public:
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
