#ifndef CGOGN_SIMULATION_SIMULATION_CONSTRAINT_H
#define CGOGN_SIMULATION_SIMULATION_CONSTRAINT_H

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class Simulation_constraint
{
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MAP>::Vertex;

public:
	virtual void solve_constraint(const MAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces,
								  double time_step) = 0;
	virtual void update_topo(const MAP&, const std::vector<Vertex>&){};
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_CONSTRAINT_H
