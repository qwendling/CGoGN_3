#ifndef CGOGN_SIMULATION_LINKED_VOLUMES_CUTTING_TOOLS_H
#define CGOGN_SIMULATION_LINKED_VOLUMES_CUTTING_TOOLS_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{
namespace simulation
{
class Cutting_Tools
{
	uint32 size_;
	using Vec3 = geometry::Vec3;
	using Mat4d = geometry::Mat4d;
	using Mat3d = geometry::Mat3d;
	std::vector<Vec3> vertices_;
	Eigen::Affine3d tr_;
	Eigen::Affine3d tr_inv_;

public:
	Cutting_Tools(uint32 size) : size_(size)
	{
		Eigen::Transform<double, 3, Eigen::Affine> rot(Eigen::AngleAxis(2 * M_PI / double(size), Vec3(0, 0, 1)));
		Eigen::Vector3d v(1, 0, 0);
		for (int i = 0; i < size; i++)
		{
			vertices_.push_back(v);
			v = rot * v;
		}
	}

	uint32 intersect(const Vec3& v)
	{
		Vec3 v2 = tr_inv_ * v;
		uint32 result = 0;
		if (v2(2) > 0)
		{
			result = 1;
		}
		for (int i = 0; i < size_ - 1; i++)
		{
			Vec3 a = vertices_[i + 1] - vertices_[i];
			Vec3 b = v2 - vertices_[i];
			if (a(0) * b(1) - a(1) * b(0) < 0)
			{
				result = 2;
			}
		}
		// retourne 1 si au dessus
		// retourne 0 si en dessous
		// retourne 2 si ni l'un ni l'autre
		return result;
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_LINKED_VOLUMES_CUTTING_TOOLS_H
