#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAP3_ADAPTATIVE_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAP3_ADAPTATIVE_H_

#include <cgogn/core/types/cmap/EMR_Map3.h>

namespace cgogn
{

struct EMR_Map3_Adaptative : EMR_Map3
{

	using MAP = CMap3;
	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Cell<PHI2>;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	template <typename T>
	using Attribute = typename CMap3::template Attribute<T>;

	std::shared_ptr<Attribute<uint32>> dart_visibility_;

	EMR_Map3_Adaptative(EMR_Map3_T<CMap3>& m) : EMR_Map3(m)
	{
	}

	bool check_integrity() const;
	uint32 get_dart_visibility(Dart d) const;

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	bool edge_is_subdivided(Dart d) const;
	uint32 edge_level(Dart d) const;

	/***************************************************
	 *                  FACE INFO                      *
	 ***************************************************/

	bool face_is_subdivided(Dart d) const;
	uint32 face_level(Dart d) const;

	/***************************************************
	 *                 VOLUME INFO                     *
	 ***************************************************/
	bool volume_is_subdivided(Dart d) const;
	uint32 volume_level(Dart d) const;
};

} // namespace cgogn

#endif
