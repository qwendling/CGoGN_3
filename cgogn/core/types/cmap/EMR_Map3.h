#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAP3_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAP3_H_

#include <cgogn/core/types/cmap/EMR_Map2.h>
#include <cgogn/core/types/cmap/cmap3.h>

namespace cgogn
{
template <typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map3_T : public EMR_Map2_T<CMAP>
{
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;
	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Cell<PHI2>;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;

	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi3_;

	EMR_Map3_T(CMAP& m) : EMR_Map2_T<CMAP>(m)
	{
		MR_phi3_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi3_->push_back(this->m_.phi3_);
	}

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	Dart edge_youngest_dart(Dart d) const;

	/***************************************************
	 *                  FACE INFO                      *
	 ***************************************************/

	Dart face_youngest_dart(Dart d) const;

	/***************************************************
	 *                 VOLUME INFO                     *
	 ***************************************************/
	Dart volume_youngest_dart(Dart d) const;
};

using EMR_Map3 = EMR_Map3_T<CMap3>;
} // namespace cgogn

#endif
