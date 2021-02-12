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

	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi3_;

	EMR_Map3_T() : EMR_Map2_T<CMAP>()
	{
		MR_phi3_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi3_->push_back(this->phi3_);
	}

	virtual void add_resolution()
	{
		EMR_Map2_T<CMAP>::add_resolution();
		this->phi3_ = (*MR_phi3_)[this->maximum_level_];
	}
};

struct EMR_Map3 : EMR_MapBase<EMR_Map3_T<CMap3>>
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

	EMR_Map3(EMR_Map3_T<CMap3>& m) : EMR_MapBase<EMR_Map3_T<CMap3>>(m)
	{
	}

	CMap3* get_map()
	{
		CMap3* result = static_cast<CMap3&>(this->m_).get_copy();
		result->phi1_ = (*m_.MR_phi1_)[current_level_];
		result->phi_1_ = (*m_.MR_phi_1_)[current_level_];
		result->phi2_ = (*m_.MR_phi2_)[current_level_];
		result->phi3_ = (*m_.MR_phi3_)[current_level_];
		return result;
	}

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	Dart edge_youngest_dart(Dart d) const;
	bool edge_is_subdivided(Dart d) const;
	uint32 edge_level(Dart d) const;

	/***************************************************
	 *                  FACE INFO                      *
	 ***************************************************/

	Dart face_youngest_dart(Dart d) const;
	Dart face_oldest_dart(Dart d) const;
	bool face_is_subdivided(Dart d) const;
	uint32 face_level(Dart d) const;

	/***************************************************
	 *                 VOLUME INFO                     *
	 ***************************************************/
	Dart volume_youngest_dart(Dart d) const;
	Dart volume_oldest_dart(Dart d) const;
	bool volume_is_subdivided(Dart d) const;
	uint32 volume_level(Dart d) const;
};

} // namespace cgogn

#endif
