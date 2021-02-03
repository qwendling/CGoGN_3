#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAP2_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAP2_H_

#include <cgogn/core/types/cmap/EMR_Map1.h>
#include <cgogn/core/types/cmap/cmap2.h>

namespace cgogn
{

template <typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map2_T : public EMR_Map1_T<CMAP>
{
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;

	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi2_;

	EMR_Map2_T() : EMR_Map1_T<CMAP>()
	{
		MR_phi2_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi2_->push_back(this->phi2_);
	}
};

struct EMR_Map2 : EMR_MapBase<EMR_Map2_T<CMap2>>
{

	EMR_Map2(EMR_Map2_T<CMap2>& m) : EMR_MapBase<EMR_Map2_T<CMap2>>(m)
	{
	}
};

} // namespace cgogn

#endif
