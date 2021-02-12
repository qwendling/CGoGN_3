#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAP1_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAP1_H_

#include <cgogn/core/types/cmap/EMR_MapBase.h>
#include <cgogn/core/types/cmap/cmap1.h>

namespace cgogn
{

template <typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map1_T : public EMR_MapBase_T<CMAP>
{
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;

	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi1_;
	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi_1_;

	EMR_Map1_T() : EMR_MapBase_T<CMAP>()
	{
		MR_phi1_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi_1_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi1_->push_back(this->phi1_);
		MR_phi_1_->push_back(this->phi_1_);
	}

	virtual void add_resolution()
	{
		EMR_MapBase_T<CMAP>::add_resolution();
		this->phi1_ = (*MR_phi1_)[this->maximum_level_];
		this->phi_1_ = (*MR_phi_1_)[this->maximum_level_];
	}
};

struct EMR_Map1 : EMR_MapBase<EMR_Map1_T<CMap1>>
{

	EMR_Map1(EMR_Map1_T<CMap1>& m) : EMR_MapBase<EMR_Map1_T<CMap1>>(m)
	{
	}
};

} // namespace cgogn

#endif
