#include <cgogn/core/types/cmap/EMR_MapBase.h>
#include <cgogn/core/types/cmap/cmap1.h>


namespace cgogn
{

template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map1_T : public EMR_MapBase_T<CMAP>
{
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;

	std::vector<std::shared_ptr<Attribute<Dart>>> MR_phi1_;
	std::vector<std::shared_ptr<Attribute<Dart>>> MR_phi_1_;

	EMR_Map1_T() : EMR_MapBase_T<CMAP>()
	{
		MR_phi1_.push_back(CMAP::phi1_);
		MR_phi_1_.push_back(CMAP::phi_1_);
	}

	virtual void change_resolution_level(uint32 new_level){
		cgogn_message_assert(0 <= new_level && new_level < MR_phi1_.size(), "Access to an undefined level");
		CMAP::phi1_ = MR_phi1_[this->current_level_];
		CMAP::phi_1_ = MR_phi_1_[this->current_level_];
	}
};

using EMR_Map1 = EMR_Map1_T<CMap1>;

} // namespace cgogn
