#include <cgogn/core/types/cmap/EMR_Map2.h>
#include <cgogn/core/types/cmap/cmap3.h>

namespace cgogn
{
template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map3_T : public EMR_Map2_T<CMAP>
{
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;
	std::vector<std::shared_ptr<Attribute<Dart>>> MR_phi3_;

	EMR_Map3_T() : EMR_Map2_T<CMAP>()
	{
		MR_phi3_.push_back(CMAP::phi3_);
	}

	virtual void change_resolution_level(uint32 new_level){
		EMR_Map2_T<CMAP>::change_resolution_level(new_level);
		CMAP::phi3_ = MR_phi3_[this->current_level_];
	}
};

using EMR_Map3 = EMR_Map3_T<CMap3>;
} // namespace cgogn
