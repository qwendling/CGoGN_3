#include <cgogn/core/types/cmap/EMR_Map1.h>
#include <cgogn/core/types/cmap/cmap2.h>


namespace cgogn
{

template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map2_T : public EMR_Map1_T<CMAP>
{
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;

	std::vector<std::shared_ptr<Attribute<Dart>>> MR_phi2_;

	EMR_Map2_T() : EMR_Map1_T<CMAP>()
	{
		MR_phi2_.push_back(CMAP::phi2_);
	}

	virtual void change_resolution_level(uint32 new_level){
		EMR_Map1_T<CMAP>::change_resolution_level(new_level);
		CMAP::phi2_ = MR_phi2_[this->current_level_];
	}
};

using EMR_Map2 = EMR_Map2_T<CMap2>;

} // namespace cgogn
