#include <cgogn/core/cgogn_core_export.h>
#include <cgogn/core/utils/numerics.h>
#include <vector>

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/container/vector.h>

#include <cgogn/core/types/cmap/cell.h>

#include <any>
#include <array>
#include <condition_variable>
#include <unordered_map>
#include <deque>

namespace cgogn
{

template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map1_T : public CMAP
{
	std::deque<std::shared_ptr<Attribute<Dart>>> MR_phi1_;
	std::deque<std::shared_ptr<Attribute<Dart>>> MR_phi_1_;

	uint32 current_level_;

	EMR_Map1_T() : CMAP()
	{
		MR_phi1_.push_back(phi1_);
		MR_phi_1_.push_back(phi_1_);
	}

	virtual void change_resolution_level(uint32 new_level){
		cgogn_message_assert(0 <= new_level && new_level < MR_phi1_.size(), "Access to an undefined level");
		current_level_ = new_level;
		phi1_ = MR_phi1_[current_level_];
		phi_1_ = MR_phi_1_[current_level_];
	}
};

using EMR_Map1 = EMR_Map1_T<CMap1>;

template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map2_T : public EMR_Map1_T<CMAP>
{
	std::deque<std::shared_ptr<Attribute<Dart>>> MR_phi2_;

	EMR_Map2_T() : EMR_Map1_T<CMAP>()
	{
		MR_phi2_.push_back(phi2_);
	}

	virtual void change_resolution_level(uint32 new_level){
		EMR_Map1_T<CMAP>::change_resolution_level(new_level);
		phi2_ = MR_phi2_[current_level_];
	}
};

using EMR_Map2 = EMR_Map2_T<CMap2>;

template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_Map3_T : public EMR_Map2_T<CMAP>
{
	std::deque<std::shared_ptr<Attribute<Dart>>> MR_phi3_;

	EMR_Map3_T() : EMR_Map2_T<CMAP>()
	{
		MR_phi3_.push_back(phi3_);
	}

	virtual void change_resolution_level(uint32 new_level){
		EMR_Map2_T<CMAP>::change_resolution_level(new_level);
		phi3_ = MR_phi3_[current_level_];
	}
};

using EMR_Map3 = EMR_Map3_T<CMap3>;
} // namespace cgogn
