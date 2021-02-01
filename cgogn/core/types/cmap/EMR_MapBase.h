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

namespace cgogn
{

template<typename CMAP>
struct CGOGN_CORE_EXPORT EMR_MapBase_T : public CMAP
{

	using MAP = CMAP;
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;

	std::shared_ptr<Attribute<uint32>> dart_level_;

	uint32 current_level_;

	EMR_MapBase_T() : CMAP(),current_level_(0)
	{
		dart_level_ = CMAP::darts_->template get_attribute<uint32>("dart_level");
		if (!dart_level_)
			dart_level_ = CMAP::darts_->template add_attribute<uint32>("dart_level");
	}

	virtual void change_resolution_level(uint32 new_level){
		current_level_ = new_level;
	}

	uint32 dart_level(Dart d) const{
		return (*dart_level_)[d.index];
	}
	void set_dart_level(Dart d, uint32 l){
		(*dart_level_)[d.index] = l;
	}

	inline Dart begin() const
	{
		Dart d(CMAP::darts_->first_index());
		uint32 lastidx = CMAP::darts_->last_index();
		while (dart_level(d) > current_level_ && d.index < lastidx)
			d = Dart(CMAP::darts_->next_index(d.index));
		return d;
	}

	inline Dart end() const
	{
		return Dart(CMAP::darts_->last_index());
	}

	inline Dart next(Dart d) const
	{
		uint32 lastidx = CMAP::darts_->last_index();
		do
		{
			d = Dart(CMAP::darts_->next_index(d.index));
		} while (dart_level(d) > current_level_ && d.index < lastidx);
		return d;
	}
};

} // namespace cgogn
