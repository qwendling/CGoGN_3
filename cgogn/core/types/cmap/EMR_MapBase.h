#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAPBASE_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAPBASE_H_

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

template <typename CMAP>
struct CGOGN_CORE_EXPORT EMR_MapBase_T
{

	using MAP = CMAP;
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;
	using AttributeGen = typename MAP::AttributeGen;
	using MarkAttribute = typename MAP::MarkAttribute;

	std::shared_ptr<Attribute<uint32>> dart_level_;
	std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>> MR_relation_;

	uint32 current_level_;

	CMAP& m_;

	EMR_MapBase_T(CMAP& m) : m_(m), current_level_(0)
	{
		MR_relation_ = std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>>(
			new std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>());
		dart_level_ = m_.darts_->template get_attribute<uint32>("dart_level");
		if (!dart_level_)
			dart_level_ = m_.darts_->template add_attribute<uint32>("dart_level");
	}
	virtual ~EMR_MapBase_T()
	{
	}

	operator MAP&()
	{
		return m_;
	}
	operator const MAP&() const
	{
		return m_;
	}

	void change_resolution_level(uint32 new_level)
	{
		cgogn_message_assert(0 <= new_level && new_level < this->max_level(), "Access to an undefined level");
		current_level_ = new_level;
	}

	uint32 max_level() const
	{
		return (*MR_relation_)[0]->size();
	}

	void add_resolution()
	{
		uint32 max = max_level();
		for (auto& r : *MR_relation_)
		{
			auto new_rel = m_.add_relation((*r)[0]->name() + "_" + std::to_string(max));
			r->push_back(new_rel);
			new_rel->copy((*r)[max - 1].get());
		}
	}

	uint32 dart_level(Dart d) const
	{
		return (*dart_level_)[d.index];
	}
	void set_dart_level(Dart d, uint32 l)
	{
		(*dart_level_)[d.index] = l;
	}

	inline Dart begin() const
	{
		Dart d(m_.darts_->first_index());
		uint32 lastidx = m_.darts_->last_index();
		while (dart_level(d) > current_level_ && d.index < lastidx)
			d = Dart(m_.darts_->next_index(d.index));
		return d;
	}

	inline Dart end() const
	{
		return Dart(m_.darts_->last_index());
	}

	inline Dart next(Dart d) const
	{
		uint32 lastidx = m_.darts_->last_index();
		do
		{
			d = Dart(m_.darts_->next_index(d.index));
		} while (dart_level(d) > current_level_ && d.index < lastidx);
		return d;
	}
};

} // namespace cgogn

#endif
