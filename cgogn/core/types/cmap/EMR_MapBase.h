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
struct CGOGN_CORE_EXPORT EMR_MapBase_T : public CMAP
{

	using MAP = CMAP;
	template <typename T>
	using Attribute = typename CMAP::template Attribute<T>;
	using AttributeGen = typename MAP::AttributeGen;
	using MarkAttribute = typename MAP::MarkAttribute;

	std::shared_ptr<Attribute<uint32>> dart_level_;
	std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>> MR_relation_;

	EMR_MapBase_T() : CMAP()
	{
		MR_relation_ = std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>>(
			new std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>());
		dart_level_ = CMAP::darts_->template get_attribute<uint32>("dart_level");
		if (!dart_level_)
			dart_level_ = CMAP::darts_->template add_attribute<uint32>("dart_level");
	}
	virtual ~EMR_MapBase_T()
	{
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
			auto new_rel = CMAP::add_relation((*r)[0]->name() + "_" + std::to_string(max));
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
};

template <typename EMR>
struct EMR_MapBase
{

	using MAP = EMR;
	using MarkAttribute = typename MAP::MarkAttribute;
	EMR& m_;
	uint32 current_level_;

	EMR_MapBase(EMR& m) : m_(m), current_level_(0)
	{
	}

	operator EMR&()
	{
		return m_;
	}
	operator const EMR&() const
	{
		return m_;
	}

	void change_resolution_level(uint32 new_level)
	{
		cgogn_message_assert(0 <= new_level && new_level < m_.max_level(), "Access to an undefined level");
		current_level_ = new_level;
	}

	inline Dart begin() const
	{
		Dart d(m_.darts_->first_index());
		uint32 lastidx = m_.darts_->last_index();
		while (m_.dart_level(d) > current_level_ && d.index < lastidx)
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
		} while (m_.dart_level(d) > current_level_ && d.index < lastidx);
		return d;
	}
};

} // namespace cgogn

#endif
