#ifndef CGOGN_CORE_TYPES_CMAP_EMR_COMPACT_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_COMPACT_H_

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
struct CGOGN_CORE_EXPORT EMR_Map3_Compact : public CMap3
{

	using MAP = CMap3;
	template <typename T>
	using Attribute = typename MAP::template Attribute<T>;
	using AttributeGen = typename MAP::AttributeGen;
	using MarkAttribute = typename MAP::MarkAttribute;

	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Cell<PHI2>;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;

	std::shared_ptr<Attribute<uint32>> dart_level_;
	std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>> MR_relation_;
	uint32& maximum_level_;
	uint32 current_level_;
	uint32& clock_;
	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi1_;
	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi_1_;
	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi2_;
	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> MR_phi3_;

	EMR_Map3_Compact()
		: MAP(), maximum_level_(MAP::template get_attribute<uint32>("emr_maximum_level")), current_level_(0),
		  clock_(MAP::template get_attribute<uint32>("emr_clock"))
	{
		MR_relation_ = std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>>(
			new std::vector<std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>>());
		dart_level_ = MAP::darts_->template get_attribute<uint32>("dart_level");
		if (!dart_level_)
			dart_level_ = MAP::darts_->template add_attribute<uint32>("dart_level");
		if (clock_ == 0)
			clock_ = 1;
		MR_phi1_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi_1_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi1_->push_back(this->phi1_);
		MR_phi_1_->push_back(this->phi_1_);
		MR_phi2_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi2_->push_back(this->phi2_);
		MR_phi3_ = this->MR_relation_->emplace_back(new std::vector<std::shared_ptr<Attribute<Dart>>>());
		MR_phi3_->push_back(this->phi3_);
	}

	bool check_integrity() const;

	void add_resolution()
	{
		uint32 max = maximum_level_;
		for (auto& r : *MR_relation_)
		{
			auto new_rel = MAP::add_relation((*r)[0]->name() + "_" + std::to_string(max));
			r->push_back(new_rel);
			for (Dart d = this->begin(), end = this->end(); d != end; d = this->next(d))
			{
				(*new_rel)[d.index] = (*(*r)[max])[d.index];
			}
			// new_rel->copy((*r)[max].get());
		}
		maximum_level_++;
		this->phi1_ = (*MR_phi1_)[this->maximum_level_];
		this->phi_1_ = (*MR_phi_1_)[this->maximum_level_];
		this->phi2_ = (*MR_phi2_)[this->maximum_level_];
		this->phi3_ = (*MR_phi3_)[this->maximum_level_];
	}

	uint32 dart_level(Dart d) const
	{
		return (*dart_level_)[d.index];
	}
	void set_dart_level(Dart d, uint32 l)
	{
		(*dart_level_)[d.index] = l;
	}

	void change_resolution_level(uint32 new_level)
	{
		cgogn_message_assert(0 <= new_level && new_level <= maximum_level_, "Access to an undefined level");
		current_level_ = new_level;
	}

	inline Dart begin() const
	{
		Dart d(darts_->first_index());
		uint32 lastidx = darts_->last_index();
		while (d.index < lastidx && dart_level(d) > current_level_)
			d = Dart(darts_->next_index(d.index));
		return d;
	}

	inline Dart end() const
	{
		return Dart(darts_->last_index());
	}

	inline Dart next(Dart d) const
	{
		uint32 lastidx = darts_->last_index();
		do
		{
			d = Dart(darts_->next_index(d.index));
		} while (d.index < lastidx && dart_level(d) > current_level_);
		return d;
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

template <>
struct mesh_traits<EMR_Map3_Compact> : public mesh_traits<CMap3>
{
	static constexpr const char* name = "EMR_Map3_Compact";
};

} // namespace cgogn

#endif
