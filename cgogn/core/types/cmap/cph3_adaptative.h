#ifndef CGOGN_CORE_TYPES_CMAP_CPH3_ADAPTATIVE_H
#define CGOGN_CORE_TYPES_CMAP_CPH3_ADAPTATIVE_H

#include <cgogn/core/types/cmap/cph3.h>
#include <set>

namespace cgogn
{

struct CPH3_adaptative : public CPH3
{
	using Self = CPH3_adaptative;

	CPH3_adaptative* father_;
	std::shared_ptr<Attribute<std::set<uint32>>> dart_visibility_level_;
	std::shared_ptr<Attribute<std::set<uint32>>> representative_visibility_level_;
	std::shared_ptr<Attribute<Dart>> representative_;
	mutable std::shared_ptr<Attribute<std::pair<uint32, Dart>>> phi1_buffer_;
	mutable std::shared_ptr<Attribute<std::pair<uint32, Dart>>> phi2_buffer_;
	mutable std::shared_ptr<Attribute<std::pair<uint32, Dart>>> phi3_buffer_;
	mutable uint32 clock_;
	mutable uint32 level_clock_;
	static int id;

	CPH3_adaptative() = delete;
	CPH3_adaptative(CMAP& m) : CPH3(m), father_(nullptr), clock_(1), level_clock_(0)
	{
		dart_visibility_level_ =
			m_.darts_->add_attribute<std::set<uint32>>("dart_visibility_level" + std::to_string(id));
		representative_visibility_level_ =
			m_.darts_->add_attribute<std::set<uint32>>("representative_visibility_level" + std::to_string(id));
		representative_ = m_.darts_->get_attribute<Dart>("representative");
		if (!representative_)
			representative_ = m_.darts_->add_attribute<Dart>("representative");
		phi1_buffer_ = m_.darts_->get_attribute<std::pair<uint32, Dart>>("phi1_buffer" + std::to_string(id));
		if (!phi1_buffer_)
			phi1_buffer_ = m_.darts_->add_attribute<std::pair<uint32, Dart>>("phi1_buffer" + std::to_string(id));
		phi2_buffer_ = m_.darts_->get_attribute<std::pair<uint32, Dart>>("phi2_buffer" + std::to_string(id));
		if (!phi2_buffer_)
			phi2_buffer_ = m_.darts_->add_attribute<std::pair<uint32, Dart>>("phi2_buffer" + std::to_string(id));
		phi3_buffer_ = m_.darts_->get_attribute<std::pair<uint32, Dart>>("phi3_buffer" + std::to_string(id));
		if (!phi3_buffer_)
			phi3_buffer_ = m_.darts_->add_attribute<std::pair<uint32, Dart>>("phi3_buffer" + std::to_string(id));
		id++;
	}
	CPH3_adaptative(const CPH3& cph) : CPH3(cph), father_(nullptr), clock_(1), level_clock_(0)
	{
		dart_visibility_level_ =
			m_.darts_->add_attribute<std::set<uint32>>("dart_visibility_level" + std::to_string(id));
		representative_visibility_level_ =
			m_.darts_->add_attribute<std::set<uint32>>("representative_visibility_level" + std::to_string(id));
		representative_ = m_.darts_->get_attribute<Dart>("representative");
		if (!representative_)
			representative_ = m_.darts_->add_attribute<Dart>("representative");
		phi1_buffer_ = m_.darts_->get_attribute<std::pair<uint32, Dart>>("phi1_buffer" + std::to_string(id));
		if (!phi1_buffer_)
			phi1_buffer_ = m_.darts_->add_attribute<std::pair<uint32, Dart>>("phi1_buffer" + std::to_string(id));
		phi2_buffer_ = m_.darts_->get_attribute<std::pair<uint32, Dart>>("phi2_buffer" + std::to_string(id));
		if (!phi2_buffer_)
			phi2_buffer_ = m_.darts_->add_attribute<std::pair<uint32, Dart>>("phi2_buffer" + std::to_string(id));
		phi3_buffer_ = m_.darts_->get_attribute<std::pair<uint32, Dart>>("phi3_buffer" + std::to_string(id));
		if (!phi3_buffer_)
			phi3_buffer_ = m_.darts_->add_attribute<std::pair<uint32, Dart>>("phi3_buffer" + std::to_string(id));
		id++;
	}
	CPH3_adaptative(const CPH3_adaptative& other)
		: CPH3(other), father_(other.father_), dart_visibility_level_(other.dart_visibility_level_),
		  representative_visibility_level_(other.representative_visibility_level_),
		  representative_(other.representative_), phi1_buffer_(nullptr), phi2_buffer_(nullptr), phi3_buffer_(nullptr),
		  clock_(1), level_clock_(0)
	{
	}

	inline Dart begin() const
	{
		Dart d(m_.darts_->first_index());
		uint32 lastidx = m_.darts_->last_index();
		while (!dart_is_visible(d) && d.index < lastidx)
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
		} while (!dart_is_visible(d) && d.index < lastidx);
		return d;
	}

	bool get_phi1_buffer(Dart d, Dart& result) const;
	void set_phi1_buffer(Dart d, Dart d2) const;
	bool get_phi2_buffer(Dart d, Dart& result) const;
	void set_phi2_buffer(Dart d, Dart d2) const;
	bool get_phi3_buffer(Dart d, Dart& result) const;
	void set_phi3_buffer(Dart d, Dart d2) const;
	CPH3_adaptative* get_copy();

	CPH3_adaptative* get_child();

	void set_representative(Dart d, Dart r);
	Dart get_representative(Dart d) const;
	uint32 get_dart_visibility_level(Dart d) const;
	uint32 get_representative_visibility_level(Dart d) const;
	bool representative_is_visible(Dart d) const;
	bool dart_is_visible(Dart d) const;

	void set_representative_visibility_level(Dart d, uint32 l);
	void unset_representative_visibility_level(Dart d, uint32 l);
	void set_visibility_level(Dart d, uint32 l);
	void unset_visibility_level(Dart d, uint32 l);

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	uint32 edge_level(Dart d) const;
	Dart edge_oldest_dart(Dart d) const;
	Dart edge_youngest_dart(Dart d) const;
	bool edge_is_subdivided(Dart d) const;

	/***************************************************
	 *                  FACE INFO                      *
	 ***************************************************/

	uint32 face_level(Dart d) const;
	Dart face_oldest_dart(Dart d) const;
	Dart face_youngest_dart(Dart d) const;
	Dart face_origin(Dart d) const;
	bool face_is_subdivided(Dart d) const;

	/***************************************************
	 *                 VOLUME INFO                     *
	 ***************************************************/

	uint32 volume_level(Dart d) const;
	Dart volume_oldest_dart(Dart d) const;
	Dart volume_youngest_dart(Dart d) const;
	Dart volume_origin(Dart d) const;
	bool volume_is_subdivided(Dart d) const;
	std::vector<Dart> get_dart_volume_level(Dart d) const;
	std::vector<Dart> get_volume_hierarchical_darts(Dart d) const;
	void change_volume_level(Dart d, int32 l);
	bool raise_volume_level(CMAP::Volume v);

	/***************************************************
	 *            ADAPTATIVE RESOLUTION                *
	 ***************************************************/

	void activate_edge_subdivision(CMAP::Edge e);
	void activate_face_subdivision(CMAP::Face f);
	bool activate_volume_subdivision(CMAP::Volume v);

	bool disable_edge_subdivision(CMAP::Edge e, bool disable_neighbor = false);
	bool disable_face_subdivision(CMAP::Face f, bool disable_edge = false, bool disable_subface = false);
	bool disable_volume_subdivision(CMAP::Volume v, bool disable_face = false);
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CPH3_ADAPTATIVE_H
