#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAP3_ADAPTATIVE_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAP3_ADAPTATIVE_H_

#include <cgogn/core/types/cmap/EMR_Map3.h>

namespace cgogn
{

struct EMR_Map3_Adaptative : EMR_Map3
{

	using Inherit = EMR_Map3;
	using MAP = CMap3;
	using CMAP = CMap3;
	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Cell<PHI2>;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	template <typename T>
	using Attribute = typename CMap3::template Attribute<T>;

	std::shared_ptr<Attribute<std::pair<bool, uint32>>> dart_visibility_;
	mutable std::shared_ptr<Attribute<std::pair<bool, Dart>>> dart_representative_;
	mutable std::shared_ptr<Attribute<std::tuple<uint32, uint32, uint32, Dart>>> phi3_buffer_;
	mutable std::shared_ptr<Attribute<std::tuple<uint32, uint32, uint32, Dart>>> phi2_buffer_;
	mutable std::shared_ptr<Attribute<std::tuple<uint32, uint32, uint32, Dart>>> phi1_buffer_;
	mutable std::shared_ptr<Attribute<std::tuple<uint32, uint32, uint32, Dart>>> volume_dart_buffer_;
	mutable std::shared_ptr<Attribute<std::tuple<uint32, uint32, uint32, uint32>>> dart_visibility_buffer_;

	static uint32 nb_views;
	EMR_Map3_Adaptative* parent;

	mutable uint32 clock_views_;
	mutable uint32 clock_parent_;
	EMR_Map3_Adaptative*& topology_;
	std::vector<EMR_Map3_Adaptative*>& list_view_;

	EMR_Map3_Adaptative(EMR_Map3_T<CMap3>& m)
		: EMR_Map3(m), parent(nullptr), clock_views_(0), clock_parent_(0),
		  topology_(m.template get_attribute<EMR_Map3_Adaptative*>("emr_topology")),
		  list_view_(m.template get_attribute<std::vector<EMR_Map3_Adaptative*>>("list_view"))
	{
		nb_views++;
		dart_visibility_ =
			m_.darts_->get_attribute<std::pair<bool, uint32>>("dart_visibility" + std::to_string(nb_views));
		if (!dart_visibility_)
		{
			dart_visibility_ =
				m_.darts_->add_attribute<std::pair<bool, uint32>>("dart_visibility" + std::to_string(nb_views));
		}
		dart_representative_ = m_.darts_->get_attribute<std::pair<bool, Dart>>("dart_representative");
		if (!dart_representative_)
		{
			dart_representative_ = m_.darts_->add_attribute<std::pair<bool, Dart>>("dart_representative");
		}
		phi1_buffer_ = m_.darts_->get_attribute<std::tuple<uint32, uint32, uint32, Dart>>("phi1_buffer" +
																						  std::to_string(nb_views));
		if (!phi1_buffer_)
		{
			phi1_buffer_ = m_.darts_->add_attribute<std::tuple<uint32, uint32, uint32, Dart>>("phi1_buffer" +
																							  std::to_string(nb_views));
		}
		phi2_buffer_ = m_.darts_->get_attribute<std::tuple<uint32, uint32, uint32, Dart>>("phi2_buffer" +
																						  std::to_string(nb_views));
		if (!phi2_buffer_)
		{
			phi2_buffer_ = m_.darts_->add_attribute<std::tuple<uint32, uint32, uint32, Dart>>("phi2_buffer" +
																							  std::to_string(nb_views));
		}
		phi3_buffer_ = m_.darts_->get_attribute<std::tuple<uint32, uint32, uint32, Dart>>("phi3_buffer" +
																						  std::to_string(nb_views));
		if (!phi3_buffer_)
		{
			phi3_buffer_ = m_.darts_->add_attribute<std::tuple<uint32, uint32, uint32, Dart>>("phi3_buffer" +
																							  std::to_string(nb_views));
		}

		volume_dart_buffer_ = m_.darts_->get_attribute<std::tuple<uint32, uint32, uint32, Dart>>(
			"volume_dart_buffer" + std::to_string(nb_views));
		if (!volume_dart_buffer_)
		{
			volume_dart_buffer_ = m_.darts_->add_attribute<std::tuple<uint32, uint32, uint32, Dart>>(
				"volume_dart_buffer" + std::to_string(nb_views));
		}
		dart_visibility_buffer_ = m_.darts_->get_attribute<std::tuple<uint32, uint32, uint32, uint32>>(
			"dart_visibility_buffer" + std::to_string(nb_views));
		if (!dart_visibility_buffer_)
		{
			dart_visibility_buffer_ = m_.darts_->add_attribute<std::tuple<uint32, uint32, uint32, uint32>>(
				"dart_visibility_buffer" + std::to_string(nb_views));
		}
		list_view_.push_back(this);
	}

	virtual ~EMR_Map3_Adaptative()
	{
		m_.darts_->remove_attribute(dart_visibility_);
		m_.darts_->remove_attribute(phi1_buffer_);
		m_.darts_->remove_attribute(phi2_buffer_);
		m_.darts_->remove_attribute(phi3_buffer_);
		m_.darts_->remove_attribute(volume_dart_buffer_);
		for (uint32 i = 0; i < list_view_.size(); i++)
		{
			if (list_view_[i] == this)
			{
				std::swap(list_view_[i], list_view_[list_view_.size() - 1]);
				list_view_.pop_back();
				break;
			}
		}
	}

	virtual bool check_integrity() const;

	uint32 get_dart_visibility(Dart d) const;
	uint32 get_dart_visibility_fast(Dart d) const;
	void set_dart_visibility(Dart d, uint32 v);
	bool dart_is_visible(Dart d) const;

	Dart get_representative(Dart d) const;
	EMR_Map3_Adaptative* get_parent() const;

	Dart get_phi1_buffer(Dart d) const;
	Dart get_phi2_buffer(Dart d) const;
	Dart get_phi3_buffer(Dart d) const;

	Dart begin() const;

	Dart end() const;

	Dart next(Dart d) const;

	EMR_Map3_Adaptative* get_child();
	EMR_Map3_Adaptative* get_copy();
	void copy_visibility(const EMR_Map3_Adaptative& other);

	bool vertex_is_visible(Dart d) const;

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	Dart edge_youngest_dart(Dart d) const;
	Dart edge_oldest_dart(Dart d) const;
	bool edge_is_subdivided(Dart d) const;
	uint32 edge_level(Dart d) const;
	bool is_topologycal_edge(Dart d) const;

	/***************************************************
	 *                  FACE INFO                      *
	 ***************************************************/

	Dart face_youngest_dart(Dart d) const;
	Dart face_oldest_dart(Dart d) const;
	bool face_is_subdivided(Dart d) const;
	uint32 face_level(Dart d) const;
	bool is_topologycal_face(Dart d) const;

	/***************************************************
	 *                 VOLUME INFO                     *
	 ***************************************************/
	Dart volume_youngest_dart(Dart d) const;
	Dart volume_oldest_dart(Dart d) const;
	bool volume_is_subdivided(Dart d) const;
	uint32 volume_level(Dart d) const;
	bool is_topologycal_volume(Dart d) const;

	/***************************************************
	 *            ADAPTATIVE SUBDIVISION               *
	 ***************************************************/

	void activate_edge_subdivision(Edge e);
	void activate_face_subdivision(Face f);
	bool activate_volume_subdivision(Volume v);
	bool activate_volume_subdivision_fast(Volume v);

	bool disable_edge_subdivision(Edge e);
	bool disable_face_subdivision(Face f, bool disable_edge = false, bool disable_subface = false);
	bool disable_volume_subdivision(Volume v, bool disable_face = false);
	bool disable_volume_subdivision_fast(Volume v, bool disable_face = false);
};

template <>
struct mesh_traits<EMR_Map3_Adaptative> : public mesh_traits<CMap3>
{
	static constexpr const char* name = "EMR_Map3_Adaptative";
};

} // namespace cgogn

#endif
