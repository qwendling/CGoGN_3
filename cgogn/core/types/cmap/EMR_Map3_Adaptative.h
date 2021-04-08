#ifndef CGOGN_CORE_TYPES_CMAP_EMR_MAP3_ADAPTATIVE_H_
#define CGOGN_CORE_TYPES_CMAP_EMR_MAP3_ADAPTATIVE_H_

#include <cgogn/core/types/cmap/EMR_Map3.h>

namespace cgogn
{

struct EMR_Map3_Adaptative : EMR_Map3
{

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
	static uint32 nb_views;
	EMR_Map3_Adaptative* parent;

	EMR_Map3_Adaptative(EMR_Map3_T<CMap3>& m) : EMR_Map3(m), parent(nullptr)
	{
		nb_views++;
		dart_visibility_ =
			m_.darts_->get_attribute<std::pair<bool, uint32>>("dart_visibility" + std::to_string(nb_views));
		if (!dart_visibility_)
		{
			dart_visibility_ =
				m_.darts_->add_attribute<std::pair<bool, uint32>>("dart_visibility" + std::to_string(nb_views));
		}
	}

	virtual bool check_integrity() const;

	uint32 get_dart_visibility(Dart d) const;
	void set_dart_visibility(Dart d, uint32 v);
	bool dart_is_visible(Dart d) const;

	Dart begin() const;

	Dart end() const;

	Dart next(Dart d) const;

	EMR_Map3_Adaptative* get_child();
	EMR_Map3_Adaptative* get_copy();

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	Dart edge_youngest_dart(Dart d) const;
	Dart edge_oldest_dart(Dart d) const;
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

	/***************************************************
	 *            ADAPTATIVE SUBDIVISION               *
	 ***************************************************/

	void activate_edge_subdivision(Edge e);
	void activate_face_subdivision(Face f);
	bool activate_volume_subdivision(Volume v);

	bool disable_edge_subdivision(Edge e);
	bool disable_face_subdivision(Face f, bool disable_edge = false, bool disable_subface = false);
	bool disable_volume_subdivision(Volume v, bool disable_face = false);
};

} // namespace cgogn

#endif
