#include "cph3_adaptative.h"
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cmap/phi.h>

namespace cgogn
{

int CPH3_adaptative::id = 0;

std::shared_ptr<CPH3_adaptative> CPH3_adaptative::get_child()
{
	std::shared_ptr<CPH3_adaptative> result = std::make_shared<CPH3_adaptative>(CPH3_adaptative(CPH3(*this)));
	result->father_ = this;
	return result;
}

uint32 CPH3_adaptative::get_dart_visibility_level(Dart d) const
{
	uint32 result = dart_level(d);
	if (result == 0)
		return result;
	if (representative_is_visible(d))
	{
		auto it = (*dart_visibility_level_)[d.index].begin();
		if (it != (*dart_visibility_level_)[d.index].end() && *it < result)
			result = *it;
	}
	if (father_ != nullptr)
	{
		return std::min(father_->get_dart_visibility_level(d), result);
	}
	return result;
}

void CPH3_adaptative::set_representative(Dart d, Dart r)
{
	(*representative_)[d.index] = r;
}

Dart CPH3_adaptative::get_representative(Dart d) const
{
	if (dart_level(d) == 0)
	{
		return d;
	}
	return (*representative_)[d.index];
}

uint32 CPH3_adaptative::get_representative_visibility_level(Dart d) const
{
	if (dart_level(d) == 0)
		return 0;
	Dart r = get_representative(d);
	auto it = (*representative_visibility_level_)[r.index].begin();
	uint32 result = dart_level(r);
	if (it != (*representative_visibility_level_)[r.index].end() && *it < result)
		result = *it;
	if (father_ != nullptr)
	{
		return std::min(father_->get_representative_visibility_level(d), result);
	}
	return result;
}

void CPH3_adaptative::set_representative_visibility_level(Dart d, uint32 l)
{
	(*representative_visibility_level_)[d.index].insert(l);
}
void CPH3_adaptative::unset_representative_visibility_level(Dart d, uint32 l)
{
	(*representative_visibility_level_)[d.index].erase(l);
}

bool CPH3_adaptative::representative_is_visible(Dart d) const
{
	return get_representative_visibility_level(d) <= current_level_;
}

bool CPH3_adaptative::dart_is_visible(Dart d) const
{
	return get_dart_visibility_level(d) <= current_level_;
}

void CPH3_adaptative::set_visibility_level(Dart d, uint32 l)
{
	(*dart_visibility_level_)[d.index].insert(l);
}

void CPH3_adaptative::unset_visibility_level(Dart d, uint32 l)
{
	(*dart_visibility_level_)[d.index].erase(l);
}

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

uint32 CPH3_adaptative::edge_level(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	return std::max(dart_level(d), dart_level(phi2(*this, d)));
}

Dart CPH3_adaptative::edge_oldest_dart(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	Dart d2 = phi2(*this, d);
	if (dart_level(d) < dart_level(d2))
		return d;
	return d2;
}

Dart CPH3_adaptative::edge_youngest_dart(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	Dart d2 = phi2(*this, d);
	if (dart_level(d) > dart_level(d2))
		return d;
	return d2;
}

bool CPH3_adaptative::edge_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	if (current_level_ == maximum_level_)
		return false;

	Dart d1 = phi1(*this, d);
	CPH3_adaptative m2(*this);
	m2.current_level_ = edge_level(d) + 1;
	Dart d1_l = phi1(m2, d);
	if (d1 != d1_l)
		return true;
	else
		return false;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

uint32 CPH3_adaptative::face_level(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	return dart_level(face_youngest_dart(d));
}

Dart CPH3_adaptative::face_oldest_dart(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	Dart it = d;
	Dart oldest = it;
	uint32 l_old = dart_level(oldest);
	do
	{
		uint32 l = dart_level(it);
		if (l == 0)
			return it;

		if (l < l_old)
		//		if(l < l_old || (l == l_old && it < oldest))
		{
			oldest = it;
			l_old = l;
		}
		it = phi1(*this, it);
	} while (it != d);

	return oldest;
}

Dart CPH3_adaptative::face_youngest_dart(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	Dart it = d;
	Dart youngest = it;
	uint32 l_young;
	Dart it2 = phi_1(*this, it);
	while (get_representative(it) == get_representative(it2))
	{
		it2 = it;
		it = phi1(*this, it);
	}
	l_young = dart_level(it);
	youngest = it;
	it2 = it;
	it = phi1(*this, it);
	while (get_representative(it) == get_representative(it2))
	{
		it2 = it;
		it = phi1(*this, it);
	}
	if (dart_level(it) > l_young)
		youngest = it;

	return youngest;
}

Dart CPH3_adaptative::face_origin(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	Dart p = d;
	uint32 pLevel = dart_level(p);
	CPH3 m2(CPH3(*this));
	do
	{
		m2.current_level_ = pLevel;
		p = m2.face_oldest_dart(p);
		pLevel = dart_level(p);
	} while (pLevel > 0);
	return p;
}

bool CPH3_adaptative::face_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	d = face_youngest_dart(d);
	uint32 fLevel = dart_level(d);
	if (fLevel < current_level_)
		return false;

	bool subd = false;

	CPH3 m2(CPH3(*this));
	m2.current_level_ = fLevel + 1;
	if (dart_level(phi1(m2, d)) == m2.current_level_ && edge_id(phi1(m2, d)) != edge_id(d))
		subd = true;

	return subd;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

uint32 CPH3_adaptative::volume_level(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");

	return dart_level(volume_youngest_dart(d));
}

Dart CPH3_adaptative::volume_oldest_dart(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");

	Dart oldest = d;
	uint32 l_old = dart_level(oldest);
	foreach_incident_face(*this, CPH3_adaptative::CMAP::Volume(oldest), [&](CPH3_adaptative::CMAP::Face f) -> bool {
		Dart old = face_oldest_dart(f.dart);
		uint32 l = dart_level(old);
		if (l < l_old)
		{
			oldest = old;
			l_old = l;
		}
		return true;
	});

	return oldest;
}

Dart CPH3_adaptative::volume_youngest_dart(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");
	// find signifiant edge

	auto is_signifiant_edge = [&](Dart d) {
		Dart tmp = phi2(*this, d);
		return dart_level(get_representative(d)) == 0 || face_id(tmp) != face_id(d);
	};

	Dart it = d;
	Dart it2;
	while (1)
	{
		do
		{
			it2 = it;
			it = phi1(*this, it);
		} while (get_representative(it) == get_representative(it2));
		if (is_signifiant_edge(it))
		{
			break;
		}
		else
		{
			it = phi<21>(*this, it);
		}
	}

	// find significant vertex

	Dart eRep = get_representative(it);
	do
	{
		do
		{
			it2 = it;
			it = phi1(*this, it);
		} while (get_representative(it) == get_representative(it2));
		if (is_signifiant_edge(it))
			break;
		do
		{
			it = phi<21>(*this, it);
		} while (get_representative(it) != eRep);
	} while (1);

	Dart y = it;
	eRep = get_representative(it);
	do
	{
		do
		{
			it2 = it;
			it = phi1(*this, it);
		} while (get_representative(it) == get_representative(it2));
		if (is_signifiant_edge(it))
			break;
		do
		{
			it = phi<21>(*this, it);
		} while (get_representative(it) != eRep);
	} while (1);
	if (dart_level(it) > dart_level(y))
		y = it;
	return y;
}

bool CPH3_adaptative::volume_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d), "Access to a dart not visible at this level");

	uint vLevel = volume_level(d);
	if (vLevel < current_level_)
		return false;

	CPH3 m2(CPH3(*this));
	d = volume_youngest_dart(d);
	m2.current_level_ = vLevel;
	Dart it = phi1(m2, d);
	m2.current_level_ = vLevel + 1;
	Dart it2 = phi1(m2, it);

	if (it == it2)
		return false;

	if (get_representative(it) == get_representative(it2) || face_id(phi2(m2, it2)) == face_id(it2))
		return false;

	return true;
}

/***************************************************
 *            ADAPTATIVE RESOLUTION                *
 ***************************************************/

void CPH3_adaptative::activate_edge_subdivision(CMAP::Edge e)
{
	if (!edge_is_subdivided(e.dart))
		return;
	CMAP::Edge e2 = CMAP::Edge(phi2(*this, e.dart));
	CPH3 m2(CPH3(*this));
	m2.current_level_ = edge_level(e.dart) + 1;
	foreach_dart_of_orbit(m2, e, [&](Dart d) -> bool {
		set_visibility_level(d, current_level_);
		return true;
	});
	foreach_dart_of_orbit(m2, e2, [&](Dart d) -> bool {
		set_visibility_level(d, current_level_);
		return true;
	});
}

void CPH3_adaptative::activate_face_subdivision(CMAP::Face f)
{
	if (!face_is_subdivided(f.dart))
	{
		return;
	}
	Dart d = face_oldest_dart(f.dart);
	CPH3 m2(CPH3(*this));
	m2.current_level_ = face_level(f.dart);
	std::vector<Dart> vect_edges;
	Dart it = d;
	do
	{
		vect_edges.push_back(it);
		it = phi1(m2, it);
	} while (it != d);
	m2.current_level_++;

	for (auto e : vect_edges)
	{
		if (edge_level(e) < m2.current_level_)
			activate_edge_subdivision(CMAP::Edge(e));
		foreach_dart_of_orbit(m2, CMAP::Edge(phi1(m2, e)), [&](Dart dd) -> bool {
			set_visibility_level(dd, current_level_);
			return true;
		});
		Dart it = phi1(m2, e);
		set_representative_visibility_level(it, current_level_);
		set_representative_visibility_level(phi3(m2, it), current_level_);
		it = phi1(m2, it);
		set_representative_visibility_level(it, current_level_);
		set_representative_visibility_level(phi3(m2, it), current_level_);
	}
}

bool CPH3_adaptative::activate_volume_subdivision(CMAP::Volume v)
{
	if (!volume_is_subdivided(v.dart))
	{
		return false;
	}
	Dart d = volume_oldest_dart(v.dart);
	CPH3 m2(CPH3(*this));
	m2.current_level_ = volume_level(v.dart);
	std::vector<Vertex> vect_vertices;
	foreach_incident_vertex(m2, CPH3_adaptative::CMAP::Volume(d), [&](CPH3_adaptative::CMAP::Vertex w) -> bool {
		vect_vertices.push_back(w);
		return true;
	});

	m2.current_level_++;
	for (Vertex w : vect_vertices)
	{
		foreach_dart_of_orbit(m2, Volume(w.dart), [&](Dart d) -> bool {
			foreach_dart_of_orbit(m2, Edge(d), [&](Dart dd) -> bool {
				set_visibility_level(dd, current_level_);
				return true;
			});
			set_representative_visibility_level(d, current_level_);
			set_representative_visibility_level(phi3(m2, d), current_level_);
			return true;
		});
	}
	return true;
}

bool CPH3_adaptative::disable_edge_subdivision(CMAP::Edge e, bool disable_neighbor)
{
	uint32 eLevel = edge_level(e.dart);
	if (eLevel <= current_level_)
		return false;
	Dart d = edge_oldest_dart(e.dart);
	Dart it = d;
	do
	{
		Dart it2 = phi1(*this, it);
		if (get_representative(it) != get_representative(it2))
			return false;
		it = phi<23>(*this, it);
	} while (it != d);
	Dart tmp = phi1(*this, d);
	if (edge_level(tmp) != eLevel)
	{
		if (disable_neighbor)
			disable_edge_subdivision(CMAP::Edge(tmp), true);
		else
			return false;
	}
	CPH3 m2(CPH3(*this));
	m2.current_level_ = eLevel;
	foreach_dart_of_orbit(m2, CMAP::Edge(d), [&](Dart dd) -> bool {
		if (dart_level(dd) == eLevel)
			unset_visibility_level(dd, current_level_);
		return true;
	});
	foreach_dart_of_orbit(m2, CMAP::Edge(tmp), [&](Dart dd) -> bool {
		if (dart_level(dd) == eLevel)
			unset_visibility_level(dd, current_level_);
		return true;
	});
	return true;
}

bool CPH3_adaptative::disable_face_subdivision(CMAP::Face f, bool disable_edge, bool disable_subface)
{
	Dart y = face_youngest_dart(f.dart);
	if (dart_level(y) <= current_level_)
		return false;
	CPH3 m2(CPH3(*this));
	m2.current_level_ = dart_level(y);
	Dart old = m2.face_oldest_dart(y);
	if (face_id(phi2(*this, phi1(m2, old))) != face_id(old))
		return false;
	Dart tmp = phi<13>(m2, old);
	if (face_id(phi2(*this, tmp)) != face_id(old))
		return false;
	m2.current_level_--;
	std::vector<Dart> list_dart;
	Dart it = old;
	do
	{
		list_dart.push_back(it);
		it = phi1(m2, it);
		while (face_level(it) != dart_level(y))
		{
			if (disable_subface)
			{
				disable_face_subdivision(CMAP::Face(it), disable_edge, true);
			}
			else
			{
				return false;
			}
		}
	} while (it != old);
	m2.current_level_++;
	for (Dart d : list_dart)
	{
		Dart dd = phi1(m2, d);
		while (m2.edge_level(dd) != dart_level(y))
			disable_edge_subdivision(CMAP::Edge(dd), true);
		unset_representative_visibility_level(dd, current_level_);
		unset_representative_visibility_level(phi3(m2, dd), current_level_);
		unset_visibility_level(dd, current_level_);
		unset_visibility_level(phi3(m2, dd), current_level_);
		dd = phi1(m2, dd);
		unset_representative_visibility_level(dd, current_level_);
		unset_representative_visibility_level(phi3(m2, dd), current_level_);
		unset_visibility_level(dd, current_level_);
		unset_visibility_level(phi3(m2, dd), current_level_);
		if (disable_edge)
			disable_edge_subdivision(CMAP::Edge(d), true);
	}
	return true;
}

bool CPH3_adaptative::disable_volume_subdivision(CMAP::Volume v, bool disable_face)
{
	Dart y = volume_youngest_dart(v.dart);
	if (dart_level(y) <= current_level_)
		return false;
	CPH3 m2(CPH3(*this));
	m2.current_level_ = dart_level(y) - 1;
	std::vector<Vertex> vect_vertices;
	foreach_incident_vertex(m2, CMAP::Volume(volume_oldest_dart(v.dart)), [&](CMAP::Vertex w) -> bool {
		vect_vertices.push_back(w);
		return true;
	});
	for (auto w : vect_vertices)
	{
		if (volume_level(w.dart) != dart_level(y))
			return false;
	}
	m2.current_level_++;
	bool test = false;
	foreach_incident_face(m2, CMAP::Vertex(phi<12111>(m2, vect_vertices[0].dart)), [&](Face f) -> bool {
		bool tmp = disable_face_subdivision(f, true);
		test = test || tmp;
		return true;
	});
	foreach_dart_of_orbit(m2, CMAP::Vertex(phi<12111>(m2, vect_vertices[0].dart)), [&](Dart d) -> bool {
		Dart d2 = phi2(m2, d);
		Dart d1 = phi1(m2, d);
		Dart d11 = phi1(m2, d1);
		unset_representative_visibility_level(d, current_level_);
		unset_representative_visibility_level(d2, current_level_);
		unset_representative_visibility_level(d1, current_level_);
		unset_representative_visibility_level(d11, current_level_);
		unset_visibility_level(d, current_level_);
		unset_visibility_level(d2, current_level_);
		unset_visibility_level(d1, current_level_);
		unset_visibility_level(d11, current_level_);
		return true;
	});
	if (disable_face)
	{
		m2.current_level_--;
		foreach_incident_face(m2, CPH3_adaptative::CMAP::Volume(volume_oldest_dart(v.dart)), [&](CMAP::Face f) -> bool {
			while (face_level(f.dart) != m2.current_level_ && disable_face_subdivision(f, true, true))
				;
			return true;
		});
	}
	return true;
}

} // namespace cgogn
