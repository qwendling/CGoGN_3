#include "EMR_Map3_Adaptative.h"
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

bool EMR_Map3_Adaptative::check_integrity() const
{
	for (Dart d = this->begin(), end = this->end(); d != end; d = this->next(d))
	{
		int limit = INT_MAX;
		int i = 0;
		Dart it = phi1(*this, d);
		bool is_permutation = false;
		while (i < limit && !is_permutation)
		{
			is_permutation = it == d;
			it = phi1(*this, it);
			i++;
		}
		if (i == 0)
		{
			std::cerr << "phi1 have a fix point : " << d.index << std::endl;
			return false;
		}
		if (!is_permutation)
		{
			std::cerr << "phi1 must be a permutation" << std::endl;
			return false;
		}
		if (phi2(*this, phi2(*this, d)) != d)
		{
			std::cerr << "phi2 must be an involution" << std::endl;
			return false;
		}
		if (phi3(*this, phi3(*this, d)) != d)
		{
			std::cerr << "phi3 must be an involution" << std::endl;
			return false;
		}
		if (phi1(*this, phi3(*this, phi1(*this, phi3(*this, d)))) != d)
		{
			std::cerr << "phi1(phi3(d)) must be an involution" << std::endl;
			return false;
		}
	}
	return true;
}

uint32 EMR_Map3_Adaptative::get_dart_visibility(Dart d) const
{
	return std::min(this->dart_level(d), (*dart_visibility_)[d.index]);
}

void EMR_Map3_Adaptative::set_dart_visibility(Dart d, uint32 v)
{
	(*dart_visibility_)[d.index] = v;
}

Dart EMR_Map3_Adaptative::begin() const
{
	Dart d(m_.darts_->first_index());
	uint32 lastidx = m_.darts_->last_index();
	while (d.index < lastidx && get_dart_visibility(d) > current_level_)
		d = Dart(m_.darts_->next_index(d.index));
	return d;
}

Dart EMR_Map3_Adaptative::end() const
{
	return Dart(m_.darts_->last_index());
}

Dart EMR_Map3_Adaptative::next(Dart d) const
{
	uint32 lastidx = m_.darts_->last_index();
	do
	{
		d = Dart(m_.darts_->next_index(d.index));
	} while (d.index < lastidx && get_dart_visibility(d) > current_level_);
	return d;
}

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

Dart EMR_Map3_Adaptative::edge_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = phi2(*this, d);
	if (m_.dart_level(d) > m_.dart_level(it))
		return d;
	return it;
}

Dart EMR_Map3_Adaptative::edge_oldest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = phi2(*this, d);
	if (m_.dart_level(d) < m_.dart_level(it))
		return d;
	return it;
}

uint32 EMR_Map3_Adaptative::edge_level(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	return m_.dart_level(edge_youngest_dart(d));
}

bool EMR_Map3_Adaptative::edge_is_subdivided(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	uint32 e_level = edge_level(d);
	if (e_level == maximum_level_)
		return false;
	EMR_Map3_Adaptative m2(m_);
	m2.current_level_ = e_level + 1;
	if (phi1(*this, d) == phi1(m2, d))
	{
		return false;
	}
	return true;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

Dart EMR_Map3_Adaptative::face_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (edge_level(d) == 0)
	{
		return d;
	}
	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	Dart it, it2;
	it = d;
	it2 = phi_1(*this, d);
	do
	{
		if (dart_level(it) == dart_level(it2))
			return it;
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		it2 = it;
		it = phi1(*this, it);
	} while (it != d);
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = true;
		Dart result_young = phi1(m2, old);
		result = marker.is_marked(result_young);
		if (result)
		{
			return result_young;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo face level");
	} while (!result);
	return d;
}

Dart EMR_Map3_Adaptative::face_oldest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = phi1(*this, d);
	Dart result = d;
	while (it != d)
	{
		if (m_.dart_level(it) < m_.dart_level(result))
			result = it;
		it = phi1(*this, it);
	}

	return result;
}

bool EMR_Map3_Adaptative::face_is_subdivided(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = face_level(d);
	return m2.face_is_subdivided(d);
}

uint32 EMR_Map3_Adaptative::face_level(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (edge_level(d) == 0)
	{
		return 0;
	}
	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	foreach_dart_of_orbit(*this, Face2(d), [&](Dart it) -> bool {
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		return true;
	});
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = marker.is_marked(phi1(m2, old));
		if (result)
		{
			return m2.current_level_;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo face level");
	} while (!result);
	return m2.current_level_;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

Dart EMR_Map3_Adaptative::volume_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");

	if (edge_level(d) == 0)
	{
		return d;
	}
	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		return true;
	});
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = true;
		Dart result_young = phi1(m2, old);
		result = marker.is_marked(result_young);
		if (result)
		{
			return result_young;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo volume youngest");
	} while (!result);
	return d;
}

Dart EMR_Map3_Adaptative::volume_oldest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart result = d;
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		if (m_.dart_level(it) < m_.dart_level(result))
		{
			result = it;
		}
		return true;
	});
	return result;
}
bool EMR_Map3_Adaptative::volume_is_subdivided(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = volume_level(d);
	return m2.volume_is_subdivided(d);
}

uint32 EMR_Map3_Adaptative::volume_level(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (edge_level(d) == 0)
	{
		return 0;
	}

	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		return true;
	});
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = marker.is_marked(phi1(m2, old));
		if (result)
		{
			return m2.current_level_;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo volume level");
	} while (!result);
	return m2.current_level_;
}

/***************************************************
 *            ADAPTATIVE SUBDIVISION               *
 ***************************************************/

void EMR_Map3_Adaptative::activate_edge_subdivision(Edge e)
{
	if (!edge_is_subdivided(e.dart))
		return;

	Edge e2 = Edge(phi2(*this, e.dart));
	EMR_Map3 m2(m_);
	m2.current_level_ = edge_level(e.dart) + 1;
	foreach_dart_of_orbit(m2, e, [&](Dart d) -> bool {
		if (get_dart_visibility(phi3(m2, d)) <= current_level_)
			set_dart_visibility(d, current_level_);
		return true;
	});
	foreach_dart_of_orbit(m2, e2, [&](Dart d) -> bool {
		if (get_dart_visibility(phi3(m2, d)) <= current_level_)
			set_dart_visibility(d, current_level_);
		return true;
	});
	return;
}
void EMR_Map3_Adaptative::activate_face_subdivision(Face f)
{
	if (!face_is_subdivided(f.dart))
	{
		return;
	}

	Dart d = face_oldest_dart(f.dart);
	EMR_Map3 m2(m_);
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
			activate_edge_subdivision(Edge(e));
		Dart it2 = phi1(m2, e);
		while (get_dart_visibility(it2) > current_level_)
		{
			set_dart_visibility(it2, current_level_);
			set_dart_visibility(phi3(m2, it2), current_level_);
			it2 = phi1(m2, it2);
		}
	}
}
bool EMR_Map3_Adaptative::activate_volume_subdivision(Volume v)
{
	if (!volume_is_subdivided(v.dart))
	{
		return false;
	}
	Dart d = volume_oldest_dart(v.dart);
	uint32 v_level = volume_level(v.dart);
	EMR_Map3 m2(m_);
	m2.current_level_ = v_level;
	std::vector<Vertex> vect_vertices;
	foreach_incident_vertex(m2, Volume(d), [&](Vertex w) -> bool {
		vect_vertices.push_back(w);
		return true;
	});

	m2.current_level_++;
	uint32 tmp = m2.current_level_;
	std::vector<Dart> intern_darts;
	for (Vertex w : vect_vertices)
	{
		foreach_dart_of_orbit(m2, Volume(w.dart), [&](Dart d) -> bool {
			if (get_dart_visibility(d) > current_level_)
			{
				Dart d2 = phi2(m2, d);
				if (get_dart_visibility(d2) <= current_level_)
				{
					uint32 e_level = edge_level(d2);
					uint32 tmp = m2.current_level_;
					std::vector<Dart> vec_dart;
					vec_dart.push_back(d2);
					while (m2.current_level_ <= e_level)
					{
						uint32 size_vec = vec_dart.size();
						for (uint32 i = 0; i < size_vec; i++)
						{
							vec_dart.push_back(phi2(m2, vec_dart[i]));
							intern_darts.push_back(phi2(m2, vec_dart[i]));
						}
						m2.current_level_++;
					}
					m2.current_level_ = tmp;
				}
			}
			return true;
		});
	}
	for (auto id : intern_darts)
	{
		set_dart_visibility(id, current_level_);
	}
	m2.current_level_ = tmp;

	for (Vertex w : vect_vertices)
	{
		foreach_dart_of_orbit(m2, Volume(w.dart), [&](Dart d) -> bool {
			if (get_dart_visibility(d) <= current_level_)
			{
				foreach_dart_of_orbit(m2, Edge(d), [&](Dart dd) -> bool {
					if (get_dart_visibility(dd) <= current_level_)
					{
						set_dart_visibility(phi3(m2, dd), current_level_);
					}
					return true;
				});
			}
			else
			{
				set_dart_visibility(d, current_level_);
				set_dart_visibility(phi3(m2, d), current_level_);
			}
			return true;
		});
	}
	return true;
}

bool EMR_Map3_Adaptative::disable_edge_subdivision(Edge e)
{
	uint32 e_level = edge_level(e.dart);
	if (e_level <= current_level_)
		return false;
	Dart old = edge_oldest_dart(e.dart);
	if (dart_level(old) == e_level)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = e_level - 1;
	Dart d2 = phi<23>(*this, old);
	d2 = phi3(m2, d2);
	while (edge_level(d2) != e_level)
		if (!disable_edge_subdivision(Edge(d2)))
			return false;
	Dart test = old;
	do
	{
		if (phi2(*this, phi1(*this, test)) != phi2(m2, test))
			return false;
		test = phi2(m2, phi3(m2, test));
	} while (test != old);
	m2.current_level_ = e_level;
	Dart it, it2;
	it = old;
	it2 = d2;
	do
	{
		Dart tmp = phi2(m2, it);
		Dart tmp2 = phi2(m2, it2);
		set_dart_visibility(tmp, UINT_MAX);
		set_dart_visibility(tmp2, UINT_MAX);
		it = phi3(m2, tmp);
		it2 = phi3(m2, tmp2);
	} while (it != old);

	return true;
}
bool EMR_Map3_Adaptative::disable_face_subdivision(Face f, bool disable_edge, bool disable_subface)
{
	uint32 f_level = face_level(f.dart);
	if (f_level == 0)
		return false;
	uint32 v_level = volume_level(f.dart);
	if (v_level == f_level)
		return false;
	uint32 v_level3 = volume_level(phi3(*this, f.dart));
	if (v_level3 == f_level)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = f_level - 1;
	std::vector<Dart> vec_vertices;
	Dart old = face_oldest_dart(f.dart);
	Dart it = old;
	do
	{
		vec_vertices.push_back(it);
		if (face_level(it) != f_level)
		{
			if (!disable_subface)
			{
				return false;
			}
			if (!disable_face_subdivision(Face(it), disable_edge, disable_subface))
				return false;
		}
		it = phi1(m2, it);
	} while (it != old);

	m2.current_level_ = f_level;
	for (Dart d : vec_vertices)
	{
		Dart tmp = phi1(m2, d);
		while (edge_level(tmp) != f_level)
			disable_edge_subdivision(Edge(tmp));
	}
	return true;
}
bool EMR_Map3_Adaptative::disable_volume_subdivision(Volume v, bool disable_face)
{
	return false;
}

} // namespace cgogn
