#include "EMR_Map3_Adaptative.h"
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{
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
	uint32 f_level = face_level(d);
	Dart it = d;
	while (dart_level(it) != f_level)
	{
		it = phi1(*this, it);
	}

	return it;
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
		return 0;
	uint32 min_e_lvl = INT_MAX;
	Dart it = phi1(*this, d);
	Dart it2 = it;
	while (it != d)
	{
		uint32 tmp = edge_level(it);
		if (tmp < min_e_lvl)
		{
			min_e_lvl = tmp;
		}
		if (dart_level(it) < dart_level(it2))
			it2 = it;
		it = phi1(*this, it);
	}
	EMR_Map3 m2(m_);
	m2.current_level_ = min_e_lvl + 1;
	Dart d1, d2;
	do
	{
		m2.current_level_--;
		d1 = it2;
		d2 = it2;
		do
		{
			if (d1 == d2)
				d2 = phi1(m2, d2);
			d1 = phi1(*this, d1);
		} while (d1 != it2 && d2 != it2);
		if (d2 == it2 && m2.current_level_ > 0)
		{
			return m2.face_level(d2);
		}
	} while (d2 != it2);

	return m2.current_level_;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

Dart EMR_Map3_Adaptative::volume_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	uint32 v_level = volume_level(d);
	Dart result;
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		if (dart_level(it) == v_level)
		{
			result = it;
			return false;
		}
		return true;
	});

	return result;
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
		return 0;
	Dart old = volume_oldest_dart(d);
	EMR_Map3 m2(m_);
	m2.current_level_ = std::max(current_level_, dart_level(old));
	if (current_level_ == maximum_level_)
		return m2.volume_level(d);

	std::vector<Dart> v1, v2;
	foreach_dart_of_orbit(*this, Volume(old), [&](Dart it) -> bool {
		v1.push_back(it);
		return true;
	});

	auto fn_sort = [](Dart d, Dart dd) -> bool { return d.index > dd.index; };
	std::sort(v1.begin(), v1.end(), fn_sort);

	uint32 result = UINT_MAX;
	uint32 i = 0, j = 0;
	do
	{
		i = 0;
		j = 0;
		v2.clear();
		foreach_dart_of_orbit(m2, Volume(old), [&](Dart it) -> bool {
			v2.push_back(it);
			return true;
		});
		std::sort(v2.begin(), v2.end(), fn_sort);
		while (i < v1.size() && j < v2.size())
		{
			if (v1[i] == v2[j])
			{
				j++;
			}
			i++;
		}
		if (j == v2.size())
			result = m2.volume_level(old);
		m2.current_level_++;
	} while (j != v2.size());

	return result;
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

bool EMR_Map3_Adaptative::disable_edge_subdivision(Edge e, bool disable_neighbor)
{
	return false;
}
bool EMR_Map3_Adaptative::disable_face_subdivision(Face f, bool disable_edge, bool disable_subface)
{
	return false;
}
bool EMR_Map3_Adaptative::disable_volume_subdivision(Volume v, bool disable_face)
{
	return false;
}

} // namespace cgogn
