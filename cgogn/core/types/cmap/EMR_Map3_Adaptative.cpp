#include "EMR_Map3_Adaptative.h"
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
	if (current_level_ == 0)
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
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == 0)
		return 0;
	Dart old = volume_oldest_dart(d);
	EMR_Map3 m2(m_);
	m2.current_level_ = current_level_ - 1;
	while (m2.current_level_ > dart_level(old) && !m2.volume_is_subdivided(old))
	{
		m2.current_level_--;
	}
	if (m2.current_level_ == dart_level(old) && !m2.volume_is_subdivided(old))
	{
		return m2.current_level_;
	}
	return m2.current_level_ + 1;
}

} // namespace cgogn
