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

uint32 EMR_Map3_Adaptative::get_dart_lookup(Dart d) const
{
	return std::max(this->dart_level(d), (*dart_lookup_)[d.index]);
}

void EMR_Map3_Adaptative::set_dart_lookup(Dart d, uint32 v)
{
	(*dart_lookup_)[d.index] = v;
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
	Dart it = phi1(*this, d);
	Dart result = d;
	while (it != d)
	{
		if (m_.dart_level(it) > m_.dart_level(result))
			result = it;
		it = phi1(*this, it);
	}

	return result;
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
	EMR_Map3_Adaptative m2(m_);
	m2.current_level_ = current_level_ + 1;
	Dart it = phi1(*this, d);
	Dart it2 = phi1(m2, d);
	while (it != d && it2 != d)
	{
		if (it == it2)
		{
			it = phi1(*this, it);
		}
		it2 = phi1(m2, it2);
	}
	return it == d;
}

uint32 EMR_Map3_Adaptative::face_level(Dart d) const
{
	uint32 visibility = get_dart_visibility(d);
	cgogn_message_assert(visibility <= current_level_, "Access to a dart introduced after current level");
	if (get_dart_visibility(d) == 0)
		return 0;
	Dart old = face_oldest_dart(d);
	EMR_Map3_Adaptative m2(m_);
	m2.current_level_ = current_level_ - 1;
	while (m2.current_level_ > 0 && !m2.face_is_subdivided(old))
	{
		m2.current_level_--;
	}
	return m2.current_level_ + 1;
}

} // namespace cgogn
