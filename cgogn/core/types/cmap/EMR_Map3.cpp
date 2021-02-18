#include "EMR_Map3.h"
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

bool EMR_Map3::check_integrity() const
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

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

Dart EMR_Map3::edge_youngest_dart(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = phi2(*this, d);
	if (m_.dart_level(d) > m_.dart_level(it))
		return d;
	return it;
}

uint32 EMR_Map3::edge_level(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	return m_.dart_level(edge_youngest_dart(d));
}

bool EMR_Map3::edge_is_subdivided(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_++;
	if (phi1(*this, d) == phi1(m2, d))
	{
		return false;
	}
	return true;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

Dart EMR_Map3::face_youngest_dart(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
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

Dart EMR_Map3::face_oldest_dart(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
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

bool EMR_Map3::face_is_subdivided(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
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
	return it != d;
}

uint32 EMR_Map3::face_level(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == 0)
		return 0;
	Dart old = face_oldest_dart(d);
	EMR_Map3 m2(m_);
	m2.current_level_ = current_level_ - 1;
	while (m2.current_level_ > dart_level(old) && !m2.face_is_subdivided(old))
	{
		m2.current_level_--;
	}
	if (m2.current_level_ == dart_level(old) && !m2.face_is_subdivided(old))
	{
		return m2.current_level_;
	}
	return m2.current_level_ + 1;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

Dart EMR_Map3::volume_youngest_dart(Dart d) const
{
	Dart result = d;
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		if (m_.dart_level(it) > m_.dart_level(result))
		{
			result = it;
		}
		return true;
	});
	return result;
}

Dart EMR_Map3::volume_oldest_dart(Dart d) const
{
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
bool EMR_Map3::volume_is_subdivided(Dart d) const
{
	cgogn_message_assert(m_.dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = current_level_ + 1;
	std::vector<Dart> v1, v2;
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		v1.push_back(it);
		return true;
	});
	foreach_dart_of_orbit(m2, Volume(d), [&](Dart it) -> bool {
		v2.push_back(it);
		return true;
	});
	auto fn_sort = [](Dart d, Dart dd) -> bool { return d.index > dd.index; };
	std::sort(v1.begin(), v1.end(), fn_sort);
	std::sort(v2.begin(), v2.end(), fn_sort);

	uint32 i = 0, j = 0;
	while (i < v1.size() && j < v2.size())
	{
		if (v1[i] == v2[j])
		{
			i++;
		}
		j++;
	}

	return i < v1.size();
}

uint32 EMR_Map3::volume_level(Dart d) const
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
