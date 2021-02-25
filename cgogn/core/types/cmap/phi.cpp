/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/core/types/cmap/cph3.h>
#include <cgogn/core/types/cmap/phi.h>

namespace cgogn
{

Dart phi1(const EMR_Map3& m, Dart d)
{
	cgogn_message_assert(m.m_.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	return (*((*m.m_.MR_phi1_)[m.current_level_]))[d.index];
}
Dart phi_1(const EMR_Map3& m, Dart d)
{
	cgogn_message_assert(m.m_.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	return (*((*m.m_.MR_phi_1_)[m.current_level_]))[d.index];
}
Dart phi2(const EMR_Map3& m, Dart d)
{
	cgogn_message_assert(m.m_.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	return (*((*m.m_.MR_phi2_)[m.current_level_]))[d.index];
}
Dart phi3(const EMR_Map3& m, Dart d)
{
	cgogn_message_assert(m.m_.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	return (*((*m.m_.MR_phi3_)[m.current_level_]))[d.index];
}

Dart phi1(const EMR_Map3_Adaptative& m, Dart d)
{
	cgogn_message_assert(m.get_dart_visibility(d) <= m.current_level_, "Access to a dart visible at a higher level");
	EMR_Map3 emr = EMR_Map3(m);

	emr.current_level_ = m.dart_level(d);
	Dart result = phi1(emr, d);
	bool is_found = true;
	while (emr.current_level_ != m.maximum_level_ && is_found)
	{
		is_found = false;
		emr.current_level_++;
		Dart tmp = phi1(emr, d);
		if (m.get_dart_visibility(tmp) > m.current_level_)
		{
			tmp = phi_1(emr, result);
			if (m.get_dart_visibility(tmp) <= m.current_level_)
			{
				if (tmp != result)
					is_found = true;
			}
		}
		else
		{
			if (tmp != result)
				is_found = true;
		}
		if (is_found)
			result = tmp;
	}

	return result;
}
Dart phi_1(const EMR_Map3_Adaptative& m, Dart d)
{
	return phi3(m, phi1(m, phi3(m, d)));
}
Dart phi2(const EMR_Map3_Adaptative& m, Dart d)
{
	cgogn_message_assert(m.get_dart_visibility(d) <= m.current_level_, "Access to a dart visible at a higher level");
	EMR_Map3 emr = EMR_Map3(m);
	emr.current_level_ = std::max(m.current_level_, m.dart_level(d));
	Dart result = phi2(emr, d);
	bool is_found = true;
	while (emr.current_level_ != m.maximum_level_ && is_found)
	{
		is_found = false;
		emr.current_level_++;
		Dart tmp = phi3(emr, d);
		Dart it = tmp;
		do
		{
			it = phi2(emr, phi3(emr, it));
		} while (it != tmp && m.get_dart_visibility(it) > m.current_level_);
		if (m.get_dart_visibility(it) <= m.current_level_ && tmp != result)
		{
			result = it;
			is_found = true;
		}
	}

	return result;
}
Dart phi3(const EMR_Map3_Adaptative& m, Dart d)
{
	cgogn_message_assert(m.get_dart_visibility(d) <= m.current_level_, "Access to a dart visible at a higher level");
	EMR_Map3 emr = EMR_Map3(m);
	Dart result;

	for (uint32 i = m.maximum_level_; i >= m.dart_level(d); --i)
	{
		emr.current_level_ = i;
		result = phi3(emr, d);
		if (m.get_dart_visibility(result) <= m.current_level_)
		{
			return result;
		}
	}
	return result;
}

Dart phi2bis(const CPH3& m, Dart d)
{
	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	uint32 face_id = m.face_id(d);
	Dart it = d;

	it = phi2(map, it);

	if (m.face_id(it) == face_id)
		return it;
	else
	{
		do
		{
			it = phi2(map, phi3(map, it));
		} while (m.face_id(it) != face_id);

		return it;
	}
}

Dart phi1(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	if (m.current_level_ == m.maximum_level_)
		return phi1(map, d);

	bool finished = false;
	uint32 edge_id = m.edge_id(d);
	Dart it = d;
	do
	{
		it = phi1(map, it);
		if (m.dart_level(it) <= m.current_level_)
			finished = true;
		else
			while (m.edge_id(it) != edge_id)
				it = phi1(map, phi2bis(m, it));
	} while (!finished);
	return it;
}

Dart phi_1(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	bool finished = false;
	Dart it = phi_1(map, d);
	uint32 edge_id = m.edge_id(it);

	do
	{
		if (m.dart_level(it) <= m.current_level_)
			finished = true;
		else
		{
			it = phi_1(map, it);
			while (m.edge_id(it) != edge_id)
				it = phi_1(map, phi2bis(m, it));
		}
	} while (!finished);
	return it;
}

Dart phi2(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);
	return phi2(map, phi_1(map, phi1(m, d)));
}

Dart phi3(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);
	if (phi3(map, d) == d)
		return d;
	return phi3(map, phi_1(map, phi1(m, d)));
}

Dart phi1(const CPH3_adaptative& m, Dart d)
{
	cgogn_message_assert(m.dart_is_visible(d), "Access to a dart not visible at this level");

	Dart it;
	if (m.get_phi1_buffer(d, it))
	{
		return it;
	}

	const CPH3_adaptative::CMAP& map = static_cast<const CPH3_adaptative::CMAP&>(m);
	if (m.maximum_level_ == m.current_level_)
		return phi1(map, d);

	uint32 edge_id = m.edge_id(d);
	it = phi1(map, d);
	do
	{
		if (m.dart_is_visible(it))
			break;
		else
			while (m.edge_id(it) != edge_id)
				it = phi1(map, phi2bis(m, it));
		if (m.dart_is_visible(it))
			break;
		it = phi1(map, it);
	} while (1);
	m.set_phi1_buffer(d, it);
	return it;
}

Dart phi_1(const CPH3_adaptative& m, Dart d)
{
	cgogn_message_assert(m.dart_is_visible(d), "Access to a dart not visible at this level");
	Dart it = phi1(m, d);
	Dart it2 = phi1(m, it);
	while (it2 != d)
	{
		it = it2;
		it2 = phi1(m, it2);
	}
	return it;
}

Dart phi2(const CPH3_adaptative& m, Dart d)
{
	cgogn_message_assert(m.dart_is_visible(d), "Access to a dart not visible at this level");
	Dart it;
	if (m.get_phi2_buffer(d, it))
	{
		return it;
	}
	it = phi1(m, d);
	CPH3 m2(m);
	m2.current_level_ = std::max(m2.dart_level(it), m2.dart_level(d));
	it = phi2(m2, d);
	while (!m.dart_is_visible(it))
	{
		it = phi<32>(m2, it);
	}
	m.set_phi2_buffer(d, it);
	return it;
	/*it = phi_1(map,phi1(m2,d));
	Dart it2 = phi2(map,it);
	if(m.dart_is_visible(it2))
		return it2;
	return phi2bis(m,it);*/
}

Dart phi3(const CPH3_adaptative& m, Dart d)
{
	cgogn_message_assert(m.dart_is_visible(d), "Access to a dart not visible at this level");
	const CPH3_adaptative::CMAP& map = static_cast<const CPH3_adaptative::CMAP&>(m);
	if (phi3(map, d) == d)
		return d;

	Dart it;
	if (m.get_phi3_buffer(d, it))
	{
		return it;
	}
	it = phi1(m, d);
	CPH3_adaptative m2(m);
	if (m2.current_level_ < m2.dart_level(it))
		m2.current_level_ = m2.dart_level(it);
	it = phi3(map, phi_1(map, phi1(m2, d)));
	m.set_phi3_buffer(d, it);
	return it;
}

} // namespace cgogn
