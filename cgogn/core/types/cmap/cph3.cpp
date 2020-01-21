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

#include <cgogn/core/functions/traversals/face.h>

namespace cgogn
{

/***************************************************
 *              LEVELS MANAGEMENT                  *
 ***************************************************/

uint32 CPH3::dart_level(Dart d) const
{
	return (*dart_level_)[d.index];
}

void CPH3::set_dart_level(Dart d, uint32 l)
{
	if(nb_darts_per_level_.size() > dart_level(d))
		nb_darts_per_level_[dart_level(d)]--;
	if(nb_darts_per_level_.size() < l)
		nb_darts_per_level_.resize(l);
	nb_darts_per_level_[l]++;
	if (l > dart_level(d) && l > maximum_level_)
		maximum_level_ = l;
	if (l < dart_level(d))
	{
		while (nb_darts_per_level_[maximum_level_] == 0u)
		{
			--maximum_level_;
			nb_darts_per_level_.pop_back();
		}
	}
	(*dart_level_)[d.index] = l;
}

/***************************************************
 *             EDGE ID MANAGEMENT                  *
 ***************************************************/

uint32 CPH3::edge_id(Dart d) const
{
	return (*edge_id_)[d.index];
}

void CPH3::set_edge_id(Dart d, uint32 i)
{
	(*edge_id_)[d.index] = i;
}

uint32 CPH3::refinement_edge_id(Dart d, Dart e) const
{
	uint32 d_id = edge_id(d);
	uint32 e_id = edge_id(e);

	uint32 id = d_id + e_id;

	if (id == 0u)
		return 1u;
	else if (id == 1u)
		return 2u;
	else if (id == 2u)
	{
		if (d_id == e_id)
			return 0u;
		else
			return 1u;
	}
	// else if (id == 3)
	return 0u;
}

/***************************************************
 *             FACE ID MANAGEMENT                  *
 ***************************************************/

uint32 CPH3::face_id(Dart d) const
{
	return (*face_id_)[d.index];
}

void CPH3::set_face_id(Dart d, uint32 i)
{
	(*face_id_)[d.index] = i;
}

uint32 CPH3::refinement_face_id(const std::vector<Dart>& cut_path) const
{
	std::unordered_set<uint32> set_fid;
	for (Dart d : cut_path)
		set_fid.insert(face_id(d));
	uint32 result = 0;
	while (set_fid.find(result) != set_fid.end())
		++result;
	return result;
}

} // namespace cgogn
