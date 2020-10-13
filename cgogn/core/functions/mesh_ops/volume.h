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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_OPS_VOLUME_H_
#define CGOGN_CORE_FUNCTIONS_MESH_OPS_VOLUME_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// add_pyramid(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Volume CGOGN_CORE_EXPORT add_pyramid(CMap2& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// add_prism(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Volume CGOGN_CORE_EXPORT add_prism(CMap2& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// cut_volume(MESH& m, const std::vector<Dart>& path, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap3 //
///////////

CMap3::Face cut_volume(CMap3& m, const std::vector<Dart>& path, bool set_indices = true);

//////////
// CPH3 //
//////////

CPH3::CMAP::Face cut_volume(CPH3& m, const std::vector<Dart>& path, bool set_indices = true);

/////////////////////
// CPH3_adaptative //
/////////////////////
CPH3_adaptative::CMAP::Face cut_volume(CPH3_adaptative& m, const std::vector<Dart>& path, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// close_hole(MESH& m, Dart d, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap3 //
///////////

CMap3::Volume close_hole(CMap3& m, Dart d, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// uint32
// close(MESH& m, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap3 //
///////////

uint32 close(CMap3& m, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// void
// unsew_volume(MESH& m);

/*****************************************************************************/

///////////
// CMap3 //
///////////

template <typename FUNC>
void unsew_volume(CMap3& m, const mesh_traits<CMap3>::Face f, const FUNC& callback_vertices, bool set_indices = true)
{
	using Vertex = typename mesh_traits<CMap3>::Vertex;
	using Face = typename mesh_traits<CMap3>::Face;
	using Face2 = typename mesh_traits<CMap3>::Face2;

	static_assert(is_func_parameter_same<FUNC, std::pair<Vertex, Vertex>>::value,
				  "Function must have std::pair<Vertex, Vertex> as a parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_incident_to_boundary(m, f))
	{
		return;
	}

	auto same_orbit = [&](auto v1, auto v2) -> bool {
		bool result = false;
		foreach_dart_of_orbit(m, v1, [&](Dart d) -> bool {
			if (v2.dart == d)
			{
				result = true;
				return false;
			}
			return true;
		});
		return result;
	};

	std::vector<std::pair<Vertex, Vertex>> list_pair_vertex;
	foreach_dart_of_orbit(m, Face2(f.dart), [&](Dart d) -> bool {
		list_pair_vertex.push_back({Vertex(d), Vertex(phi<31>(m, d))});
		Dart tmp = d;
		do
		{
			tmp = phi<23>(m, tmp);
			if (is_boundary(m, tmp))
			{
				phi2_unsew(m, tmp);
				break;
			}
		} while (tmp != d);
		phi3_unsew(m, d);
		return true;
	});

	close_hole(m, f.dart, false);
	foreach_dart_of_orbit(m, Face2(phi3(m, f.dart)), [&](Dart d) -> bool {
		set_boundary(m, d, true);
		return true;
	});
	if (list_pair_vertex[0].second.dart == phi3(m, list_pair_vertex[0].second.dart))
	{
		close_hole(m, list_pair_vertex[0].second.dart, false);
	}
	foreach_dart_of_orbit(m, Face2(phi3(m, list_pair_vertex[0].second.dart)), [&](Dart d) -> bool {
		set_boundary(m, d, true);
		return true;
	});
	if (set_indices)
	{
		std::pair<Vertex, Vertex> pv = list_pair_vertex[0];
		Face f1 = Face(pv.first.dart);
		Face f2 = Face(pv.second.dart);
		Dart d = f1.dart;
		Dart it = d;
		Dart it2 = f2.dart;
		set_index(m, f2, new_index<Face>(m));
		do
		{

			Dart it_3 = phi3(m, it);
			Dart it2_3 = phi3(m, it2);
			if (is_indexed<CMap3::Vertex>(m))
			{
				if (!same_orbit(Vertex(it), Vertex(it2)))
				{
					auto tmp = new_index<Vertex>(m);
					set_index(m, Vertex(it2), tmp);
				}
				else
				{
					copy_index<CMap3::Vertex>(m, it2_3, phi1(m, it2));
				}
				copy_index<CMap3::Vertex>(m, it_3, phi1(m, it));
			}
			if (is_indexed<CMap3::Edge>(m))
			{
				if (!same_orbit(CMap3::Edge(it), CMap3::Edge(it2)))
				{
					set_index(m, CMap3::Edge(it2), new_index<CMap3::Edge>(m));
				}
				else
				{
					copy_index<CMap3::Edge>(m, it2_3, phi1(m, it2));
				}
				copy_index<CMap3::Edge>(m, it_3, phi1(m, it));
			}
			if (is_indexed<CMap3::Face>(m))
			{
				copy_index<CMap3::Face>(m, it_3, phi1(m, it));
			}
			it = phi1(m, it);
			it2 = phi_1(m, it2);
		} while (it != d);
	}
	for (auto p : list_pair_vertex)
	{
		callback_vertices(p);
	}
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_VOLUME_H_
