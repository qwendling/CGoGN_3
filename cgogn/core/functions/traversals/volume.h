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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_volume(const MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_incident_volume(m, c, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func, CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Volume = typename mesh_traits<MESH>::Volume;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		func(Volume(c.dart));
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
	{
		Dart d = c.dart;
		do
		{
			if (!is_boundary(m, d))
			{
				if (!func(Volume(d)))
					break;
			}
			d = phi3(m, phi2(m, d));
		} while (d != c.dart);
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		Dart d = c.dart;
		if (!is_boundary(m, d))
			if (!func(Volume(d)))
				return;
		d = phi3(m, d);
		if (!is_boundary(m, d))
			func(Volume(d));
	}
	else
	{
		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Volume>(m))
		{
			CellMarkerStore<MESH, Volume> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Volume v(d);
				if constexpr (mesh_traits<MESH>::dimension == 3) // volumes can be boundary cells
				{
					if (!is_boundary(m, d) && !marker.is_marked(v))
					{
						marker.mark(v);
						return func(v);
					}
				}
				else
				{
					if (!marker.is_marked(v))
					{
						marker.mark(v);
						return func(v);
					}
				}
				return true;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				if constexpr (mesh_traits<MESH>::dimension == 3) // volumes can be boundary cells
				{
					if (!is_boundary(m, d) && !marker.is_marked(d))
					{
						Volume v(d);
						foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(v);
					}
				}
				else
				{
					if (!marker.is_marked(d))
					{
						Volume v(d);
						foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(v);
					}
				}
				return true;
			});
		}
	}
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Volume> incident_volumes(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::Volume> incident_volumes(const MESH& m, CELL c)
{
	using Volume = typename mesh_traits<MESH>::Volume;
	if constexpr (mesh_traits<MESH>::dimension == 2)
		return {Volume(c.dart)};
	else
	{
		std::vector<Volume> volumes;
		volumes.reserve(32u);
		foreach_incident_volume(m, c, [&](Volume v) -> bool {
			volumes.push_back(v);
			return true;
		});
		return volumes;
	}
}

template <typename MESH, typename CELL, typename FUNC>
void foreach_adjacent_volume_through_vertex(const MESH& m, CELL v, const FUNC& func)
{
	using Volume = typename mesh_traits<MESH>::Volume;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	if (is_indexed<Volume>(m))
	{
		CellMarkerStore<MESH, Volume> marker_volume(m);
		marker_volume.mark(v);
		foreach_incident_vertex(m, v, [&](Vertex inc_vert) -> bool {
			bool res_nested_lambda = true;
			foreach_incident_volume(m, inc_vert, [&](Volume inc_vol) -> bool {
				if (!marker_volume.is_marked(inc_vol) && !is_boundary(m, inc_vol.dart))
				{
					marker_volume.mark(inc_vol);
					res_nested_lambda = func(inc_vol);
				}
				return res_nested_lambda;
			});
			return res_nested_lambda;
		});
	}
	else
	{
		/*DartMarkerStore<MESH> marker_volume(m);
		marker_volume.
		marker_volume.mark_orbit(v);
		foreach_incident_vertex(v, [&] (Vertex inc_vert)->bool
		{
			bool res_nested_lambda = true;
			foreach_incident_volume(inc_vert, [&](Volume inc_vol)->bool
			{
				if (!marker_volume.is_marked(inc_vol.dart) && !is_boundary(inc_vol.dart))
				{
					marker_volume.mark_orbit(inc_vol);
					res_nested_lambda = func(inc_vol);
				}
				return res_nested_lambda;
			});
			return res_nested_lambda;
		});*/
	}
}

template <typename MESH, typename CELL, typename FUNC>
inline void foreach_adjacent_volume_through_edge(const MESH& m, CELL v, const FUNC& func)
{
	using Volume = typename mesh_traits<MESH>::Volume;
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	if (is_indexed<Volume>(m))
	{
		CellMarkerStore<MESH, Volume> marker_volume(m);
		marker_volume.mark(v);
		foreach_incident_edge(m, v, [&](Edge inc_edge) -> bool {
			bool res_nested_lambda = true;
			foreach_incident_volume(m, inc_edge, [&](Volume inc_vol) -> bool {
				if (!marker_volume.is_marked(inc_vol) && !is_boundary(m, inc_vol.dart))
				{
					marker_volume.mark(inc_vol);
					res_nested_lambda = func(inc_vol);
				}
				return res_nested_lambda;
			});
			return res_nested_lambda;
		});
	}
	else
	{
		/*DartMarkerStore marker_volume(m);
		marker_volume.mark_orbit(v);
		foreach_incident_edge(v, [&] (Edge inc_edge)->bool
		{
			bool res_nested_lambda = true;
			foreach_incident_volume(inc_edge, [&](Volume inc_vol)->bool
			{
				if (!marker_volume.is_marked(inc_vol.dart) && !is_boundary(inc_vol.dart))
				{
					marker_volume.mark_orbit(inc_vol);
					res_nested_lambda = func(inc_vol);
				}
				return res_nested_lambda;
			});
			return res_nested_lambda;
		});*/
	}
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_
