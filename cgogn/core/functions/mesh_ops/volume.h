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

#include <cgogn/core/types/cmap/EMR_Map3.h>
#include <cgogn/core/types/cmap/EMR_Map3_Adaptative.h>
#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/cmap/cph3.h>
#include <cgogn/core/types/cmap/cph3_adaptative.h>
#include <cgogn/core/types/cmap/graph.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>

#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/functions/mesh_info.h>

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
// typename mesh_traits<MESH>::Volume
// remove_volume(MESH& m, typename mesh_traits<MESH>::Volume v);

/*****************************************************************************/

///////////
// CMap2 //
///////////

void CGOGN_CORE_EXPORT remove_volume(CMap2& m, CMap2::Volume v);

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

//////////////
// EMR_Map3 //
//////////////

EMR_Map3::Face cut_volume(EMR_Map3& m, const std::vector<Dart>& path, bool set_indices = true);

/////////////////////////
// EMR_Map3_Adaptative //
/////////////////////////

EMR_Map3_Adaptative::Face cut_volume(EMR_Map3_Adaptative& m, const std::vector<Dart>& path, bool set_indices = true);

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
		list_pair_vertex.push_back({Vertex(d), Vertex(phi<3, 1>(m, d))});
		Dart tmp = d;
		do
		{
			tmp = phi<2, 3>(m, tmp);
			if (is_boundary(m, tmp))
			{
				phi2_unsew(m, tmp);
				break;
			}
		} while (tmp != d);
		phi3_unsew(m, d);
		return true;
	});
	// close(m, false);
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

//////////////
// EMR_Map3 //
//////////////

template <typename FUNC>
void unsew_volume(EMR_Map3& m, const mesh_traits<EMR_Map3>::Face f, const FUNC& callback_vertices,
				  bool set_indices = true)
{
	using Face = typename mesh_traits<EMR_Map3>::Face;
	EMR_Map3 m2(m);
	uint32 tmp_f_level = m.face_level(f.dart);
	Dart tmp_old = m.face_oldest_dart(f.dart);
	if (m.dart_level(tmp_old) != tmp_f_level)
	{
		m2.current_level_ = m.dart_level(tmp_old);
		unsew_volume(m2, Face(tmp_old), callback_vertices, set_indices);
		return;
	}
	unsew_volume_aux(m, f, callback_vertices, set_indices);
}

template <typename FUNC>
void unsew_volume_aux(EMR_Map3& m, const mesh_traits<EMR_Map3>::Face f, const FUNC& callback_vertices,
					  bool set_indices = true)
{
	using Vertex = typename mesh_traits<EMR_Map3>::Vertex;
	using Edge = typename mesh_traits<EMR_Map3>::Edge;
	using Face = typename mesh_traits<EMR_Map3>::Face;

	static_assert(is_func_parameter_same<FUNC, std::pair<Vertex, Vertex>>::value,
				  "Function must have std::pair<Vertex, Vertex> as a parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_incident_to_boundary(m, f))
	{
		return;
	}

	Dart f_rep = f.dart;
	Dart f3_rep = phi3(m, f_rep);
	std::pair<Vertex, Vertex> p_rep = {Vertex(f_rep), Vertex(phi3(m, f_rep))};
	EMR_Map3 m2(m);
	m2.current_level_ = std::min(m.current_level_ + 1, m.maximum_level_);

	auto same_orbit = [&](auto mo, auto v1, auto v2) -> bool {
		bool result = false;
		foreach_dart_of_orbit(mo, v1, [&](Dart d) -> bool {
			if (v2.dart == d)
			{
				result = true;
				return false;
			}
			return true;
		});
		return result;
	};

	if (m.current_level_ == m.maximum_level_)
	{
		std::vector<std::pair<Vertex, Vertex>> list_pair_vertex;
		auto get_list_pair_vertex = [&list_pair_vertex](std::pair<Vertex, Vertex> p) -> bool {
			list_pair_vertex.push_back(p);
			return true;
		};
		unsew_volume(*m.get_map(), f, get_list_pair_vertex, false);
		if (set_indices && is_indexed<CMap3::Vertex>(m))
		{
			Dart it = f_rep;
			Dart it2 = phi1(m, f3_rep);
			do
			{
				Dart it_3 = phi3(m, it);
				Dart it2_3 = phi3(m, it2);
				if (!same_orbit(m, Vertex(it), Vertex(it2)))
				{
					auto tmp = new_index<Vertex>(m);
					set_index(m, Vertex(it2), tmp);
				}
				else
				{
					copy_index<CMap3::Vertex>(m, phi1(m, it2_3), it2);
				}
				copy_index<CMap3::Vertex>(m, phi1(m, it_3), it);
				it = phi1(m, it);
				it2 = phi_1(m, it2);
			} while (it != f_rep);
		}

		for (auto p : list_pair_vertex)
		{
			callback_vertices(p);
		}
		Dart it = f_rep;
		Dart it3 = f3_rep;
		do
		{
			m.set_dart_level(phi3(m, it), m.dart_level(it3));
			m.set_dart_level(phi3(m, it3), m.dart_level(it));
			it = phi1(m, it);
			it3 = phi_1(m, it3);
		} while (it != f_rep);
	}
	else
	{
		std::vector<Dart> vect_dart_unsew;
		Dart it = f_rep;
		do
		{
			vect_dart_unsew.push_back(it);
			it = phi1(m, it);
		} while (it != f_rep);
		for (auto it : vect_dart_unsew)
		{
			unsew_volume_aux(m2, Face(it), callback_vertices, set_indices);
		}
	}

	std::array<Dart, 2> ar_dart = {m.face_oldest_dart(f_rep), m.face_oldest_dart(f3_rep)};
	int nb_dart = 0;
	Dart it;
	while (nb_dart < 2)
	{
		it = ar_dart[nb_dart];
		// phi3
		do
		{
			Dart d3 = phi<23>(m2, phi2(m, it));
			(*((*m.m_.MR_phi3_)[m.current_level_]))[it.index] = d3;
			(*((*m.m_.MR_phi3_)[m.current_level_]))[d3.index] = it;
			it = phi1(m, it);
		} while (it != ar_dart[nb_dart]);
		// phi2
		do
		{

			Dart tmp = phi<32>(m2, it);
			Dart d3 = phi3(m, it);

			(*((*m.m_.MR_phi2_)[m.current_level_]))[tmp.index] = d3;

			(*((*m.m_.MR_phi2_)[m.current_level_]))[d3.index] = tmp;

			it = phi1(m, it);
		} while (it != ar_dart[nb_dart]);
		// phi1
		do
		{
			Dart d3_1 = phi<13>(m, it);
			Dart d3 = phi3(m, it);
			(*((*m.m_.MR_phi1_)[m.current_level_]))[d3_1.index] = d3;
			(*((*m.m_.MR_phi_1_)[m.current_level_]))[d3.index] = d3_1;
			it = phi1(m, it);
		} while (it != ar_dart[nb_dart]);
		nb_dart++;
	}
	if (set_indices)
	{
		Face f1 = Face(p_rep.first.dart);
		Face f2 = Face(p_rep.second.dart);
		Dart d = f1.dart;
		Dart it = d;
		Dart it2 = f2.dart;
		uint32 new_id = new_index<Face>(m);
		uint32 f1_id = index_of(m, Face(d));
		uint32 f_level = m.face_level(f.dart);
		do
		{
			if (is_indexed<Edge>(m))
			{
				if (!same_orbit(m, Edge(it), Edge(it2)))
				{
					set_index<Edge>(m, it2, new_index<Edge>(m));
				}
				uint32 e_level = m.edge_level(it2);
				foreach_dart_of_orbit(m, Edge(it2), [&](Dart dd) -> bool {
					if (m.dart_level(dd) == e_level)
					{
						set_index<Edge>(m, dd, index_of(m, Edge(it2)));
					}
					return true;
				});
				foreach_dart_of_orbit(m, Edge(it), [&](Dart dd) -> bool {
					if (m.dart_level(dd) == e_level)
					{
						set_index<Edge>(m, dd, index_of(m, Edge(it)));
					}
					return true;
				});
			}
			if (is_indexed<Face>(m))
			{

				if (m.dart_level(it2) >= f_level)
				{
					set_index<Face>(m, it2, new_id);
				}
				if (m.dart_level(phi3(m, it2)) >= f_level)
				{
					set_index<Face>(m, phi3(m, it2), new_id);
				}
				if (m.dart_level(phi3(m, it)) >= f_level)
				{
					set_index<Face>(m, phi3(m, it), f1_id);
				}
			}
			it = phi1(m, it);
			it2 = phi_1(m, it2);
		} while (it != d);
	}
	m.m_.clock_++;
}

/////////////////////////
// EMR_Map3_Adaptative //
/////////////////////////

template <typename FUNC>
void unsew_volume(EMR_Map3_Adaptative& m, const mesh_traits<EMR_Map3_Adaptative>::Face f, const FUNC& callback_vertices,
				  bool set_indices = true)
{
	using Vertex = typename mesh_traits<EMR_Map3_Adaptative>::Vertex;
	using Face = typename mesh_traits<EMR_Map3_Adaptative>::Face;

	static_assert(is_func_parameter_same<FUNC, std::pair<Vertex, Vertex>>::value,
				  "Function must have std::pair<Vertex, Vertex> as a parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_incident_to_boundary(m, f))
	{
		return;
	}

	uint32 f_level = m.face_level(f.dart);
	/*EMR_Map3_Adaptative m2(m.m_);
	m2.current_level_ = f_level;
	unsew_volume_aux(m2, Face(m.face_oldest_dart(f.dart)), callback_vertices, set_indices);*/
	EMR_Map3_Adaptative* topo = &m;
	if (m.topology_)
		topo = m.topology_;

	uint32 cur = topo->current_level_;
	topo->current_level_ = f_level;

	unsew_volume_aux(*topo, Face(m.face_oldest_dart(f.dart)), callback_vertices, set_indices);

	topo->current_level_ = cur;

	cgogn_message_assert(m.check_integrity(), "check integrity fail after cut");
}

template <typename FUNC>
void unsew_volume_aux(EMR_Map3_Adaptative& m, const mesh_traits<EMR_Map3>::Face f, const FUNC& callback_vertices,
					  bool set_indices = true)
{
	using Vertex = typename mesh_traits<EMR_Map3>::Vertex;
	using Edge = typename mesh_traits<EMR_Map3>::Edge;
	using Face = typename mesh_traits<EMR_Map3>::Face;

	static_assert(is_func_parameter_same<FUNC, std::pair<Vertex, Vertex>>::value,
				  "Function must have std::pair<Vertex, Vertex> as a parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_incident_to_boundary(m, f))
	{
		return;
	}

	Dart f_rep = f.dart;
	Dart f3_rep = phi3(m, f_rep);
	std::pair<Vertex, Vertex> p_rep = {Vertex(f_rep), Vertex(phi3(m, f_rep))};
	/*EMR_Map3_Adaptative m2(m.m_);
	m2.current_level_ = std::min(m.current_level_ + 1, m.maximum_level_);*/
	std::array<Dart, 2> ar_dart = {m.face_oldest_dart(f_rep), m.face_oldest_dart(f3_rep)};

	auto same_orbit = [&](auto mo, auto v1, auto v2) -> bool {
		bool result = false;
		foreach_dart_of_orbit(mo, v1, [&](Dart d) -> bool {
			if (v2.dart == d)
			{
				result = true;
				return false;
			}
			return true;
		});
		return result;
	};

	if (m.current_level_ == m.maximum_level_)
	{
		std::vector<std::pair<Vertex, Vertex>> list_pair_vertex;
		auto get_list_pair_vertex = [&list_pair_vertex](std::pair<Vertex, Vertex> p) -> bool {
			list_pair_vertex.push_back(p);
			return true;
		};
		unsew_volume(*m.get_map(), f, get_list_pair_vertex, false);
		m.m_.clock_++;
		if (set_indices && is_indexed<CMap3::Vertex>(m))
		{
			Dart it = f_rep;
			Dart it2 = phi1(m, f3_rep);
			do
			{
				Dart it_3 = phi3(m, it);
				Dart it2_3 = phi3(m, it2);
				if (!same_orbit(m, Vertex(it), Vertex(it2)))
				{
					auto tmp = new_index<Vertex>(m);
					set_index(m, Vertex(it2), tmp);
				}
				else
				{
					copy_index<CMap3::Vertex>(m, phi1(m, it2_3), it2);
				}
				copy_index<CMap3::Vertex>(m, phi1(m, it_3), it);
				it = phi1(m, it);
				it2 = phi_1(m, it2);
			} while (it != f_rep);
		}
		m.m_.clock_++;

		for (auto p : list_pair_vertex)
		{
			callback_vertices(p);
		}
		Dart it = f_rep;
		Dart it3 = f3_rep;
		do
		{
			m.set_dart_level(phi3(m, it), m.dart_level(it3));
			m.set_dart_level(phi3(m, it3), m.dart_level(it));
			// m.set_dart_visibility(phi3(m, it), m.get_dart_visibility_fast(it3));
			// m.set_dart_visibility(phi3(m, it3), m.get_dart_visibility_fast(it));
			m.m_.clock_++;
			it = phi1(m, it);
			it3 = phi_1(m, it3);
		} while (it != f_rep);
	}
	else
	{
		std::vector<Dart> vect_dart_unsew;
		Dart it = f_rep;
		do
		{
			vect_dart_unsew.push_back(it);
			it = phi1(m, it);
		} while (it != f_rep);
		uint32 cur = m.current_level_;
		m.current_level_ = std::min(cur + 1, m.maximum_level_);
		for (auto it : vect_dart_unsew)
		{
			unsew_volume_aux(m, Face(it), callback_vertices, set_indices);
		}
		m.current_level_ = cur;
	}
	m.m_.clock_++;

	EMR_Map3_Adaptative m2(m.m_);
	m2.current_level_ = std::min(m.current_level_ + 1, m.maximum_level_);

	int nb_dart = 0;
	Dart it;
	while (nb_dart < 2)
	{
		it = ar_dart[nb_dart];
		// phi3
		do
		{
			Dart d3 = phi<2, 3>(m2, phi2(m, it));
			(*((*m.m_.MR_phi3_)[m.current_level_]))[it.index] = d3;
			(*((*m.m_.MR_phi3_)[m.current_level_]))[d3.index] = it;
			m.m_.clock_++;
			it = phi1(m, it);
		} while (it != ar_dart[nb_dart]);
		// phi2
		do
		{

			// Dart tmp = phi<32>(m2, it);
			Dart tmp = phi3(m2, it);
			tmp = (*((*m.m_.MR_phi2_)[m2.current_level_]))[tmp.index];
			Dart d3 = phi3(m, it);

			(*((*m.m_.MR_phi2_)[m.current_level_]))[tmp.index] = d3;

			(*((*m.m_.MR_phi2_)[m.current_level_]))[d3.index] = tmp;

			m.m_.clock_++;

			uint32 e_level = m.edge_level(it);
			if (e_level == m2.edge_level(it))
			{
				(*((*m.m_.MR_phi2_)[e_level]))[tmp.index] = (*((*m.m_.MR_phi2_)[m.current_level_]))[tmp.index];
				(*((*m.m_.MR_phi2_)[e_level]))[d3.index] = (*((*m.m_.MR_phi2_)[m.current_level_]))[d3.index];
				// std::cout << "hello" << std::endl;
			}
			m.m_.clock_++;

			/*std::cout << "m lvl = " << m.current_level_ << std::endl;
			std::cout << "phi2 " << d3.index << " = " << tmp.index << std::endl;
			std::cout << "phi2 EMR " << d3.index << " = " << phi2(m, d3).index << std::endl;
			std::cout << "phi2 EMR " << tmp.index << " = " << phi2(m, tmp).index << std::endl;*/

			it = phi1(m, it);
		} while (it != ar_dart[nb_dart]);
		// phi1
		do
		{
			Dart d3_1 = phi<1, 3>(m, it);
			Dart d3 = phi3(m, it);
			(*((*m.m_.MR_phi1_)[m.current_level_]))[d3_1.index] = d3;
			(*((*m.m_.MR_phi_1_)[m.current_level_]))[d3.index] = d3_1;
			m.m_.clock_++;
			it = phi1(m, it);
		} while (it != ar_dart[nb_dart]);
		nb_dart++;
		// std::cout << "___________________________" << std::endl;
	}
	// std::cout << "####################" << std::endl;
	cgogn_message_assert(m.check_integrity(), "map invalid");

	if (set_indices)
	{
		Face f1 = Face(p_rep.first.dart);
		Face f2 = Face(p_rep.second.dart);
		Dart d = f1.dart;
		Dart it = d;
		Dart it2 = f2.dart;
		uint32 new_id = new_index<Face>(m);
		uint32 f1_id = index_of(m, Face(d));
		uint32 f_level = m.face_level(f.dart);
		do
		{
			if (is_indexed<Edge>(m))
			{
				if (!same_orbit(m, Edge(it), Edge(it2)))
				{
					set_index<Edge>(m, it2, new_index<Edge>(m));
				}
				uint32 e_level = m.edge_level(it2);
				foreach_dart_of_orbit(m, Edge(it2), [&](Dart dd) -> bool {
					if (m.dart_level(dd) == e_level)
					{
						set_index<Edge>(m, dd, index_of(m, Edge(it2)));
					}
					return true;
				});
				foreach_dart_of_orbit(m, Edge(it), [&](Dart dd) -> bool {
					if (m.dart_level(dd) == e_level)
					{
						set_index<Edge>(m, dd, index_of(m, Edge(it)));
					}
					return true;
				});
			}
			if (is_indexed<Face>(m))
			{

				if (m.dart_level(it2) >= f_level)
				{
					set_index<Face>(m, it2, new_id);
				}
				if (m.dart_level(phi3(m, it2)) >= f_level)
				{
					set_index<Face>(m, phi3(m, it2), new_id);
				}
				if (m.dart_level(phi3(m, it)) >= f_level)
				{
					set_index<Face>(m, phi3(m, it), f1_id);
				}
			}
			it = phi1(m, it);
			it2 = phi_1(m, it2);
		} while (it != d);
	}
	m.m_.clock_++;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_VOLUME_H_
