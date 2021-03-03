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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/traversals/global.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;

///////////
// CMap2 //
///////////

void hexagon_to_triangles(CMap2& m, CMap2::Face f)
{
	cgogn_message_assert(codegree(m, f) == 6, "hexagon_to_triangles: given face should have 6 edges");
	Dart d0 = phi1(m, f.dart);
	Dart d1 = phi<11>(m, d0);
	cut_face(m, CMap2::Vertex(d0), CMap2::Vertex(d1));
	Dart d2 = phi<11>(m, d1);
	cut_face(m, CMap2::Vertex(d1), CMap2::Vertex(d2));
	Dart d3 = phi<11>(m, d2);
	cut_face(m, CMap2::Vertex(d2), CMap2::Vertex(d3));
}

CMap2::Vertex quadrangulate_face(CMap2& m, CMap2::Face f)
{
	cgogn_message_assert(codegree(m, f) == 8, "quadrangulate_face: given face should have 8 edges");
	Dart d0 = phi1(m, f.dart);
	Dart d1 = phi<11>(m, d0);

	cut_face(m, CMap2::Vertex(d0), CMap2::Vertex(d1));
	cut_edge(m, CMap2::Edge(phi_1(m, d0)));

	Dart x = phi2(m, phi_1(m, d0));
	Dart dd = phi<1111>(m, x);
	while (dd != x)
	{
		Dart next = phi<11>(m, dd);
		cut_face(m, CMap2::Vertex(dd), CMap2::Vertex(phi1(m, x)));
		dd = next;
	}

	return CMap2::Vertex(phi2(m, x));
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void subdivide(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;

	CellCache<MESH> cache(m);
	cache.template build<Edge>();
	cache.template build<Face>();

	foreach_cell(cache, [&](Edge e) -> bool {
		std::vector<Vertex> vertices = incident_vertices(m, e);
		Vertex v = cut_edge(m, e);
		value<Vec3>(m, vertex_position, v) =
			0.5 * (value<Vec3>(m, vertex_position, vertices[0]) + value<Vec3>(m, vertex_position, vertices[1]));
		return true;
	});

	foreach_cell(cache, [&](Face f) -> bool {
		if (codegree(m, f) == 6)
			hexagon_to_triangles(m, f);
		return true;
	});
}

template <typename MESH, typename FUNC>
void cut_all_edges(MESH& m, const FUNC& on_edge_cut)
{
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;

	CellCache<MESH> cache(m);
	cache.template build<Edge>();

	foreach_cell(cache, [&](Edge e) -> bool {
		on_edge_cut(cut_edge(m, e));
		return true;
	});
}

template <typename MESH, typename FUNC1, typename FUNC2>
void quadrangulate_all_faces(MESH& m, const FUNC1& on_edge_cut, const FUNC2& on_face_cut)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;

	CellCache<MESH> cache(m);
	cache.template build<Edge>();
	cache.template build<Face>();

	cut_all_edges(m, on_edge_cut);

	foreach_cell(cache, [&](Face f) -> bool {
		on_face_cut(quadrangulate_face(m, f));
		return true;
	});
}

/* -------------------------- BUTTERFLY VOLUME MASKS -------------------------- */

// Remplit les vecteurs de p/q-points pour le volume incident à d
// Si le masque est mal défini (proche du bord) on ne garde que les p-points
// et le vecteur de q-points est vide après l'appel
template <typename MESH>
auto volumePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	// p-points : ajouter les sommets de Face(d) + Face(phi<2112>(m,d))
	Dart b = d;
	Dart t = b;
	do
	{
		p_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);
	b = phi<2112>(m, d);
	t = b;
	do
	{
		p_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);

	// q-points : pour chaque face f du volume,
	// ajouter Face(Inherit::phi<32112>(f)) (face opposée du volume adjacent)
	// si le volume adjacent n'existe pas (bord) alors le masque est mal défini,
	// donc on ignore les q-points
	Dart f = d;
	unsigned i = 0;
	do
	{
		// S'il n'y a pas de q-points pour cette face,
		// on vide le vecteur de q-points et on quitte la recherche de q-points
		if (is_boundary(m, phi3(m, f)))
		{
			q_point.clear();
			break;
		}
		// Sinon, récupération des brins de la face opposée du volume adjacent
		b = phi<32112>(m, f);
		t = b;
		do
		{
			q_point.push_back(t);
			t = phi1(m, t);
		} while (t != b);
		// face suivante sur l'hexaèdre
		f = phi2(m, (i % 2 == 0 ? phi1(m, f) : phi_1(m, f)));
		++i;
	} while (f != d);
}

// Remplit les vecteurs de p/q/r/s/t-points pour la face incidente à d
// Si le masque est mal défini (proche du bord) on ne garde que les p-points
// et les autres vecteurs sont vides après l'appel
template <typename MESH>
auto facePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point,
				   std::vector<Dart>& r_point, std::vector<Dart>& s_point, std::vector<Dart>& t_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	// p-points : ajouter les sommets de Face(d)
	Dart t = d;
	foreach_incident_vertex(m, typename MESH::Face(d), [&](typename MESH::Vertex dd) -> bool {
		p_point.push_back(dd.dart);
		return true;
	});

	// q-points : ajouter les sommets de Face(Inherit::phi<2112>(d)) + Face(Inherit::phi<32112>(d))
	Dart b = phi<2112>(m, d);
	t = b;
	do
	{
		q_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);
	b = phi<32112>(m, d);
	t = b;
	do
	{
		q_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);

	// r/t-points
	// Question : si un volume adjacent au volume sur ou sous la face n'existe pas, doit-on chercher
	// les r/s/t-points par d'autres chemins ? (nombre de config énorme...)

	// Pour chaque brin de p-point
	for (Dart p : p_point)
	{
		// Est-ce que les volumes adjacents aux volumes au dessus et au dessous existent ou est-ce en surface
		bool face_phi2_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi2(m, p))),
			 face_phi32_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi<32>(m, p)));

		Dart phi323211 = phi<323211>(m, p), phi3232111 = phi<3232111>(m, p), phi23211 = phi<23211>(m, p),
			 phi232111 = phi<232111>(m, p);

		// Si Face(phi<2>(p)) est sur la surface et Face(phi<32>(p)) aussi
		// alors il n'y a pas de r/t-points
		if (face_phi2_is_surface && face_phi32_is_surface)
		{
			q_point.clear();
			r_point.clear();
			s_point.clear();
			t_point.clear();
			break;
		}

		// Si Face(phi<2>(p)) est sur la surface et pas Face(phi<32>(p))
		// alors il y a deux r-points : Edge(phi<323211>(p))
		if (face_phi2_is_surface && !face_phi32_is_surface)
		{
			r_point.push_back(phi323211);
			r_point.push_back(phi3232111);
			continue;
		}

		// Si Face(phi<2>(p)) n'est pas sur la surface mais que Face(phi<32>(p)) l'est
		// alors il y a deux r-points : Edge(phi<23211>(p))
		if (!face_phi2_is_surface && face_phi32_is_surface)
		{
			r_point.push_back(phi23211);
			r_point.push_back(phi232111);
			continue;
		}

		// Ni Face(phi<2>(p)) ni Face(phi<32>(p)) ne sont sur la surface
		// alors soit Edge(phi<23211>(p)) == Edge(phi<323211>(p)) (topologie régulière)
		// soit Edge(phi<23211>(p)) != Edge(phi<323211>(p)) (arete extraordinaire)

		// Si Edge(phi<23211>(p)) != Edge(phi<323211>(p)) (arete extraordinaire)
		// alors il y a 4 t-points : Edge(phi<23211>(p)) + Edge(phi<323211>(p))
		if (index_of(m, typename MESH::Edge(phi23211)) != index_of(m, typename MESH::Edge(phi323211)))
		{
			t_point.push_back(phi23211);
			t_point.push_back(phi232111);
			t_point.push_back(phi323211);
			t_point.push_back(phi3232111);
		}
		// Sinon Edge(phi<23211>(p)) == Edge(phi<323211>(p))
		// et donc il y a 2 r-points
		else
		{
			r_point.push_back(phi23211);
			r_point.push_back(phi232111);
		}
	}

	// s-points : pour chaque q, ajouter Edge(phi<23211>(q))
	for (Dart q : q_point)
	{
		// S'il n'y a pas de s-points
		if (is_boundary(m, phi<23>(m, q)))
		{
			q_point.clear();
			r_point.clear();
			s_point.clear();
			t_point.clear();
			break;
		}
		s_point.push_back(phi<23211>(m, q));
		s_point.push_back(phi<232111>(m, q));
	}
}

// Remplit les vecteurs de p/q/r/s-points pour l'arete incidente à d
// Si le masque est mal défini (proche du bord) on ne garde que les p-points
// et les autres vecteurs sont vides après l'appel
// Attention, un appel récursif est fait (un seul)
template <typename MESH>
auto edgePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point,
				   std::vector<Dart>& r_point, std::vector<Dart>& s_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	if (is_incident_to_boundary(m, typename MESH::Edge(d)))
	{
		std::cerr << "Fatal error : edgePointMask called on dart which is surface" << std::endl;
		exit(EXIT_FAILURE);
	}

	// pour savoir si c'est le premier appel ou le rappel pour la deuxieme partie des r-points
	bool first_call = p_point.size() == 0;

	// p-points : ajouter les extremités de l'arete (sauf si c'est le deuxième appel
	if (first_call)
	{
		p_point.push_back(d);
		p_point.push_back(phi1(m, d));
	}

	// Pour toutes les faces incidentes à l'arete (au moins 4)
	// (arret garanti car l'arete n'est pas en surface)
	Dart t = d;
	do
	{
		// q-points : ajouter l'extremité de l'arete opposée (autre extremité : appel récursif)
		q_point.push_back(phi<11>(m, t));

		// r/s-points du coté sens direct
		// Question : s'il n'y a pas de volume adjacent à celui de la face courante,
		// doit-on chercher le point r/s par d'autres chemins (nombre de config!)

		bool face_phi12_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi<12>(m, t))),
			 face_phi132_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi<132>(m, t)));

		Dart phi1323211 = phi<1323211>(m, t), phi1232_1 = phi_1(m, phi<1232>(m, t));

		// Si Face(phi<12>(t)) et Face(phi<132>(t)) sont sur la surface
		// alors pas de r/s-points coté direct et on peux donc arreter la recherche le points
		if (face_phi12_is_surface && face_phi132_is_surface)
		{
			q_point.clear();
			r_point.clear();
			s_point.clear();
			break;
		}
		// Au moins l'une des face n'est pas sur la surface
		// Si ni Face(phi<12>(t)) ni Face(phi<132>(t)) ne sont sur la surface
		// alors soit Vertex(phi<1323211>(t)) == Vertex(phi<1232-1>(t)) (topologie régulière)
		// soit Vertex(phi<1323211>(t)) != Vertex(phi<1232-1>(t)) (arete extraordinaire)
		if (!face_phi12_is_surface && !face_phi132_is_surface)
		{
			// Si Vertex(phi<1323211>(t)) != Vertex(phi<1232-1>(t)) (arete extraordinaire)
			// alors 2 s-points Vertex(phi<1323211>(t)) et Vertex(phi<1232-1>(t))
			if (index_of(m, typename MESH::Vertex(phi1323211)) != index_of(m, typename MESH::Vertex(phi1232_1)))
			{
				s_point.push_back(phi1323211);
				s_point.push_back(phi1232_1);
			}
			// Sinon Vertex(phi<1323211>(t)) == Vertex(phi<1232-1>(t)) (topologie régulière)
			// alors 1 r-point
			else
			{
				r_point.push_back(phi1323211);
			}
		}
		// sinon si Face(phi<12>(t)) est sur la surface (et pas Face(phi<132>(t)))
		// alors 1 r-point : Vertex(phi<1323211>(t))
		else if (face_phi12_is_surface)
		{
			r_point.push_back(phi1323211);
		}
		// sinon (si Face(phi<132>(t)) est sur la surface et pas Face(phi<12>(t))))
		// alors 1 r-point : Vertex(phi<1232-1>(t))
		else
		{
			r_point.push_back(phi1232_1);
		}

		// Face incidente suivante
		t = phi<23>(m, t);
	} while (t != d);

	// Sinon, faire la meme chose dans l'autre sens (si ce n'est pas déjà fait)
	// pour avoir les r/s-points de l'autre coté de l'arete
	if (q_point.size() != 0 && first_call)
	{
		edgePointMask(m, phi2(m, d), p_point, q_point, r_point, s_point);
	}
}

/* ------------------------- BUTTERFLY SURFACE MASKS ------------------------- */

template <typename MESH>
auto surfaceAdjacentFace(const MESH& m, Dart d) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Dart>
{
	if (!is_boundary(m, phi3(m, d)))
	{
		std::cerr << "Fatal error : surfaceAdjacentFace called on dart which is not surface" << std::endl;
		exit(EXIT_FAILURE);
	}
	Dart t = phi2(m, d);
	while (!is_boundary(m, phi3(m, t)))
	{
		t = phi<32>(m, t);
	}
	return t;
}

template <typename MESH>
auto surfaceFacePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	Dart t = d;
	do
	{
		p_point.push_back(t);
		Dart e = surfaceAdjacentFace(m, t);
		q_point.push_back(phi<11>(m, e));
		q_point.push_back(phi_1(m, e));
		t = phi1(m, t);
	} while (t != d);
}

template <typename MESH>
auto surfaceEdgePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point,
						  std::vector<Dart>& r_point) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	p_point.push_back(d);
	p_point.push_back(phi1(m, d));

	q_point.push_back(phi<11>(m, d));
	q_point.push_back(phi_1(m, d));
	Dart t = surfaceAdjacentFace(m, phi1(m, d));
	r_point.push_back(phi_1(m, t));
	t = surfaceAdjacentFace(m, phi_1(m, d));
	r_point.push_back(phi<11>(m, t));

	Dart e = surfaceAdjacentFace(m, d);
	q_point.push_back(phi<11>(m, e));
	q_point.push_back(phi_1(m, e));
	t = surfaceAdjacentFace(m, phi1(m, e));
	r_point.push_back(phi_1(m, t));
	t = surfaceAdjacentFace(m, phi_1(m, e));
	r_point.push_back(phi<11>(m, t));
}

/* ----------------------------- BUTTERFLY RULES ----------------------------- */

template <typename MESH, typename T>
auto sum(const MESH& m, const std::vector<Dart>& v, typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	Vec3 s = Vec3(0.f, 0.f, 0.f);
	for (Dart i : v)
		s += value<T>(m, attribute, typename MESH::Vertex(i));
	return s;
}

template <typename T, typename MESH>
auto volumePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
					 typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	if (q_point.size() == 0)
	{
		return (1.f / 8) * sum<MESH>(m, p_point, attribute);
	}
	else
	{
		return ((6 * _W_ + 1) / 8) * sum<MESH>(m, p_point, attribute) - (_W_ / 4) * sum<MESH>(m, q_point, attribute);
	}
}

template <typename T, typename MESH>
auto facePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
				   const std::vector<Dart>& r_point, const std::vector<Dart>& s_point, const std::vector<Dart>& t_point,
				   typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	if (q_point.size() == 0)
	{
		return (1.f / 4) * sum<MESH>(m, p_point, attribute);
	}
	else
	{
		float w1 = (2 * _W_ + 1) / 4, w2 = _W_ / 4, w3 = _W_ / 8;
		return w1 * sum<MESH>(m, p_point, attribute) + w2 * sum<MESH>(m, q_point, attribute) -
			   w2 * sum<MESH>(m, r_point, attribute) - w3 * sum<MESH>(m, s_point, attribute) -
			   w3 * sum<MESH>(m, t_point, attribute);
	}
}

template <typename T, typename MESH>
auto edgePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
				   const std::vector<Dart>& r_point, const std::vector<Dart>& s_point,
				   typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	int32 N = uint32(q_point.size()) / 2;
	float32 w1 = 1.f / 2;
	if (N == 0)
	{
		return w1 * sum<MESH>(m, p_point, attribute);
	}
	else
	{
		float w2 = _W_ / N, w3 = _W_ / (2 * N);
		return w1 * sum<MESH>(m, p_point, attribute) + w2 * sum<MESH>(m, q_point, attribute) -
			   w2 * sum<MESH>(m, r_point, attribute) - w3 * sum<MESH>(m, s_point, attribute);
	}
}

template <typename T, typename MESH>
auto surfaceFacePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
						  typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	Vec3 sum_p = sum<MESH>(m, p_point, attribute), sum_q = sum<MESH>(m, q_point, attribute);
	float wp = 1.f / 4.f + _W_, wq = -_W_ / 2.f;
	return wp * sum_p + wq * sum_q;
}

template <typename T, typename MESH>
auto surfaceEdgePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
						  const std::vector<Dart>& r_point,
						  typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	Vec3 sum_p = sum<MESH>(m, p_point, attribute), sum_q = sum<MESH>(m, q_point, attribute),
		 sum_r = sum<MESH>(m, r_point, attribute);
	float wp = 1.f / 2.f, wq = _W_ / 2.f, wr = -_W_ / 2.f;
	return wp * sum_p + wq * sum_q + wr * sum_r;
}

template <typename MESH, typename FUNC>
auto subdivideEdge(MESH& m, Dart d, const FUNC& after_cut) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	after_cut(cut_edge(m, typename MESH::Edge(d)));
}

// /!\ The edges of the face are already cut
template <typename MESH, typename FUNC>
auto subdivideFace(MESH& m, Dart d, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Dart>
{
	using Vertex = typename MESH::Vertex;
	Dart e = phi1(m, d);
	Dart f = phi<1111>(m, e);
	cut_face(m, Vertex(e), Vertex(f));
	subdivideEdge(m, phi_1(m, e), after_cut);

	Dart result = phi_1(m, e);

	cut_face(m, Vertex(phi_1(m, e)), Vertex(phi<11>(m, e)));
	cut_face(m, Vertex(phi_1(m, f)), Vertex(phi<11>(m, f)));

	return result;
}

// /!\ The faces of the volume are already cut
template <typename MESH, typename FUNC>
auto subdivideVolume(MESH& m, Dart d, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Vertex = typename MESH::Vertex;
	Dart first_cut_dir = phi1(m, d);
	Dart left_cut_dir = phi<1>(m, first_cut_dir), right_cut_dir = phi<21112>(m, first_cut_dir);

	Dart cut_direction[4];
	cut_direction[0] = phi<1211>(m, left_cut_dir);
	cut_direction[1] = phi<1212111>(m, left_cut_dir);
	cut_direction[2] = phi<1211>(m, right_cut_dir);
	cut_direction[3] = phi<1212111>(m, right_cut_dir);

	auto cutVolume = [&](Dart d) -> typename MESH::Face {
		std::vector<Dart> cut_path;
		Dart t = d;
		do
		{
			cut_path.push_back(t);
			t = phi<121>(m, t);
		} while (t != d);
		return cut_volume(m, cut_path);
	};

	typename MESH::Face F = cutVolume(first_cut_dir);
	subdivideFace(m, F.dart, after_cut);

	F = cutVolume(left_cut_dir);
	cut_face(m, Vertex(F.dart), Vertex(phi<111>(m, F.dart)));
	F = cutVolume(right_cut_dir);
	cut_face(m, Vertex(F.dart), Vertex(phi<111>(m, F.dart)));
	for (int i = 0; i < 4; ++i)
	{
		cutVolume(cut_direction[i]);
	}
}

template <typename MESH, typename FUNC>
auto subdivideListEdges(MESH& m, std::vector<Dart>& edges, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&> && !(std::is_convertible_v<MESH&, CPH3&>)>
{
	for (Dart d : edges)
	{
		subdivideEdge(m, d, after_cut);
	}
}

template <typename MR_MESH, typename FUNC>
auto subdivideListEdges(MR_MESH& m, std::vector<Dart>& edges, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MR_MESH&, CPH3&>>
{
	MR_MESH m2(m);
	for (Dart d : edges)
	{
		m2.current_level_ = m.edge_level(d) + 1;
		subdivideEdge(m2, d, after_cut);
	}
}

template <typename MESH, typename FUNC>
auto subdivideListFaces(MESH& m, std::vector<Dart>& faces, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&> && !(std::is_convertible_v<MESH&, CPH3&>)>
{
	for (Dart d : faces)
	{
		subdivideFace(m, d, after_cut);
	}
}

template <typename MR_MESH, typename FUNC>
auto subdivideListFaces(MR_MESH& m, std::vector<Dart>& faces, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MR_MESH&, CPH3&>>
{
	MR_MESH m2(m);
	for (Dart d : faces)
	{
		m2.current_level_ = m.face_level(d) + 1;
		subdivideFace(m2, d, after_cut);
	}
}

template <typename MESH, typename FUNC>
auto subdivideListVolumes(MESH& m, std::vector<Dart>& volumes, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&> && !(std::is_convertible_v<MESH&, CPH3&>)>
{
	for (Dart d : volumes)
	{
		subdivideVolume(m, d, after_cut);
	}
}

template <typename MR_MESH, typename FUNC>
auto subdivideListVolumes(MR_MESH& m, std::vector<Dart>& volumes, const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MR_MESH&, CPH3&>>
{
	MR_MESH m2(m);
	for (Dart d : volumes)
	{
		m2.current_level_ = m.volume_level(d) + 1;
		subdivideVolume(m2, d, after_cut);
	}
}

template <typename MESH>
auto isHexaVolume(MESH& m, Dart d) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, bool>
{
	using Face = typename MESH::Face;
	using Volume = typename MESH::Volume;
	bool r = true;
	unsigned v = 0;
	foreach_incident_face(m, Volume(d), [&](Face f) -> bool {
		if (!(codegree(m, f) == 4))
		{
			r = false;
		}
		++v;
		return r && v <= 6;
	});
	return (r ? (v == 6) : false);
}

template <typename MESH>
auto isValidForSubdivision(MESH& m, Dart d) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, bool>
{
	using Face = typename MESH::Face;
	using Volume = typename MESH::Volume;

	CellMarker<MESH, Face> fm(m);
	bool v = false;
	if (isHexaVolume(m, d))
	{
		v = true;
		foreach_incident_face(m, Volume(d), [&](Face f) -> bool {
			if (!is_boundary(m, phi3(m, f.dart)) && !isHexaVolume(m, phi3(m, f.dart)))
			{
				v = false;
			}
			return v;
		});
	}
	return v;
}

template <typename MESH>
auto butterflySubdivisionVolumeAdaptative(MESH& m, double angle_threshold,
										  std::vector<typename mesh_traits<MESH>::template Attribute<Vec3>*> attributes)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	butterflySubdivisionVolumeAdaptative(m, angle_threshold, attributes, [&](typename MESH::Vertex) {});
}

template <typename MESH, typename FUNC>
auto butterflySubdivisionVolumeAdaptative(MESH& m, double angle_threshold,
										  std::vector<typename mesh_traits<MESH>::template Attribute<Vec3>*> attributes,
										  const FUNC& after_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	butterflySubdivisionVolumeAdaptative(m, angle_threshold, attributes, after_cut, after_cut, after_cut);
}

template <typename MESH, typename FUNC>
auto butterflySubdivisionVolumeAdaptative(MESH& m, double angle_threshold,
										  std::vector<typename mesh_traits<MESH>::template Attribute<Vec3>*> attributes,
										  const FUNC& after_cut_edge, const FUNC& after_cut_face,
										  const FUNC& after_cut_volume)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Volume = typename MESH::Volume;
	using Face = typename MESH::Face;
	using Face2 = typename MESH::Face2;
	using Edge = typename MESH::Edge;
	using Vertex = typename MESH::Vertex;

	CellMarker<MESH, Edge> cm_surface(m);
	CellMarker<MESH, Volume> cm_cell(m);

	std::vector<typename mesh_traits<MESH>::template Attribute<Vec3>*> attrs;
	for (auto a : attributes)
		if (a)
			attrs.push_back(a);

	for (Dart d = m.begin(); d != m.end(); d = m.next(d))
	{
		if (is_boundary(m, phi3(m, d)) && !cm_surface.is_marked(Edge(d)))
		{
			Dart d2 = phi2(m, d);
			while (!is_boundary(m, phi3(m, d2)))
			{
				d2 = phi<32>(m, d2);
			}
			// First attribute is the one watch for adaptive subdivision
			auto edge_angle = geometry::angle(m, Face2(d), Face2(d2), attrs.front());
			if (std::abs(edge_angle) > angle_threshold)
			{
				if (!cm_cell.is_marked(Volume(d)) && isValidForSubdivision(m, d))
					cm_cell.mark(Volume(d));
				if (!cm_cell.is_marked(Volume(d2)) && isValidForSubdivision(m, d2))
					cm_cell.mark(Volume(d2));
			}
			cm_surface.mark(Edge(d));
		}
	}

	CellMarker<MESH, Edge> cm_edge(m);
	CellMarker<MESH, Face> cm_face(m);
	CellMarker<MESH, Volume> cm_volume(m);
	std::queue<std::queue<Vec3>> volume_points, face_points, edge_points;
	std::vector<Dart> edges, faces, volumes;
	std::vector<Dart> p_point, q_point, r_point, s_point, t_point;

	// computing the new vertices's embedding
	for (Dart d = m.begin(); d != m.end(); d = m.next(d))
	{
		if (!is_boundary(m, d) && cm_cell.is_marked(Volume(d)))
		{
			foreach_dart_of_orbit(m, Volume(d), [&](Dart t) -> bool {
				// edges vertices
				if (!cm_edge.is_marked(Edge(t)))
				{
					if (is_incident_to_boundary(m, Edge(t)) && !is_boundary(m, phi3(m, t)))
					{
						// we ignore the surfaces darts which are in junction of two volumes
					}
					// surface case
					else if (is_boundary(m, phi3(m, t)))
					{
						p_point.clear();
						q_point.clear();
						r_point.clear();
						surfaceEdgePointMask(m, t, p_point, q_point, r_point);
						std::queue<Vec3> list_points;
						for (auto a : attrs)
							list_points.push(surfaceEdgePointRule<Vec3>(m, p_point, q_point, r_point, a));

						edge_points.push(list_points);
						edges.push_back(t);
						cm_edge.mark(Edge(t));
					}
					// volume case
					else
					{
						p_point.clear();
						q_point.clear();
						r_point.clear();
						s_point.clear();
						edgePointMask(m, t, p_point, q_point, r_point, s_point);

						std::queue<Vec3> list_points;
						for (auto a : attrs)
							list_points.push(edgePointRule<Vec3>(m, p_point, q_point, r_point, s_point, a));

						edge_points.push(list_points);

						edges.push_back(t);
						cm_edge.mark(Edge(t));
					}
				}
				// faces vertices
				if (!cm_face.is_marked(Face(t)))
				{
					// cas surface
					if (is_incident_to_boundary(m, Face(t)))
					{
						p_point.clear();
						q_point.clear();
						surfaceFacePointMask(m, t, p_point, q_point);

						std::queue<Vec3> list_points;
						for (auto a : attrs)
							list_points.push(surfaceFacePointRule<Vec3>(m, p_point, q_point, a));

						face_points.push(list_points);
					}
					// volume case
					else
					{
						p_point.clear();
						q_point.clear();
						r_point.clear();
						s_point.clear();
						t_point.clear();
						facePointMask(m, t, p_point, q_point, r_point, s_point, t_point);

						std::queue<Vec3> list_points;
						for (auto a : attrs)
							list_points.push(facePointRule<Vec3>(m, p_point, q_point, r_point, s_point, t_point, a));

						face_points.push(list_points);
					}
					faces.push_back(t);
					cm_face.mark(Face(t));
				}
				// volumes vertices
				if (!cm_volume.is_marked(Volume(t)))
				{
					p_point.clear();
					q_point.clear();
					volumePointMask(m, t, p_point, q_point);

					std::queue<Vec3> list_points;
					for (auto a : attrs)
						list_points.push(volumePointRule<Vec3>(m, p_point, q_point, a));

					volume_points.push(list_points);
					volumes.push_back(t);
					cm_volume.mark(Volume(t));
				}
				return true;
			});
			cm_cell.unmark(Volume(d));
		}
	}

	subdivideListEdges<MESH>(m, edges, [&](Vertex v) {
		for (auto a : attrs)
		{
			value<Vec3>(m, a, v) = edge_points.front().front();
			edge_points.front().pop();
		}
		edge_points.pop();
		after_cut_edge(v);
	});

	subdivideListFaces(m, faces, [&](Vertex v) {
		for (auto a : attrs)
		{
			value<Vec3>(m, a, v) = face_points.front().front();
			face_points.front().pop();
		}
		face_points.pop();
		after_cut_face(v);
	});

	subdivideListVolumes(m, volumes, [&](Vertex v) {
		for (auto a : attrs)
		{
			value<Vec3>(m, a, v) = volume_points.front().front();
			volume_points.front().pop();
		}
		volume_points.pop();
		after_cut_volume(v);
	});
}

template <typename MRMESH>
auto butterflyMultiresolution(
	MRMESH& m, double angle_threshold, std::vector<typename mesh_traits<MRMESH>::template Attribute<Vec3>*> attributes,
	typename mesh_traits<MRMESH>::template Attribute<std::array<typename mesh_traits<MRMESH>::Vertex, 3>>* parents =
		nullptr,
	typename mesh_traits<MRMESH>::template Attribute<Vec3>* position_relative = nullptr)
{
	using Vertex = typename mesh_traits<MRMESH>::Vertex;
	using Volume = typename mesh_traits<MRMESH>::Volume;
	std::vector<Vertex> new_vertices;

	if (!parents || !position_relative)
	{
		butterflySubdivisionVolumeAdaptative(m, angle_threshold, attributes);
		return;
	}

	MRMESH m2(m);

	m2.current_level_ = m.maximum_level_ + 1;

	auto after_cut = [&](Vertex v) {
		Dart old = m2.volume_oldest_dart(v.dart);
		m2.current_level_--;
		int i = 0;
		std::array<Vertex, 3> p;
		foreach_incident_vertex(m2, Volume(old), [&](Vertex w) -> bool {
			p[i++] = w;
			return i < 3;
		});
		value<std::array<Vertex, 3>>(m2, parents, v) = p;
		new_vertices.push_back(v);
		m2.current_level_++;
	};

	butterflySubdivisionVolumeAdaptative(m, angle_threshold, attributes, after_cut);
	m2.current_level_--;
	for (auto v : new_vertices)
	{
		auto pos = attributes[0];
		std::array<Vertex, 3> p = value<std::array<Vertex, 3>>(m, parents, v);
		Vec3 A = value<Vec3>(m, pos, p[0]);
		Vec3 B = value<Vec3>(m, pos, p[1]);
		Vec3 C = value<Vec3>(m, pos, p[2]);

		Vec3 X = (B - A).normalized();
		Vec3 Y = (C - A).normalized();
		Vec3 Z = (X.cross(Y));
		Z = Z.normalized();
		Eigen::Matrix3d mb;
		mb.col(0) = X;
		mb.col(1) = Y;
		mb.col(2) = Z;
		Eigen::Matrix3d inv;
		bool is_inversible;
		mb.computeInverseWithCheck(inv, is_inversible);
		assert(is_inversible);
		value<Vec3>(m, position_relative, v) = inv * (value<Vec3>(m, pos, v) - A);
	}
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_H_
