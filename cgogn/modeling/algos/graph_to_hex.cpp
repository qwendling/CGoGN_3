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

#include <cgogn/modeling/algos/graph_to_hex.h>

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/halfedge.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>
#include <cgogn/io/surface/surface_import.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/intersection.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/modeling/algos/subdivision.h>

#include <iostream>
#include <fstream>

namespace cgogn
{

namespace modeling
{

///////////
// CMaps //
///////////

bool graph_to_hex(Graph& g, CMap2& m2, CMap3& m3)
{
	bool okay;

	GData gData;
	GAttributes gAttribs;
	M2Attributes m2Attribs;
	M3Attributes m3Attribs;
	//okay = subdivide_graph(g);
	//if (!okay)
	//{
	//	std::cout << "error graph_to_hex: subdivide_graph" << std::endl;
	//	return false;
	//}
	//else
	//	std::cout << "graph_to_hex (/): sudivided graph" << std::endl;

	okay = get_graph_data(g, gData);
	std::cout << uint32(gData.intersections.size()) << " inters" << std::endl;
	std::cout << uint32(gData.branches.size()) << " branches" << std::endl;
	for (auto b : gData.branches)
		std::cout << index_of(g, Graph::Vertex(b.first.dart)) << " - " << index_of(g, Graph::Vertex(b.second.dart))
				  << std::endl;
	if (!okay)
	{
		std::cout << "error graph_to_hex: get_graph_data" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): got graph data" << std::endl;

	okay = add_graph_attributes(g, gAttribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: add_graph_attributes" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): added graph attributes" << std::endl;

	okay = add_cmap2_attributes(m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: add_cmap2_attributes" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): added cmap2 attributes" << std::endl;

	okay = build_contact_surfaces(g, gAttribs, m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: build_contact_surfaces" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): contact surfaces built" << std::endl;

	okay = create_intersection_frames(g, gAttribs, m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: create_intersections_frames" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): create_intersections_frames completed" << std::endl;

	okay = propagate_frames(g, gAttribs, gData, m2);
	if (!okay)
	{
		std::cout << "error graph_to_hex: propagate_frames" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): propagate_frames completed" << std::endl;

	okay = set_contact_surfaces_geometry(g, gAttribs, m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: set_contact_surfaces_geometry" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): set_contact_surfaces_geometry completed" << std::endl;

	okay = build_branch_sections(g, gAttribs, m2, m2Attribs, m3);
	if (!okay)
	{
		std::cout << "error graph_to_hex: build_branch_sections" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): build_branch_sections completed" << std::endl;

	okay = sew_branch_sections(m2, m2Attribs, m3);
	if (!okay)
	{
		std::cout << "error graph_to_hex: sew_sections" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): sew_sections completed" << std::endl;

	okay = set_volumes_geometry(m2, m2Attribs, m3);
	if (!okay)
	{
		std::cout << "error graph_to_hex: set_volumes_geometry" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): set_volumes_geometry completed" << std::endl;
	

	okay = add_quality_attributes(m3, m3Attribs);
	okay = set_hex_frames(m3, m3Attribs);
	okay = compute_scaled_jacobians(m3, m3Attribs, true);
	okay = compute_jacobians(m3, m3Attribs, true);
	okay = compute_maximum_aspect_frobenius(m3, m3Attribs, true);
	okay = compute_mean_aspect_frobenius(m3, m3Attribs, true);

	return okay;
}

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

void index_volume_cells(CMap2& m, CMap2::Volume vol)
{
	DartMarkerStore<CMap2> vertex_marker(m);
	DartMarkerStore<CMap2> edge_marker(m);
	DartMarkerStore<CMap2> face_marker(m);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (is_indexed<CMap2::Vertex>(m))
		{
			if (!vertex_marker.is_marked(d))
			{
				CMap2::Vertex v(d);
				foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
					vertex_marker.mark(d);
					return true;
				});
				if (index_of<CMap2::Vertex>(m, v) == INVALID_INDEX)
					set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			}
		}
		if (is_indexed<CMap2::Edge>(m) || is_indexed<CMap2::HalfEdge>(m))
		{
			if (!edge_marker.is_marked(d))
			{
				CMap2::Edge e(d);
				foreach_dart_of_orbit(m, e, [&](Dart d) -> bool {
					edge_marker.mark(d);
					return true;
				});
				if (is_indexed<CMap2::Edge>(m))
					if (index_of<CMap2::Edge>(m, e) == INVALID_INDEX)
						set_index(m, e, new_index<CMap2::Edge>(m));
				if (is_indexed<CMap2::HalfEdge>(m))
				{
					if (index_of<CMap2::HalfEdge>(m, CMap2::HalfEdge(e.dart)) == INVALID_INDEX)
						set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
					if (index_of<CMap2::HalfEdge>(m, CMap2::HalfEdge(phi2(m, e.dart))) == INVALID_INDEX)
						set_index(m, CMap2::HalfEdge(phi2(m, e.dart)), new_index<CMap2::HalfEdge>(m));
				}
				return true;
			}
		}
		if (is_indexed<CMap2::Face>(m))
		{
			if (!face_marker.is_marked(d))
			{
				CMap2::Face f(d);
				foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
					face_marker.mark(d);
					return true;
				});
				if (index_of<CMap2::Face>(m, f) == INVALID_INDEX)
					set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			}
		}
		return true;
	});
	if (is_indexed<CMap2::Volume>(m))
		if (index_of<CMap2::Volume>(m, vol) == INVALID_INDEX)
			set_index(m, vol, new_index<CMap2::Volume>(m));
}

void sew_volumes(CMap3& m, Dart d0, Dart d1)
{
	cgogn_message_assert(codegree(m, CMap3::Face(d0)) == codegree(m, CMap3::Face(d1)),
						 "The faces to sew do not have the same codegree");
	Dart it0 = d0;
	Dart it1 = d1;
	do
	{
		cgogn_message_assert(phi3(m, it0) == it0 && phi3(m, it1) == it1, "The faces to sew are already sewn");
		phi3_sew(m, it0, it1);
		it0 = phi1(m, it0);
		it1 = phi_1(m, it1);
	} while (it0 != d0);
}

Graph::HalfEdge branch_extremity(const Graph& g, Graph::HalfEdge v1, CellMarker<Graph, Graph::Edge>& cm)
{
	Dart d = v1.dart;
	while (degree(g, Graph::Vertex(d)) == 2)
	{
		d = alpha0(g, alpha1(g, d));
		cm.mark(Graph::Edge(d));
	}
	return Graph::HalfEdge(d);
}

Dart add_branch_section(CMap3& m3)
{
	std::vector<Dart> D = {add_prism(static_cast<CMap2&>(m3), 4).dart, add_prism(static_cast<CMap2&>(m3), 4).dart,
						   add_prism(static_cast<CMap2&>(m3), 4).dart, add_prism(static_cast<CMap2&>(m3), 4).dart};

	sew_volumes(m3, phi2(m3, D[0]), phi2(m3, phi_1(m3, (D[1]))));
	sew_volumes(m3, phi2(m3, D[1]), phi2(m3, phi_1(m3, (D[2]))));
	sew_volumes(m3, phi2(m3, D[2]), phi2(m3, phi_1(m3, (D[3]))));
	sew_volumes(m3, phi2(m3, D[3]), phi2(m3, phi_1(m3, (D[0]))));

	return D[0];
}

void project_on_sphere(Vec3& P, const Vec3& C, Scalar R)
{
	P = C + (P - C).normalized() * R;
}

void shift_frame(Mat3& frame, uint32 nb_shifts)
{
	for (uint32 i = 0; i < nb_shifts; ++i)
	{
		Vec3 R = frame.col(1);
		Vec3 S = -frame.col(0);
		frame.col(0) = R;
		frame.col(1) = S;
	}
}

void dualize_volume(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs, const Graph& g, GAttributes& gAttribs)
{
	// set the new phi1
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		Dart dd = phi2(m, phi_1(m, d));
		(*(m.phi1_))[d.index] = dd;
		return true;
	});

	// set the new phi_1
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		Dart next = phi1(m, d);
		(*(m.phi_1_))[next.index] = d;
		return true;
	});

	DartMarkerStore<CMap2> face_marker(m);
	const Graph::Vertex gv = value<Graph::Vertex>(m, m2Attribs.volume_gvertex, vol);
	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, gv);
	const Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, gv);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (!face_marker.is_marked(d))
		{
			CMap2::Face f(d);
			foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
				face_marker.mark(d);
				return true;
			});
			if (is_indexed<CMap2::Face>(m))
				set_index(m, f, new_index<CMap2::Face>(m)); // give a new index to the face
			// darts of the new face orbit are darts of the old dual vertex
			// get the outgoing graph branch from the old dual vertex
			Dart branch = value<Dart>(m, m2Attribs.dual_vertex_graph_branch, CMap2::Vertex(d));
			// store the contact surface face on the outgoing branch
			value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(branch)) = d;
			return true;
		}
		return true;
	});

	DartMarkerStore<CMap2> vertex_marker(m);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (!vertex_marker.is_marked(d))
		{
			CMap2::Vertex v(d);

			std::vector<Vec3> points;
			 foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
				vertex_marker.mark(d);
				points.push_back(value<Vec3>(m, m2Attribs.vertex_position, CMap2::Vertex(d)) - center);
				return true;
			});
			Vec3 b;
			 if (points.size() == 2)
				 b = slerp(points[0], points[1], 0.5, true) + center;
			 else
				b = spherical_barycenter(points, 10) + center;

			set_index(m, v, new_index<CMap2::Vertex>(m));	  // give a new index to the vertex
			value<Vec3>(m, m2Attribs.vertex_position, v) = b; // set the position to the computed position
			return true;
		}
		return true;
	});
}

void extract_volume_surface(CMap3& m3, CMap2& m2)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position3 = get_attribute<Vec3, CMap3::Vertex>(m3, "position");
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position2 = add_attribute<Vec3, CMap2::Vertex>(m2, "position");

	auto vertex_m2_embeddings = cgogn::add_attribute<uint32, CMap3::Vertex>(m3, "vertex_m2_embeddings");
	
	cgogn::io::SurfaceImportData surface_data;

	foreach_cell(m3, [&](CMap3::Vertex v3) -> bool {
		if (is_incident_to_boundary(m3, v3))
		{
			uint32 vertex_id = new_index<CMap2::Vertex>(m2);
			value<uint32>(m3, vertex_m2_embeddings, v3) = vertex_id;
			Vec3 p = value<Vec3>(m3, vertex_position3, v3);
			(*vertex_position2)[vertex_id] = p;
			surface_data.vertices_id_.push_back(vertex_id);
		}
		return true;
	});

	std::vector<uint32> indices;
	indices.reserve(4);
	foreach_cell(m3, [&](CMap3::Face f) -> bool {
		if (is_incident_to_boundary(m3, f))
		{
			foreach_incident_vertex(m3, f, [&](CMap3::Vertex v3) -> bool {
				//uint32 vertex_id = value<uint32>(m3, vertex_m2_embeddings, v3);
				indices.push_back(value<uint32>(m3, vertex_m2_embeddings, v3));
				return true;
			});

			surface_data.faces_nb_vertices_.push_back(4);
			surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
													  indices.end());
			indices.clear();
		}
		return true;
	});

	import_surface_data(m2, surface_data);

	remove_attribute<CMap3::Vertex>(m3, vertex_m2_embeddings);
}

void catmull_clark_approx(CMap2& m, uint32 iterations)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position = get_attribute<Vec3, CMap2::Vertex>(m, "position");
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_delta = add_attribute<Vec3, CMap2::Vertex>(m, "delta");
	std::shared_ptr<CMap2::Attribute<Vec3>> incident_faces_mid = add_attribute<Vec3, CMap2::Vertex>(m, "incident_faces_mid");

	for (uint32 i = 0; i < iterations; ++i)
	{
		CellCache<CMap2> cache_init_vertices(m);
		cache_init_vertices.template build<CMap2::Vertex>();
		CellCache<CMap2> cache_edge_vertices(m);
		CellCache<CMap2> cache_face_vertices(m);

		cgogn::modeling::quadrangulate_all_faces(
			m,
			[&](CMap2::Vertex v) -> void {
				cache_edge_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				value<Vec3>(m, incident_faces_mid, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += 0.5 * value<Vec3>(m, vertex_position, v2);
					value<Vec3>(m, vertex_delta, v) -= 0.25 * value<Vec3>(m, vertex_position, v2);
					return true;
				});
			},
			[&](CMap2::Vertex v) -> void {
				cache_face_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_position, v2);
					return true;
				});
				value<Vec3>(m, vertex_position, v) /= degree(m, v);

				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_delta, v2) += 0.25 * value<Vec3>(m, vertex_position, v);
					value<Vec3>(m, incident_faces_mid, v2) += 0.5 * value<Vec3>(m, vertex_position, v);
					return true;
				});
			});

		foreach_cell(cache_init_vertices, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
			Vec3 F = {0, 0, 0};
			Vec3 R = {0, 0, 0};
			Scalar n = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				F += value<Vec3>(m, incident_faces_mid, v2);
				R += value<Vec3>(m, vertex_position, v2);
				n += 1;
				return true;
			});

			value<Vec3>(m, vertex_delta, v) = (-3 * n * value<Vec3>(m, vertex_position, v) + F + 2 * R) / (n * n);
			return true;
		});

		foreach_cell(m, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_delta, v);
			return true;
		});
	}
	remove_attribute<CMap2::Vertex>(m, vertex_delta);
	remove_attribute<CMap2::Vertex>(m, incident_faces_mid);
}

void catmull_clark_inter(CMap2& m, uint32 iterations)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position = get_attribute<Vec3, CMap2::Vertex>(m, "position");
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position2 = add_attribute<Vec3, CMap2::Vertex>(m, "position2");
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_delta = add_attribute<Vec3, CMap2::Vertex>(m, "delta");
	std::shared_ptr<CMap2::Attribute<Vec3>> incident_faces_mid =
		add_attribute<Vec3, CMap2::Vertex>(m, "incident_faces_mid");
	std::shared_ptr<CMap2::Attribute<Vec3>> incident_init_mid =
		add_attribute<Vec3, CMap2::Vertex>(m, "incident_init_mid");

	for (uint32 i = 0; i < iterations; ++i)
	{
		CellCache<CMap2> cache_init_vertices(m);
		cache_init_vertices.template build<CMap2::Vertex>();
		CellCache<CMap2> cache_edge_vertices(m);
		CellCache<CMap2> cache_face_vertices(m);

		foreach_cell(m, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
			return true;
		});

		cgogn::modeling::quadrangulate_all_faces(
			m,
			[&](CMap2::Vertex v) -> void {
				cache_edge_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				value<Vec3>(m, incident_faces_mid, v) = {0, 0, 0};
				value<Vec3>(m, incident_init_mid, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += 0.5 * value<Vec3>(m, vertex_position, v2);
					return true;
				});
			},
			[&](CMap2::Vertex v) -> void {
				cache_face_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_position, v2);
					return true;
				});
				value<Vec3>(m, vertex_position, v) /= degree(m, v);

				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_delta, v2) += 0.25 * value<Vec3>(m, vertex_position, v);
					value<Vec3>(m, incident_faces_mid, v2) += 0.5 * value<Vec3>(m, vertex_position, v);
					return true;
				});
			});

		foreach_cell(cache_init_vertices, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, incident_faces_mid, v) = {0, 0, 0};
			Scalar n = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				value<Vec3>(m, incident_faces_mid, v) += value<Vec3>(m, incident_faces_mid, v2);
				n += 1;
				return true;
			});
			value<Vec3>(m, incident_faces_mid, v) /= n;

			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				value<Vec3>(m, vertex_delta, v2) -= 0.25 * value<Vec3>(m, incident_faces_mid, v);
				value<Vec3>(m, incident_init_mid, v2) += 0.5 * value<Vec3>(m, incident_faces_mid, v);
				return true;
			});
			return true; 
		});

		foreach_cell(cache_face_vertices, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_delta, v);
			Scalar n = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				value<Vec3>(m, vertex_delta, v) -= 2 * value<Vec3>(m, incident_faces_mid, v2);
				value<Vec3>(m, vertex_delta, v) -= value<Vec3>(m, incident_init_mid, v2);
				n += 1;
				return true;
			});
			value<Vec3>(m, vertex_delta, v) += 3 * n * value<Vec3>(m, vertex_position, v);
			value<Vec3>(m, vertex_delta, v) /= n * n;
			return true;
		});

		foreach_cell(m, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_delta, v);
			return true;
		});
	}
	remove_attribute<CMap2::Vertex>(m, vertex_delta);
	remove_attribute<CMap2::Vertex>(m, incident_faces_mid);
	remove_attribute<CMap2::Vertex>(m, incident_init_mid);
}

bool intersection_surface(const CMap2& m, const CMap2::Attribute<Vec3>* vertex_position, const Vec3& P, const Vec3& Dir,
						  Vec3* inter)
{
	Scalar min_dist = std::numeric_limits<Scalar>::max();
	Vec3 closest_inter;
	bool found_point = false;

	foreach_cell(m, [&](CMap2::Face f) -> bool { 
		std::vector<CMap2::Vertex> vertices = incident_vertices(m, f);
		if (vertices.size() > 2)
		{
			for (uint32 i = 2; i < vertices.size(); ++i)
			{
				Vec3 inter_t;
				Vec3 a = value<Vec3>(m, vertex_position, vertices[0]);
				Vec3 b = value<Vec3>(m, vertex_position, vertices[i - 1]);
				Vec3 c = value<Vec3>(m, vertex_position, vertices[i]);
				bool hit = cgogn::geometry::intersection_ray_triangle(P, Dir, a, b, c, &inter_t);
				
				if (hit && (inter_t - P).norm() < min_dist)
				{
					found_point = true;
					min_dist = (inter_t - P).norm();
					closest_inter = inter_t;
					break;
				}
			}
		}

		return true; 
	});

	if (found_point)
		*inter = closest_inter;

	return found_point;
}


/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool subdivide_graph(Graph& g)
{
	auto vertex_position = get_attribute<Vec3, Graph::Vertex>(g, "position");
	if (!vertex_position)
	{
		std::cout << "The graph has no vertex position attribute" << std::endl;
		return false;
	}

	auto vertex_radius = get_attribute<Scalar, Graph::Vertex>(g, "radius");
	if (!vertex_radius)
	{
		std::cout << "The graph has no vertex radius attribute" << std::endl;
		return false;
	}

	CellCache<Graph> cache(g);
	cache.template build<Graph::Edge>();

	foreach_cell(cache, [&](Graph::Edge eg) -> bool {
		const Vec3& P0 = value<Vec3>(g, vertex_position, Graph::Vertex(eg.dart));
		const Vec3& Pn = value<Vec3>(g, vertex_position, Graph::Vertex(alpha0(g, eg.dart)));
		const Scalar R0 = value<Scalar>(g, vertex_radius, Graph::Vertex(eg.dart));
		const Scalar Rn = value<Scalar>(g, vertex_radius, Graph::Vertex(alpha0(g, eg.dart)));

		Scalar avg_radius = (R0 + Rn) / 2;
		Vec3 edge = Pn - P0;
		Scalar D = edge.norm();
		edge = edge.normalized();

		uint32 n = uint32(D / avg_radius);
		//if (D > 2 * avg_radius)
		//{
		//	n = 2;
			if (n > 1)
			{
				//Scalar y = (Rn - R0) / D;
				Scalar ratio = pow(Rn / R0, 1.0 / Scalar(n));

				Scalar sum_ratio = 0;
				Scalar alpha = 1;
				std::vector<Scalar> ratios_sums;
				for (uint32 i = 0; i < n; ++i)
				{
					sum_ratio += alpha;
					ratios_sums.push_back(sum_ratio);
					alpha *= ratio;
				}
				ratios_sums.pop_back();

				Vec3 D0 = D / sum_ratio * edge;
				std::vector<Vec3> Pi;
				for (Scalar r : ratios_sums)
				{
					Pi.push_back(P0 + r * D0);
				}

				Dart d = eg.dart;
				alpha = ratio * R0;
				for (Vec3 P : Pi)
				{
					d = cut_edge(g, Graph::Edge(d), true).dart;
					value<Vec3>(g, vertex_position, Graph::Vertex(d)) = P;
					value<Scalar>(g, vertex_radius, Graph::Vertex(d)) = alpha;
					alpha *= ratio;
					d = alpha1(g, d);
				}
			}
		//}
		return true;
	});
	return true;
}

bool get_graph_data(const Graph& g, GData& gData)
{
	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		if (degree(g, v) > 2)
			gData.intersections.push_back(v);
		return true;
	});

	CellMarker<Graph, Graph::Edge> cm(g);
	foreach_cell(g, [&](Graph::Edge e) -> bool {
		if (cm.is_marked(e))
			return true;
		cm.mark(e);
		std::vector<Graph::HalfEdge> halfedges = incident_halfedges(g, e);
		gData.branches.push_back({branch_extremity(g, halfedges[0], cm), branch_extremity(g, halfedges[1], cm)});
		return true;
	});

	return true;
}

bool add_graph_attributes(Graph& g, GAttributes& gAttribs)
{
	gAttribs.vertex_position = get_attribute<Vec3, Graph::Vertex>(g, "position");
	if (!gAttribs.vertex_position)
	{
		std::cout << "The graph has no vertex position attribute" << std::endl;
		return false;
	}

	gAttribs.vertex_radius = get_attribute<Scalar, Graph::Vertex>(g, "radius");
	if (!gAttribs.vertex_radius)
	{
		std::cout << "The graph has no vertex radius attribute" << std::endl;
		return false;
	}

	gAttribs.vertex_contact_surface = add_attribute<Dart, Graph::Vertex>(g, "contact_surface");
	if (!gAttribs.vertex_contact_surface)
	{
		std::cout << "Failed to add vertex_contact_surface attribute to graph" << std::endl;
		return false;
	}

	gAttribs.halfedge_contact_surface_face = add_attribute<Dart, Graph::HalfEdge>(g, "contact_surface_face");
	if (!gAttribs.halfedge_contact_surface_face)
	{
		std::cout << "Failed to add halfedge_contact_surface_face attribute to graph" << std::endl;
		return false;
	}

	gAttribs.halfedge_frame = add_attribute<Mat3, Graph::HalfEdge>(g, "frame");
	if (!gAttribs.halfedge_frame)
	{
		std::cout << "Failed to add halfedge_frame attribute to graph" << std::endl;
		return false;
	}

	return true;
}

bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs)
{
	m2Attribs.vertex_position = add_attribute<Vec3, CMap2::Vertex>(m2, "position");
	if (!m2Attribs.vertex_position)
	{
		std::cout << "Failed to add vertex_position attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.dual_vertex_graph_branch = add_attribute<Dart, CMap2::Vertex>(m2, "graph_branch");
	if (!m2Attribs.dual_vertex_graph_branch)
	{ 
		std::cout << "Failed to add dual_vertex_graph_branch attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.volume_center = add_attribute<Vec3, CMap2::Volume>(m2, "center");
	if (!m2Attribs.volume_center)
	{
		std::cout << "Failed to add volume_center attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.volume_gvertex = add_attribute<Graph::Vertex, CMap2::Volume>(m2, "gvertex");
	if (!m2Attribs.volume_gvertex)
	{
		std::cout << "Failed to add volume_gvertex attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.edge_mid = add_attribute<Vec3, CMap2::Edge>(m2, "edge_mid");
	if (!m2Attribs.edge_mid)
	{
		std::cout << "Failed to add edge_mid attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.halfedge_volume_connection = add_attribute<Dart, CMap2::HalfEdge>(m2, "volume_connection");
	if (!m2Attribs.halfedge_volume_connection)
	{
		std::cout << "Failed to add halfedge_volume_connection attribute to map2" << std::endl;
		return false;
	}

	return true;
}

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	bool res = true;
	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		value<Dart>(g, gAttribs.vertex_contact_surface, v) = Dart();
		return true;
	});

	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		switch (degree(g, v))
		{
		case 1:
			build_contact_surface_1(g, gAttribs, m2, m2Attribs, v);
			break;
		case 2:
			build_contact_surface_2(g, gAttribs, m2, m2Attribs, v);
			break;
		case 3:
			build_contact_surface_3(g, gAttribs, m2, m2Attribs, v);
			break;
		default:
			build_contact_surface_n(g, gAttribs, m2, m2Attribs, v);
			break;
		}
		return res;
	});
	return res;
}

void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	Dart d = add_face(m2, 4, true).dart;

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d;
	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(d)) = v;

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(v.dart)) = d;
}

void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, d1);
	phi2_sew(m2, phi1(m2, d0), phi_1(m2, d1));
	phi2_sew(m2, phi<11>(m2, d0), phi<11>(m2, d1));
	phi2_sew(m2, phi_1(m2, d0), phi1(m2, d1));

	index_volume_cells(m2, CMap2::Volume(d0));

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d0;
	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(d0)) = v;

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(v.dart)) = d0;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(alpha1(g, v.dart))) = phi<12>(m2, d0);
}

void build_contact_surface_3(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d2 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, phi1(m2, d1));
	phi2_sew(m2, phi1(m2, d0), d2);
	phi2_sew(m2, phi<11>(m2, d0), phi_1(m2, d2));
	phi2_sew(m2, phi_1(m2, d0), phi<11>(m2, d1));
	phi2_sew(m2, d1, phi1(m2, d2));
	phi2_sew(m2, phi_1(m2, d1), phi<11>(m2, d2));

	index_volume_cells(m2, CMap2::Volume(d0));

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d0;

	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);

	std::vector<Vec3> Ppos;
	Ppos.reserve(3);
	std::vector<Dart> Pdart;
	Pdart.reserve(3);

	foreach_dart_of_orbit(g, v, [&](Dart d) {
		Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		project_on_sphere(p, center, radius);
		Ppos.push_back(p);
		Pdart.push_back(d);
		return true;
	});

	Vec3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
	std::vector<Vec3> Q{center + V * radius, center - V * radius};
	std::vector<Vec3> M{center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
						center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
						center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(Pdart[0])) = d0;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(Pdart[1])) = d1;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(Pdart[2])) = d2;

	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d0)) = M[0];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d1)) = M[1];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d2)) = M[2];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, d0))) = Q[0];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, d0))) = Q[1];

	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(d0)) = v;
}

void build_contact_surface_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	// compute the n points on the sphere
	// generate Delaunay mesh from the n points
	// store the graph branch on their respective delaunay vertex (m2Attribs.dual_vertex_graph_branch)
	// modify connectivity until all vertices are valence 4
	// call dualize_volume

		/// Brute force bête et méchant

	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);

	std::vector<Vec3> Ppos;
	std::vector<Dart> Pdart;
	std::vector<uint32> Pid;

	cgogn::io::SurfaceImportData surface_data;

	foreach_dart_of_orbit(g, v, [&](Dart d) -> bool{
		Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		project_on_sphere(p, center, radius);
		uint32 vertex_id = new_index<CMap2::Vertex>(m2);
		Ppos.push_back(p);
		Pid.push_back(vertex_id);
		(*m2Attribs.vertex_position)[vertex_id] = p;
		(*m2Attribs.dual_vertex_graph_branch)[vertex_id] = d;
		surface_data.vertices_id_.push_back(vertex_id);
		return true;
	});

	std::vector<uint32> indices;
	for(uint32 i = 0; i < Ppos.size() - 2; ++i)
	{
		for(uint32 j = i + 1; j < Ppos.size() - 1; ++j)
		{
			for(uint32 k = j + 1; k < Ppos.size(); ++k)
			{
				Vec3 t0 = Ppos[j] - Ppos[i];
				Vec3 t1 = Ppos[k] - Ppos[i];
				Vec3 n = t0.cross(t1);
				int sign = 0;

				for(uint32 m = 0; m < Ppos.size(); ++m)
				{
					if(m == i || m == j || m == k)
						continue;

					Vec3 vec = Ppos[m] - Ppos[i];
					Scalar d = vec.dot(n);

					if(!sign)
						sign = (d < 0? -1: 1);
					else
					{
						if(sign != (d < 0? -1 :1))
						{
							sign = 0;
							break;
						}
					}
				}

				if(sign != 0)
				{
					if(sign == 1)
					{
						indices.push_back(Pid[j]);
						indices.push_back(Pid[i]);
						indices.push_back(Pid[k]);
					}
					else
					{
						indices.push_back(Pid[i]);
						indices.push_back(Pid[j]);
						indices.push_back(Pid[k]);
					}

					surface_data.faces_nb_vertices_.push_back(3);
					surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
												indices.end());
					indices.clear();
				}

			}
		}
	}

	Dart vol_dart = convex_hull(m2, surface_data);

	index_volume_cells(m2, CMap2::Volume(vol_dart));
	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(vol_dart)) = v;

	vol_dart = remesh(m2, CMap2::Volume(vol_dart), m2Attribs);
	dualize_volume(m2, CMap2::Volume(vol_dart), m2Attribs, g, gAttribs);

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = vol_dart;

	return;
}

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	bool res = true;
	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		if (degree(g, v) > 2)
			res = create_intersection_frame_n(g, gAttribs, m2, m2Attribs, v);
		return res;
	});
	return res;
}

bool create_intersection_frame_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 Graph::Vertex v)
{
	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	foreach_dart_of_orbit(g, v, [&](Dart d) -> bool {
		Dart d0 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(d));
		Dart d1 = phi<11>(m2, d0);

		Vec3 R, S, T, diag, temp;
		T = (value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d))) - center).normalized();
		diag = (value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d1)) -
				value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d0)))
				   .normalized();
		R = diag.cross(T).normalized();
		S = T.cross(R).normalized();

		Mat3& f = value<Mat3>(g, gAttribs.halfedge_frame, Graph::HalfEdge(d));
		f.col(0) = R;
		f.col(1) = S;
		f.col(2) = T;

		return true;
	});
	return true;
}

bool propagate_frames(const Graph& g, GAttributes& gAttribs, const GData& gData, CMap2& m2)
{
	for (auto& branch : gData.branches)
	{
		if (degree(g, Graph::Vertex(branch.first.dart)) > 1)
		{
			if (degree(g, Graph::Vertex(branch.second.dart)) == 1)
				propagate_frame_n_1(g, gAttribs, branch.first);
			else
				propagate_frame_n_n(g, gAttribs, m2, branch.first);
		}
		else
			propagate_frame_n_1(g, gAttribs, branch.second);
	}
	return true;
}

void propagate_frame_n_1(const Graph& g, GAttributes& gAttribs, Graph::HalfEdge h_from_start)
{
	Dart d_from = h_from_start.dart;
	Dart d_to = alpha0(g, d_from);
	Graph::HalfEdge h_from(d_from);
	Graph::HalfEdge h_to(d_to);
	Graph::Vertex v_from(d_from);
	Graph::Vertex v_to(d_to);

	uint32 valence = degree(g, v_to);

	Mat3 U = value<Mat3>(g, gAttribs.halfedge_frame, h_from);
	Mat3 U_;

	while (valence == 2)
	{
		Dart next_d_from = alpha1(g, d_to);
		Dart next_d_to = alpha0(g, next_d_from);
		Graph::HalfEdge next_h_from(next_d_from);
		Graph::HalfEdge next_h_to(next_d_to);
		Graph::Vertex next_v_from(next_d_from);
		Graph::Vertex next_v_to(next_d_to);

		Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
		Scalar c1, c2;

		Ri = U.col(0);
		Ti = U.col(2);
		v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
		c1 = v1.dot(v1);
		RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
		TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;

		v3 = value<Vec3>(g, gAttribs.vertex_position, next_v_to) - value<Vec3>(g, gAttribs.vertex_position, v_to);
		Ti1 = (v1.normalized() + v3.normalized()).normalized();
		v2 = Ti1 - TiL;
		c2 = v2.dot(v2);
		Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
		Si1 = Ti1.cross(Ri1);

		U.col(0) = Ri1;
		U.col(1) = Si1;
		U.col(2) = Ti1;

		U_.col(0) = U.col(0);
		U_.col(1) = -U.col(1);
		U_.col(2) = -U.col(2);

		value<Mat3>(g, gAttribs.halfedge_frame, h_to) = U_;
		value<Mat3>(g, gAttribs.halfedge_frame, next_h_from) = U;

		d_from = next_d_from;
		d_to = next_d_to;
		v_from = next_v_from;
		v_to = next_v_to;
		h_from = next_h_from;
		h_to = next_h_to;

		valence = degree(g, v_to);
	}

	if (valence == 1)
	{
		Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
		Scalar c1, c2;

		Ri = U.col(0);
		Ti = U.col(2);
		v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
		c1 = v1.dot(v1);
		RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
		TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;
		Ti1 = v1.normalized();
		v2 = Ti1 - TiL;
		c2 = v2.dot(v2);
		Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
		Si1 = Ti1.cross(Ri1);

		U.col(0) = Ri1;
		U.col(1) = Si1;
		U.col(2) = Ti1;

		U_.col(0) = U.col(0);
		U_.col(1) = -U.col(1);
		U_.col(2) = -U.col(2);

		value<Mat3>(g, gAttribs.halfedge_frame, h_to) = U_;
	}
}

bool propagate_frame_n_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::HalfEdge h_from_start)
{
	Dart d_from = h_from_start.dart;
	Dart d_to = alpha0(g, d_from);
	Graph::HalfEdge h_from(d_from);
	Graph::HalfEdge h_to(d_to);
	Graph::Vertex v_from(d_from);
	Graph::Vertex v_to(d_to);

	uint32 valence = degree(g, v_to);

	Mat3 U = value<Mat3>(g, gAttribs.halfedge_frame, h_from);
	Mat3 U_;

	uint32 nb_e = 0;
	Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
	Scalar c1, c2;

	while (valence == 2)
	{
		Dart next_d_from = alpha1(g, d_to);
		Dart next_d_to = alpha0(g, next_d_from);
		Graph::HalfEdge next_h_from(next_d_from);
		Graph::HalfEdge next_h_to(next_d_to);
		Graph::Vertex next_v_from(next_d_from);
		Graph::Vertex next_v_to(next_d_to);

		nb_e++;

		Ri = U.col(0);
		Ti = U.col(2);
		v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
		c1 = v1.dot(v1);
		RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
		TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;

		v3 = value<Vec3>(g, gAttribs.vertex_position, next_v_to) - value<Vec3>(g, gAttribs.vertex_position, v_to);
		Ti1 = (v1.normalized() + v3.normalized()).normalized();
		v2 = Ti1 - TiL;
		c2 = v2.dot(v2);
		Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
		Si1 = Ti1.cross(Ri1);

		U.col(0) = Ri1;
		U.col(1) = Si1;
		U.col(2) = Ti1;

		value<Mat3>(g, gAttribs.halfedge_frame, next_h_from) = U;

		d_from = next_d_from;
		d_to = next_d_to;
		v_from = next_v_from;
		v_to = next_v_to;
		h_from = next_h_from;
		h_to = next_h_to;

		valence = degree(g, v_to);
	}

	Ri = U.col(0);
	Ti = U.col(2);
	v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
	c1 = v1.dot(v1);
	RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
	TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;
	Ti1 = v1.normalized();
	v2 = Ti1 - TiL;
	c2 = v2.dot(v2);
	Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
	Si1 = Ti1.cross(Ri1);

	U.col(0) = Ri1;
	U.col(1) = Si1;
	U.col(2) = Ti1;

	U_.col(0) = U.col(0);
	U_.col(1) = -U.col(1);
	U_.col(2) = -U.col(2);

	Vec3 X = (U_.col(0) + U_.col(1)).normalized();
	Mat3 UE = value<Mat3>(g, gAttribs.halfedge_frame, h_to);
	Vec3 RE = UE.col(0), SE = UE.col(1);
	bool A = (RE.dot(X) >= 0);
	bool B = (SE.dot(X) >= 0);

	uint32 nb_shifts = 0;
	if (!A && B)
		nb_shifts = 1;
	else if (!A && !B)
		nb_shifts = 2;
	else if (A && !B)
		nb_shifts = 3;
	if (nb_shifts > 0)
	{
		shift_frame(UE, nb_shifts);
		value<Mat3>(g, gAttribs.halfedge_frame, h_to) = UE;
		Dart csface = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h_to);
		for (uint32 i = 0; i < nb_shifts; ++i)
			csface = phi1(m2, csface);
		value<Dart>(g, gAttribs.halfedge_contact_surface_face, h_to) = csface;
	}

	if (nb_e > 0)
	{
		Scalar cos = UE.col(0).dot(U_.col(0));
		Scalar angle = cos > 1 ? std::acos(1) : std::acos(cos);
		Scalar angle_step = -angle / Scalar(nb_e + 1);

		d_from = h_from_start.dart;
		d_to = alpha0(g, d_from);
		h_from.dart = d_from;
		h_to.dart = d_to;
		v_from.dart = d_from;
		v_to.dart = d_to;

		valence = degree(g, v_to);
		uint32 step = 0;
		while (valence == 2)
		{
			step++;

			Dart next_d_from = alpha1(g, d_to);
			Dart next_d_to = alpha0(g, next_d_from);
			Graph::HalfEdge next_h_from(next_d_from);
			Graph::HalfEdge next_h_to(next_d_to);
			Graph::Vertex next_v_from(next_d_from);
			Graph::Vertex next_v_to(next_d_to);

			U = value<Mat3>(g, gAttribs.halfedge_frame, next_h_from);
			Eigen::AngleAxisd rot(-angle_step * step, U.col(2));
			U.col(0) = rot * U.col(0);
			U.col(1) = U.col(2).cross(U.col(0));

			U_.col(0) = U.col(0);
			U_.col(1) = -U.col(1);
			U_.col(2) = -U.col(2);

			value<Mat3>(g, gAttribs.halfedge_frame, h_to) = U_;
			value<Mat3>(g, gAttribs.halfedge_frame, next_h_from) = U;

			d_from = next_d_from;
			d_to = next_d_to;
			v_from = next_v_from;
			v_to = next_v_to;
			h_from = next_h_from;
			h_to = next_h_to;

			valence = degree(g, v_to);
		}
	}

	return true;
}

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const Graph& g, const GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		CMap2::Volume contact_surface(value<Dart>(g, gAttribs.vertex_contact_surface, v));

		const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
		value<Vec3>(m2, m2Attribs.volume_center, contact_surface) = center;

		Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);

		if (degree(g, v) < 3)
		{
			Graph::HalfEdge h(v.dart);
			Dart csf = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h);
			Mat3 frame = value<Mat3>(g, gAttribs.halfedge_frame, h);

			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(csf)) = center - frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, csf))) = center + frame.col(0) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi<11>(m2, csf))) =
				center + frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, csf))) = center - frame.col(0) * radius;
		}

		foreach_incident_edge(m2, contact_surface, [&](CMap2::Edge e) -> bool {
			std::vector<CMap2::Vertex> vertices = incident_vertices(m2, e);
			Vec3 mid = 0.5 * (value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) +
							  value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]));
			project_on_sphere(mid, center, radius);
			//mid = center + ((value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) - center) +
			//				  (value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]) - center));
			value<Vec3>(m2, m2Attribs.edge_mid, e) = mid;
			return true;
		});
		return true;
	});

	return true;
}

//bool set_contact_surfaces_geometry_from_surface(const Graph& g, const GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
//								   const CMap2& surface)
//{
//
//	std::shared_ptr<CMap2::Attribute<Vec3>> surface_vertex_position =
//		get_attribute<Vec3, CMap2::Vertex>(surface, "position");
//
//	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
//		CMap2::Volume contact_surface(value<Dart>(g, gAttribs.vertex_contact_surface, v));
//
//		const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
//		value<Vec3>(m2, m2Attribs.volume_center, contact_surface) = center;
//
//		Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);
//
//		if (degree(g, v) < 3)
//		{
//			Graph::HalfEdge h(v.dart);
//			Dart csf = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h);
//			Mat3 frame = value<Mat3>(g, gAttribs.halfedge_frame, h);
//			Vec3 inter0, inter1, inter2, inter3;
//			intersection_surface(surface, surface_vertex_position.get(), center, -frame.col(1), &inter0);
//			intersection_surface(surface, surface_vertex_position.get(), center, frame.col(0), &inter1);
//			intersection_surface(surface, surface_vertex_position.get(), center, frame.col(1), &inter2);
//			intersection_surface(surface, surface_vertex_position.get(), center, frame.col(1), &inter3);
//			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(csf)) = center - frame.col(1) * radius;
//			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, csf))) = center + frame.col(0) * radius;
//			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi<11>(m2, csf))) =
//				center + frame.col(1) * radius;
//			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, csf))) = center - frame.col(0) * radius;
//		}
//
//		foreach_incident_edge(m2, contact_surface, [&](CMap2::Edge e) -> bool {
//			std::vector<CMap2::Vertex> vertices = incident_vertices(m2, e);
//			Vec3 mid = 0.5 * (value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) +
//							  value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]));
//			project_on_sphere(mid, center, radius);
//			// mid = center + ((value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) - center) +
//			//				  (value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]) - center));
//			value<Vec3>(m2, m2Attribs.edge_mid, e) = mid;
//			return true;
//		});
//		return true;
//	});
//
//	return true;
//}

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/

bool build_branch_sections(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	parallel_foreach_cell(g, [&](Graph::Edge e) -> bool {
		std::vector<Graph::HalfEdge> halfedges = incident_halfedges(g, e);

		Dart m2f0 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, halfedges[0]);
		Dart m2f1 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, halfedges[1]);

		std::vector<Dart> F0 = {m2f0, phi1(m2, m2f0), phi<11>(m2, m2f0), phi_1(m2, m2f0)};
		std::vector<Dart> F1 = {m2f1, phi1(m2, m2f1), phi<11>(m2, m2f1), phi_1(m2, m2f1)};

		Dart m3d = add_branch_section(m3);
		std::vector<Dart> D0 = {m3d, phi<2321>(m3, m3d), phi<23212321>(m3, m3d), phi<111232>(m3, m3d)};
		std::vector<Dart> D1 = {phi<2112>(m3, D0[0]), phi<2112>(m3, D0[1]), phi<2112>(m3, D0[2]), phi<2112>(m3, D0[3])};

		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);

		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[0])) = phi<11>(m3, D1[1]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[1])) = phi<11>(m3, D1[0]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[2])) = phi<11>(m3, D1[3]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[3])) = phi<11>(m3, D1[2]);

		return true;
	});

	return true;
}

bool sew_branch_sections(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		if (is_incident_to_boundary(m2, e))
			return true;

		std::vector<CMap2::HalfEdge> halfedges = incident_halfedges(m2, e);
		sew_volumes(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, halfedges[0]),
					phi1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, halfedges[1])));
		return true;
	});

	close(m3, false);
	return true;
}

bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	auto m3pos = cgogn::add_attribute<Vec3, CMap3::Vertex>(m3, "position");

	parallel_foreach_cell(m2, [&](CMap2::Volume v) -> bool {
		Dart m3d = phi_1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(v.dart)));
		value<Vec3>(m3, m3pos, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.volume_center, v);
		return true;
	});

	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		Dart m3d = phi1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(e.dart)));
		value<Vec3>(m3, m3pos, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.edge_mid, e);
		return true;
	});

	CellMarker<CMap3, CMap3::Vertex> cm(m3);
	for (Dart m2d = m2.begin(), end = m2.end(); m2d != end; m2d = m2.next(m2d))
	{
		if (!is_boundary(m2, m2d))
		{
			Dart m3d = value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(m2d));
			CMap3::Vertex m3v(m3d);
			if (!cm.is_marked(m3v))
			{
				cm.mark(m3v);
				value<Vec3>(m3, m3pos, m3v) = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, m2d)));
			}
		}
	}


	return true;
}


/*****************************************************************************/
/* utils				                                                     */
/*****************************************************************************/

bool dijkstra_topo(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
				   std::shared_ptr<CMap2::Attribute<uint32>> dist)
{
	foreach_incident_vertex(m2, CMap2::Volume(v0.dart), [&](CMap2::Vertex v) -> bool {
				value<Dart>(m2, previous, v) = Dart();
				value<uint32>(m2, dist, v) = UINT_MAX;
				return true;
			});

	DartMarker visited(m2);

	std::vector<CMap2::Vertex> vertices = {v0};
	value<uint32>(m2, dist, v0) = 0;
	CMap2::Vertex v_act; 
	uint32 dist_act;
	while(vertices.size())
	{
		v_act = vertices[vertices.size() - 1];
		vertices.pop_back();
		dist_act = value<uint32>(m2, dist, v_act) + 1;

		std::vector<CMap2::Vertex> neighbors;
		foreach_dart_of_orbit(m2, v_act, [&](Dart d) -> bool {
			if(!visited.is_marked(d))
				{ 
					Dart d2 = phi2(m2, d);
					visited.mark(d);
					visited.mark(d2);

					uint32 dist_2 = value<uint32>(m2, dist, CMap2::Vertex(d2));
					if(dist_act < dist_2)
					{
						value<uint32>(m2, dist, CMap2::Vertex(d2)) = dist_act;
						value<Dart>(m2, previous, CMap2::Vertex(d2)) = d;
						neighbors.push_back(CMap2::Vertex(d2));
					}
				}

			return true;
		});


		vertices.insert(vertices.begin(), neighbors.begin(), neighbors.end());

	}

	return false;
}

Dart convex_hull(CMap2& m2, const cgogn::io::SurfaceImportData& surface_data)
{
	auto darts_per_vertex = add_attribute<std::vector<Dart>, CMap2::Vertex>(m2, "__darts_per_vertex");

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(16u);

	Dart volume_dart;
	std::vector<Dart> all_darts;
	for (std::size_t i = 0u, end = surface_data.faces_nb_vertices_.size(); i < end; ++i)
	{
		uint32 nbv = surface_data.faces_nb_vertices_[i];

		vertices_buffer.clear();

		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = surface_data.faces_vertex_indices_[faces_vertex_index++];
			vertices_buffer.push_back(idx);
		}

		if (nbv > 2u)
		{
			CMap1::Face f = add_face(static_cast<CMap1&>(m2), nbv, false);
			Dart d = f.dart;
			for (uint32 j = 0u; j < nbv; ++j)
			{
				all_darts.push_back(d);
				const uint32 vertex_index = vertices_buffer[j];
				set_index<CMap2::Vertex>(m2, d, vertex_index);
				(*darts_per_vertex)[vertex_index].push_back(d);
				d = phi1(m2, d);
			}
		}
	}

	for (Dart d : all_darts)
	{
		if (phi2(m2, d) == d)
		{
			uint32 vertex_index = index_of(m2, CMap2::Vertex(d));

			const std::vector<Dart>& next_vertex_darts =
				value<std::vector<Dart>>(m2, darts_per_vertex, CMap2::Vertex(phi1(m2, d)));
			bool phi2_found = false;

			for (auto it = next_vertex_darts.begin(); it != next_vertex_darts.end() && !phi2_found; ++it)
			{
				if (index_of(m2, CMap2::Vertex(phi1(m2, *it))) == vertex_index)
				{
					if (phi2(m2, *it) == *it)
					{
						phi2_sew(m2, d, *it);
						phi2_found = true;
					}
				}
			}
		}
	}

	remove_attribute<CMap2::Vertex>(m2, darts_per_vertex);
	return all_darts[0];
}

Vec3 slerp(Vec3 A, Vec3 B, Scalar alpha, bool in)
{
	Scalar phi = geometry::angle(A, B) - (in? 0 : 2*M_PI);
	Scalar s0 = std::sin(phi*(1-alpha));
	Scalar s1 = std::sin(phi*alpha);
	Scalar s2 = std::sin(phi);
	Vec3 sl = A * (s0 / s2);
	sl += B * (s1 / s2);
	return sl;
}


Scalar angle_on_sphere(Vec3 A, Vec3 B, Vec3 C)
{
	Vec3 sB = slerp(A, B, 0.01, true);
	Vec3 sC = slerp(A, C, 0.01, true);
	Vec3 AB = sB - A;
	Vec3 AC = sC - A;
	return geometry::angle(AB, AC);
}

Scalar edge_max_angle(CMap2& m2, CMap2::Edge e, M2Attributes& m2Attribs)
{
	Dart ed0 = e.dart;
	Dart ed1 = phi2(m2, ed0);

	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed0)));

	Scalar a0 = angle_on_sphere(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed1)));
	a0 += angle_on_sphere(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed1)));

	Scalar a1 = angle_on_sphere(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed0)));
	a1 += angle_on_sphere(A, B, C);

	return std::max(a0, a1);
}

Scalar min_cut_angle(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, M2Attributes& m2Attribs)
{
	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v0.dart)));
	Scalar a0 = angle_on_sphere(A, B, C);
	
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v1.dart)));
	Scalar a1 = angle_on_sphere(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	B = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v1.dart)));

	Scalar a2 = angle_on_sphere(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v0.dart)));
	Scalar a3 = angle_on_sphere(A, B, C);

	return std::min({a0, a1, a2, a3});
}

Vec3 spherical_barycenter(std::vector<Vec3>& points, uint32 iterations)
{
	uint32 nb_pts = uint32(points.size());
	std::vector<Vec3> pts0(nb_pts);
	std::vector<Vec3> pts1(nb_pts);

	for (uint32 i = 0; i < nb_pts; ++i)
		pts0[i] = points[i];

	std::vector<Vec3> readvec = pts0;
	std::vector<Vec3> writevec = pts1;
	std::vector<Vec3> tempvec;

	for (uint32 it = 0; it < iterations; ++it)
	{
		for (uint32 i = 0; i < nb_pts; ++i)
		{
			writevec[i] = slerp(readvec[i], readvec[(i + 1) % nb_pts], Scalar(0.5), true);
		}

		tempvec = writevec;
		writevec = readvec;
		readvec = tempvec;
	}
		
	Vec3 bary = {0, 0, 0};
	for (uint32 i = 0; i < nb_pts; ++i)
	{
		bary += readvec[i];
	}
	bary /= nb_pts;

	Vec3 normal = {0, 0, 0};
	for (uint32 i = 0; i < nb_pts; ++i)
	{
		normal += points[i].cross(points[(i + 1) % nb_pts]);
	}
	normal.normalize();

	if (normal.dot(bary) < 0)
		bary *= -1;

	return bary;
}


Dart remesh(CMap2& m2, CMap2::Volume vol, M2Attributes& m2Attribs)
{
	Dart vol_dart = vol.dart;


	std::vector<CMap2::Vertex> valence_sup4;
	std::vector<CMap2::Vertex> valence_3;

	//	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
	//	std::cout << "v3:" << v << " " << value<Vec3>(m2, m2Attribs.vertex_position, v)[0] << " "
	//			  << value<Vec3>(m2, m2Attribs.vertex_position, v)[1] << " "
	//			  << value<Vec3>(m2, m2Attribs.vertex_position, v)[2] << std::endl;
	//	return true;
	//});

	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
				uint32 valence = degree(m2, v);
				if(valence == 3) valence_3.push_back(v);
				if(valence > 4) valence_sup4.push_back(v);

				return true;
			});

	if(valence_sup4.size())
	{
		auto edge_angle_max = add_attribute<Scalar, CMap2::Edge>(m2, "angle_max");

		std::vector<CMap2::Edge> edges_n_n;
		std::vector<CMap2::Edge> edges_n_4;
		std::vector<CMap2::Edge> candidate_edges;

		do
		{
			edges_n_n.clear();
			edges_n_4.clear();
			foreach_incident_edge(m2, vol, [&](CMap2::Edge e) -> bool {
						auto vertices = incident_vertices(m2, e);

						uint32 deg_0 = degree(m2, vertices[0]); 
						uint32 deg_1 = degree(m2, vertices[1]); 
						uint32 deg_min = deg_0 < deg_1 ? deg_0 : deg_1;

						if(deg_min > 4)
						{
							edges_n_n.push_back(e);
						}
						else if(deg_min == 4 && deg_0 + deg_1 > 8)
						{
							edges_n_4.push_back(e);
						}

						return true;
					});

			candidate_edges = edges_n_n.size()? edges_n_n : edges_n_4;
			if (!candidate_edges.size())
				break;
			for(CMap2::Edge e : candidate_edges)
			{
				value<Scalar>(m2, edge_angle_max, e) = edge_max_angle(m2, e, m2Attribs);
			}

			std::sort(candidate_edges.begin(), candidate_edges.end(), [&](CMap2::Edge e0, CMap2::Edge e1) -> bool
			{
				return value<Scalar>(m2, edge_angle_max, e0) < value<Scalar>(m2, edge_angle_max, e1);
			});


			CMap2::Edge prime_edge = candidate_edges[0];
			auto neigh_vertices = incident_vertices(m2, prime_edge);


			vol_dart = phi_1(m2, prime_edge.dart);
			merge_incident_faces(m2, prime_edge, true);

		} while (candidate_edges.size() > 1);
		remove_attribute<CMap2::Edge>(m2, edge_angle_max);

	}
	
	std::vector<std::pair<std::pair<CMap2::Vertex, CMap2::Vertex>, Scalar>> candidate_vertices_pairs;

	valence_3.clear();
	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
		if (degree(m2, v) == 3)
			valence_3.push_back(v);
		return true;
	});
	while (valence_3.size())
	{
		candidate_vertices_pairs.clear();
		std::vector<CMap2::Vertex> verts_3;
		verts_3.reserve(10);
		foreach_incident_face(m2, vol, [&](CMap2::Face f) -> bool {
			verts_3.clear();
			foreach_incident_vertex(m2, f, [&](CMap2::Vertex v) -> bool {
				if (degree(m2, v) == 3) 
					verts_3.push_back(v);

				return true;
			});

			if (verts_3.size() >= 2)
				for(uint32 i = 0; i < verts_3.size() - 1; ++i)
				{
					for(uint32 j = i + 1; j < verts_3.size(); ++j)
					{
						candidate_vertices_pairs.push_back({{verts_3[i], verts_3[j]}, min_cut_angle(m2, verts_3[i], verts_3[j], m2Attribs)});
					}
				}
			return true;
		});

		if(candidate_vertices_pairs.size())
		{
			std::sort(candidate_vertices_pairs.begin(), candidate_vertices_pairs.end(), [&](auto pair0, auto pair1) -> bool
			{
				return pair0.second < pair1.second;
			});

			auto pair = candidate_vertices_pairs[0].first;

			cut_face(m2, pair.first, pair.second);


		}
		else
		{
			uint32 max_path_length = 0;
			std::shared_ptr<CMap2::Attribute<Dart>> max_previous;
			std::shared_ptr<CMap2::Attribute<uint32>> max_dist;
			std::pair<CMap2::Vertex, CMap2::Vertex> max_shortest_path;
			for (uint32 i = 0; i < valence_3.size(); ++i)
			{
				CMap2::Vertex v = valence_3[i];
				std::shared_ptr<CMap2::Attribute<Dart>> previous =
					add_attribute<Dart, CMap2::Vertex>(m2, "previous" + std::to_string(v.dart.index));
				std::shared_ptr<CMap2::Attribute<uint32>> dist =
					add_attribute<uint32, CMap2::Vertex>(m2, "dist" + std::to_string(v.dart.index));

				dijkstra_topo(m2, v, previous, dist);

				uint32 curr_min = UINT32_MAX;
				CMap2::Vertex curr_min_vert;
				for(CMap2::Vertex v2 : valence_3)
				{
					if(index_of<CMap2::Vertex>(m2, v) == index_of<CMap2::Vertex>(m2, v2))
						continue;

					uint32 new_min = value<uint32>(m2, dist, CMap2::Vertex(v2));
					if(new_min < curr_min)
					{
						
						curr_min = new_min; // value<uint32>(m2, dist, CMap2::Vertex(v2)); ????
						curr_min_vert = v2;
					}
				}	
				
				if(curr_min > max_path_length)
				{
					if(max_previous) 
						remove_attribute<CMap2::Vertex>(m2, max_previous);
					if(max_dist) 
						remove_attribute<CMap2::Vertex>(m2, max_dist);

					max_previous = previous;
					max_dist = dist;
					max_shortest_path = {v, curr_min_vert};
					max_path_length = curr_min;
				}

			}

			// building path
			std::vector<Dart> path;
			path.reserve(max_path_length);
			Dart dp = value<Dart>(m2, max_previous, max_shortest_path.second);
			do
			{
				path.push_back(dp);
				dp = value<Dart>(m2, max_previous, CMap2::Vertex(dp));
			} while (dp.index != INVALID_INDEX);

			if (!(path.size() % 2))
			{
				std::vector<CMap2::Edge> first_edges;
				first_edges.push_back(CMap2::Edge(phi_1(m2, path.front())));
				first_edges.push_back(CMap2::Edge(phi2(m2, phi1(m2, path.front()))));
				first_edges.push_back(CMap2::Edge(phi1(m2, path.back())));
				first_edges.push_back(CMap2::Edge(phi2(m2, phi_1(m2, path.back()))));

				std::vector<std::pair<Scalar, uint32>> angles;
				angles.push_back({edge_max_angle(m2, first_edges[0], m2Attribs), 0});
				angles.push_back({edge_max_angle(m2, first_edges[1], m2Attribs), 1});
				angles.push_back({edge_max_angle(m2, first_edges[2], m2Attribs), 2});
				angles.push_back({edge_max_angle(m2, first_edges[3], m2Attribs), 3});

				std::sort(angles.begin(), angles.end(),
						  [&](std::pair<Scalar, uint32> p0, std::pair<Scalar, uint32> p1) -> bool {
							  return p0.first < p1.first;
						  });

				CMap2::Vertex v0, v1;
				CMap2::Edge e0;
				switch (angles[0].second)
				{
				case 0:
					v0 = CMap2::Vertex(phi1(m2, path.front()));
					v1 = CMap2::Vertex(phi_1(m2, path.front()));
					e0 = CMap2::Edge(v1.dart);
					path.erase(path.begin());
					break;
				case 1:
					v0 = CMap2::Vertex(phi2(m2, path.front()));
					v1 = CMap2::Vertex(phi<211>(m2, path.front()));
					e0 = CMap2::Edge(phi<21>(m2, path.front()));
					path.erase(path.begin());
					break;

				case 2:
					v0 = CMap2::Vertex(path.back());
					v1 = CMap2::Vertex(phi<11>(m2, path.back()));
					e0 = CMap2::Edge(phi1(m2, path.back()));
					path.pop_back();
					break;

				case 3:
					v0 = CMap2::Vertex(phi_1(m2, phi2(m2, path.back())));
					v1 = CMap2::Vertex(phi<21>(m2, path.back()));
					e0 = CMap2::Edge(v0.dart);
					path.pop_back();
					break;
				default:
					break;
				}

				 vol_dart = cut_face(m2, v0, v1, true).dart;
				 merge_incident_faces(m2, e0, true);		
			}

			for (uint32 i = 0; i < path.size(); ++i)
			{
				if (!(i % 2))
				{

					cut_face(m2, CMap2::Vertex(path[i]), CMap2::Vertex(phi1(m2, path[i])), true).dart;
				}
				else
				{
					vol_dart = phi1(m2, path[i]);

					merge_incident_faces(m2, CMap2::Edge(path[i]), true);
				}
			}
		}

		valence_3.clear();
		foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
			if (degree(m2, v) == 3) 
				valence_3.push_back(v);
			return true;
		});
	}


	return vol_dart;
}

/*****************************************************************************/
/* mesh volume quality                                                       */
/*****************************************************************************/

bool add_quality_attributes(CMap3& m3, M3Attributes& m3Attribs)
{
	m3Attribs.vertex_position = get_attribute<Vec3, CMap3::Vertex>(m3, "position");
	if (!m3Attribs.vertex_position)
	{
		std::cout << "m3 has no vertex position attribute" << std::endl;
		return false;
	}

	m3Attribs.corner_frame = add_attribute<Mat3, CMap3::Vertex2>(m3, "corner_frame");
	if (!m3Attribs.corner_frame)
	{
		std::cout << "Failed to add corner_frame attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.hex_frame = add_attribute<Mat3, CMap3::Volume>(m3, "hex_frame");
	if (!m3Attribs.hex_frame)
	{
		std::cout << "Failed to add hex_frame attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.scaled_jacobian = add_attribute<Scalar, CMap3::Volume>(m3, "scaled_jacobian");
	if (!m3Attribs.scaled_jacobian)
	{
		std::cout << "Failed to add scaled_jacobian attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.jacobian = add_attribute<Scalar, CMap3::Volume>(m3, "jacobian");
	if (!m3Attribs.jacobian)
	{
		std::cout << "Failed to add jacobian attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.max_frobenius = add_attribute<Scalar, CMap3::Volume>(m3, "max_frobenius");
	if (!m3Attribs.jacobian)
	{
		std::cout << "Failed to add max_frobenius attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.mean_frobenius = add_attribute<Scalar, CMap3::Volume>(m3, "mean_frobenius");
	if (!m3Attribs.jacobian)
	{
		std::cout << "Failed to add mean_frobenius attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.color_scaled_jacobian = add_attribute<Vec3, CMap3::Volume>(m3, "color_scaled_jacobian");
	if (!m3Attribs.color_scaled_jacobian)
	{
		std::cout << "Failed to add color_scaled_jacobian attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.color_jacobian = add_attribute<Vec3, CMap3::Volume>(m3, "color_jacobian");
	if (!m3Attribs.color_jacobian)
	{
		std::cout << "Failed to add color_jacobian attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.color_max_frobenius = add_attribute<Vec3, CMap3::Volume>(m3, "color_max_frobenius");
	if (!m3Attribs.color_max_frobenius)
	{
		std::cout << "Failed to add color_max_frobenius attribute to cmap3" << std::endl;
		return false;
	}
	m3Attribs.color_mean_frobenius = add_attribute<Vec3, CMap3::Volume>(m3, "color_mean_frobenius");
	if (!m3Attribs.color_mean_frobenius)
	{
		std::cout << "Failed to add color_maen_frobenius attribute to cmap3" << std::endl;
		return false;
	}
	return true;
}

bool set_hex_frames(CMap3& m3, M3Attributes& m3Attribs)
{
	foreach_cell(m3, [&](CMap3::Volume vol3) -> bool {
		Dart d0 = vol3.dart;

		Dart D[8];
		D[0] = d0;
		D[1] = phi1(m3, d0);
		D[2] = phi1(m3, D[1]);
		D[3] = phi1(m3, D[2]);
		D[4] = phi<211>(m3, d0);
		D[5] = phi<211>(m3, D[1]);
		D[6] = phi<211>(m3, D[2]);
		D[7] = phi<211>(m3, D[3]);

		Vec3 P[8];
		P[0] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[0]));
		P[1] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[1]));
		P[2] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[2]));
		P[3] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[3]));
		P[4] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[4]));
		P[5] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[5]));
		P[6] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[6]));
		P[7] = value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(D[7]));

		value<Mat3>(m3, m3Attribs.hex_frame, vol3)
			<< ((P[0] + P[1] + P[2] + P[3]) / 4 - (P[4] + P[5] + P[6] + P[7]) / 4),
			((P[0] + P[3] + P[4] + P[7]) / 4 - (P[1] + P[2] + P[5] + P[6]) / 4),
			((P[0] + P[1] + P[4] + P[5]) / 4 - (P[2] + P[3] + P[6] + P[7]) / 4);

		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[0])) << (P[1] - P[0]), (P[4] - P[0]), (P[3] - P[0]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[1])) << (P[0] - P[1]), (P[2] - P[1]), (P[5] - P[1]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[2])) << (P[1] - P[2]), (P[3] - P[2]), (P[6] - P[2]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[3])) << (P[0] - P[3]), (P[7] - P[3]), (P[2] - P[3]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[4])) << (P[0] - P[4]), (P[5] - P[4]), (P[7] - P[4]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[5])) << (P[1] - P[5]), (P[6] - P[5]), (P[4] - P[5]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[6])) << (P[2] - P[6]), (P[7] - P[6]), (P[5] - P[6]);
		value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(D[7])) << (P[3] - P[7]), (P[4] - P[7]), (P[6] - P[7]);
 		return true;
	});

	return true;
}

bool compute_scaled_jacobians(CMap3& m3, M3Attributes& m3Attribs, bool add_color)
{
	foreach_cell(m3, [&](CMap3::Volume vol3) -> bool {
		Mat3 frame_h = value<Mat3>(m3, m3Attribs.hex_frame, vol3);
		frame_h.col(0).normalize();
		frame_h.col(1).normalize();
		frame_h.col(2).normalize();

		Scalar jacobian = frame_h.determinant();

		foreach_incident_vertex(m3, vol3, [&](CMap3::Vertex v3) -> bool {
			Mat3 frame = value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(v3.dart));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar temp = frame.determinant();
			jacobian = temp < jacobian ? temp : jacobian;
			
			return true;
		});
		

		value<Scalar>(m3, m3Attribs.scaled_jacobian, vol3) = jacobian;
		if (add_color) 
			value<Vec3>(m3, m3Attribs.color_scaled_jacobian, vol3) = get_quality_color(jacobian);
		return true;
	});

	return true;
}

bool compute_jacobians(CMap3& m3, M3Attributes& m3Attribs, bool add_color)
{
	foreach_cell(m3, [&](CMap3::Volume vol3) -> bool {
		Mat3 frame_h = value<Mat3>(m3, m3Attribs.hex_frame, vol3);


		Scalar jacobian = frame_h.determinant();

		foreach_incident_vertex(m3, vol3, [&](CMap3::Vertex v3) -> bool {
			Mat3 frame = value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(v3.dart));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar temp = frame.determinant();
			jacobian = temp < jacobian ? temp : jacobian;

			return true;
		});

		value<Scalar>(m3, m3Attribs.jacobian, vol3) = jacobian;
		if (add_color)
		{
			value<Vec3>(m3, m3Attribs.color_jacobian, vol3) = get_quality_color(jacobian);
		}
		return true;
	});

	return true;
}

Scalar frame_frobenius(Mat3 frame)
{
	Scalar det = frame.determinant();
	if (det <= std::numeric_limits<Scalar>::min())
		return std::numeric_limits<Scalar>::max();

	Vec3 c0 = frame.col(0);
	Vec3 c1 = frame.col(1);
	Vec3 c2 = frame.col(2);

	Scalar t1 = c0.dot(c0) + c1.dot(c1) + c2.dot(c2);
	Scalar t2 = (c0.cross(c1)).dot(c0.cross(c1)) + (c1.cross(c2)).dot(c1.cross(c2))
		+ (c2.cross(c0)).dot(c2.cross(c0));

	return sqrt(t1 * t2) / (3 * det);
}

bool compute_maximum_aspect_frobenius(CMap3& m3, M3Attributes& m3Attribs, bool add_color)
{
	foreach_cell(m3, [&](CMap3::Volume vol3) -> bool {
		Scalar frobenius;

		foreach_incident_vertex(m3, vol3, [&](CMap3::Vertex v3) -> bool {
			Mat3 frame = value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(v3.dart));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar temp = frame_frobenius(frame);
			frobenius = temp > frobenius ? temp : frobenius;

			return true;
		});

		value<Scalar>(m3, m3Attribs.max_frobenius, vol3) = frobenius;
		if (add_color)
		{
			value<Vec3>(m3, m3Attribs.color_max_frobenius, vol3) = get_quality_color(1 -(frobenius - 1)/2);
		}
		return true;
	});

	return true;
}

bool compute_mean_aspect_frobenius(CMap3& m3, M3Attributes& m3Attribs, bool add_color)
{
	foreach_cell(m3, [&](CMap3::Volume vol3) -> bool {
		Scalar frobenius;

		foreach_incident_vertex(m3, vol3, [&](CMap3::Vertex v3) -> bool {
			Mat3 frame = value<Mat3>(m3, m3Attribs.corner_frame, CMap3::Vertex2(v3.dart));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar temp = frame_frobenius(frame);
			frobenius += temp;

			return true;
		});

		value<Scalar>(m3, m3Attribs.mean_frobenius, vol3) = frobenius / 8.0;
		if (add_color)
		{
			value<Vec3>(m3, m3Attribs.color_mean_frobenius, vol3) = get_quality_color(1 - (frobenius - 1) / 2);
		}
		return true;
	});

	return true;
}


Vec3 get_quality_color(Scalar quality)
{
	return Vec3(1 - 2 * quality, quality, 1 - 2 * std::abs(0.5 - quality));
}


void export_surface_off(CMap2& m2, std::string name)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position = get_attribute<Vec3, CMap2::Vertex>(m2, "position");
	std::shared_ptr<CMap2::Attribute<uint32>> vertex_id = add_attribute<uint32 , CMap2::Vertex>(m2, "vertex_id");
	std::cout << name << std::endl;
	uint32 nb_vertices = nb_cells<CMap2::Vertex>(m2);
	uint32 nb_faces = nb_cells<CMap2::Face>(m2);

	std::ofstream out_file;
	out_file.open(name);
	out_file << "OFF\n";
	out_file << nb_vertices << " " << nb_faces << " " << 0 << "\n";

	uint32 id = 0;
	foreach_cell(m2, [&](CMap2::Vertex v) -> bool {
		Vec3 pos = value<Vec3>(m2, vertex_position, v);
		value<uint32>(m2, vertex_id, v) = id++;
		out_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
		return true;
	});
	//for (uint32 i = 0; i < nb_vertices; ++i)
	//{
	//	Vec3 pos = (*vertex_position)[i];
	//	out_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
	//}

	foreach_cell(m2, [&](CMap2::Face f) -> bool {
		out_file << codegree(m2, f);
		foreach_incident_vertex(m2, f, [&](CMap2::Vertex v) -> bool {
			out_file << " " << value<uint32>(m2, vertex_id, v);
			return true;
		});
		out_file << "\n";
		return true;
	});

	remove_attribute<CMap2::Vertex>(m2, vertex_id);

	out_file.close();
}

} // namespace modeling

} // namespace cgogn

