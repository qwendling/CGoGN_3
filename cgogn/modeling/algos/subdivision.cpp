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

#include "cgogn/modeling/algos/subdivision.h"

namespace cgogn
{

namespace modeling
{

void butterflySubdivisionVolume(CMap3& m, double,
								std::vector<typename mesh_traits<CMap3>::template Attribute<Vec3>*> attributes,
								std::vector<typename CMap3::Volume> vec_v)
{
	using Volume = typename CMap3::Volume;
	using Face = typename CMap3::Face;
	using Edge = typename CMap3::Edge;
	using Vertex = typename CMap3::Vertex;

	std::vector<typename mesh_traits<CMap3>::template Attribute<Vec3>*> attrs;
	for (auto a : attributes)
		if (a)
			attrs.push_back(a);

	CellMarkerStore<CMap3, Edge> cm_edge(m);
	CellMarkerStore<CMap3, Face> cm_face(m);
	CellMarkerStore<CMap3, Volume> cm_volume(m);
	std::queue<std::queue<Vec3>> volume_points, face_points, edge_points;
	std::vector<Dart> edges, faces, volumes;
	std::vector<Dart> p_point, q_point, r_point, s_point, t_point;

	// computing the new vertices's embedding
	for (Volume v : vec_v)
	{
		foreach_dart_of_orbit(m, v, [&](Dart t) -> bool {
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
	}

	subdivideListEdges<CMap3>(m, edges, [&](Vertex v) {
		for (auto a : attrs)
		{
			value<Vec3>(m, a, v) = edge_points.front().front();
			edge_points.front().pop();
		}
		edge_points.pop();
	});

	subdivideListFaces(m, faces, [&](Vertex v) {
		for (auto a : attrs)
		{
			value<Vec3>(m, a, v) = face_points.front().front();
			face_points.front().pop();
		}
		face_points.pop();
	});

	subdivideListVolumes(m, volumes, [&](Vertex v) {
		for (auto a : attrs)
		{
			value<Vec3>(m, a, v) = volume_points.front().front();
			volume_points.front().pop();
		}
		volume_points.pop();
	});
}

} // namespace modeling

} // namespace cgogn
