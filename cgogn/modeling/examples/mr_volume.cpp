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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>

#include <GLFW/glfw3.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_mr_modeling/volume_mr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

using MRMesh = cgogn::CPH3_adaptative;
using Mesh = MRMesh::CMAP;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

using Vec3 = cgogn::geometry::Vec3;

int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " volume_mesh_file" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("MR Volume");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);

	cgogn::ui::VolumeMRModeling vmrm(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mrmp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);

	cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mrmp);
	v2->link_module(&mrsr);
	v2->link_module(&vs);

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	vmrm.selected_vertex_parents_ = cgogn::add_attribute<std::array<Vertex, 3>, Vertex>(*m, "parents");
	vmrm.selected_vertex_relative_position_ = cgogn::add_attribute<Vec3, Vertex>(*m, "relative_position");

	MRMesh* cph1 = vmrm.create_cph3(*m, mp.mesh_name(m));
	MRMesh* cph2 = vmrm.create_cph3(*m, mp.mesh_name(m));

	vs.selected_mesh_ = cph2;
	auto md = mrmp.mesh_data(cph2);
	md->template add_cells_set<Edge>();

	cgogn::index_cells<Mesh::Volume>(*m);
	cgogn::index_cells<Mesh::Edge>(*m);
	cgogn::index_cells<Mesh::Face>(*m);

	mrsr.set_vertex_position(*v1, *cph1, position);
	mrsr.set_vertex_position(*v1, *cph2, nullptr);
	// mrsr.set_vertex_position(*v1, *cph2, position);
	// mrsr.set_vertex_position(*v2, *cph1, position);
	mrsr.set_vertex_position(*v2, *cph2, position);

	std::clock_t start;
	double duration;

	start = std::clock();

	vmrm.subdivide(*cph2, position.get());
	vmrm.subdivide(*cph2, position.get());

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps subdivide  : " << duration << std::endl;

	std::srand(std::time(nullptr));

	vs.f_keypress = [&](cgogn::ui::View* view, MRMesh* selected_mesh, std::int32_t k,
						cgogn::ui::CellsSet<MRMesh, Vertex>* selected_vertices,
						cgogn::ui::CellsSet<MRMesh, Edge>* selected_edges) {
		switch (k)
		{
		case GLFW_KEY_E:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex v) {
					cgogn::foreach_incident_edge(*m, v, [&](Edge e) -> bool {
						if (view->shift_pressed())
						{
							selected_mesh->disable_edge_subdivision(e, true);
						}
						else
						{
							selected_mesh->activate_edge_subdivision(e);
						}
						return true;
					});
				});
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			if (selected_edges != nullptr)
			{
				selected_edges->foreach_cell([&](Edge e) {
					if (view->shift_pressed())
					{
						selected_mesh->disable_edge_subdivision(e, true);
					}
					else
					{
						selected_mesh->activate_edge_subdivision(e);
					}
				});
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			break;
		case GLFW_KEY_F:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex e) {
					std::vector<Face> face_list;
					cgogn::foreach_incident_face(*selected_mesh, e, [&](Face f) -> bool {
						face_list.push_back(f);
						return true;
					});
					for (auto f : face_list)
					{
						if (view->shift_pressed())
						{
							if (selected_mesh->dart_is_visible(f.dart))
								selected_mesh->disable_face_subdivision(f, true, true);
						}
						else
						{
							selected_mesh->activate_face_subdivision(f);
						}
					}
				});

				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			if (selected_edges != nullptr)
			{
				selected_edges->foreach_cell([&](Edge e) {
					std::vector<Face> face_list;
					cgogn::foreach_incident_face(*selected_mesh, e, [&](Face f) -> bool {
						face_list.push_back(f);
						return true;
					});
					for (auto f : face_list)
					{
						if (view->shift_pressed())
						{
							if (selected_mesh->dart_is_visible(f.dart))
								selected_mesh->disable_face_subdivision(f, true, true);
						}
						else
						{
							selected_mesh->activate_face_subdivision(f);
						}
					}
				});
				int nb_rep = 0, nb_visible = 0;
				for (cgogn::Dart d = m->begin(), end = m->end(); d != end; d = m->next(d))
				{
					if (selected_mesh->dart_level(selected_mesh->get_representative(d)) > 0 &&
						selected_mesh->representative_is_visible(d))
						nb_rep++;
					if (selected_mesh->dart_level(d) > 0 &&
						selected_mesh->get_dart_visibility_level(d) < selected_mesh->dart_level(d))
						nb_visible++;
				}
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			break;
		case GLFW_KEY_S:
			cgogn::foreach_cell(*selected_mesh, [&](Face f) -> bool {
				if (is_incident_to_boundary(*selected_mesh, f))
				{

					selected_mesh->activate_face_subdivision(f);
				}
				return true;
			});
			vmrm.changed_connectivity(*selected_mesh, position.get());
			break;
		case GLFW_KEY_V:
			if (selected_vertices != nullptr)
			{

				std::vector<Volume> volume_list;
				selected_vertices->foreach_cell([&](Vertex e) {
					cgogn::foreach_incident_volume(*m, e, [&](Volume v) -> bool {
						volume_list.push_back(v);
						return true;
					});
				});
				for (Volume v : volume_list)
				{
					if (view->shift_pressed())
					{
						if (selected_mesh->dart_is_visible(v.dart))
							selected_mesh->disable_volume_subdivision(v, true);
						std::cout << "hello" << std::endl;
					}
					else
					{
						selected_mesh->activate_volume_subdivision(v);
					}
				}
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			if (selected_edges != nullptr)
			{
				std::vector<Volume> volume_list;
				selected_edges->foreach_cell([&](Edge e) {
					cgogn::foreach_incident_volume(*m, e, [&](Volume v) -> bool {
						volume_list.push_back(v);
						return true;
					});
				});
				for (Volume v : volume_list)
				{
					if (view->shift_pressed())
					{
						if (selected_mesh->dart_is_visible(v.dart))
							selected_mesh->disable_volume_subdivision(v, true);
					}
					else
					{
						selected_mesh->activate_volume_subdivision(v);
					}
				}
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			break;
		case GLFW_KEY_A: {
			std::vector<Edge> list_cut_edges;
			cgogn::foreach_cell(*selected_mesh, [&](Edge e) -> bool {
				list_cut_edges.push_back(e);
				return true;
			});
			for (Edge e : list_cut_edges)
				if (view->shift_pressed())
				{
					// selected_mesh->disable_edge_subdivision(e);
				}
				else
				{
					selected_mesh->activate_edge_subdivision(e);
				}
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		case GLFW_KEY_Q: {
			std::vector<Face> list_cut_faces;
			cgogn::foreach_cell(*selected_mesh, [&list_cut_faces](Face f) -> bool {
				list_cut_faces.push_back(f);
				return true;
			});
			for (Face f : list_cut_faces)
				selected_mesh->activate_face_subdivision(f);
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		case GLFW_KEY_Z: {
			std::vector<Volume> list_cut_volumes;
			cgogn::foreach_cell(*selected_mesh, [&list_cut_volumes](Volume v) -> bool {
				list_cut_volumes.push_back(v);
				return true;
			});
			clock_t start = clock();
			for (Volume v : list_cut_volumes)
				selected_mesh->activate_volume_subdivision(v);
			double time = (clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "time for " << list_cut_volumes.size() << " volume activation : " << time << "s" << std::endl;
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		case GLFW_KEY_L:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex v) {
					cgogn::foreach_incident_face(*selected_mesh, v, [&](Face f) -> bool {
						std::cout << "face level : " << selected_mesh->face_level(f.dart) << std::endl;
						return true;
					});
				});
			}
			if (selected_edges != nullptr)
			{
				selected_edges->foreach_cell([&](Edge e) {
					cgogn::foreach_incident_face(*selected_mesh, e, [&](Face f) -> bool {
						std::cout << "face level : " << selected_mesh->face_level(f.dart) << std::endl;
						std::cout << "face oldest : " << selected_mesh->face_youngest_dart(f.dart).index << std::endl;
						return true;
					});
					std::cout << "edge level : " << selected_mesh->edge_level(e.dart) << std::endl;
					std::cout << "_____________________________________" << std::endl;
				});
				std::cout << "######################################" << std::endl;
			}
			break;
		case GLFW_KEY_R: {
			std::vector<Volume> vec_volume;

			cgogn::foreach_cell(*selected_mesh, [&](Volume v) -> bool {
				if ((rand() / (double)RAND_MAX) * 100 < 10)
				{
					vec_volume.push_back(v);
				}
				return true;
			});
			std::clock_t start;
			double duration;

			start = std::clock();
			for (auto v : vec_volume)
			{
				selected_mesh->activate_volume_subdivision(v);
			}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps activate " << vec_volume.size() << " volume : " << duration << std::endl;
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}

		break;
		case GLFW_KEY_P: {
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex v) { cph2->raise_volume_level(Volume(v.dart)); });
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
		}
		break;
		}
	};

	return app.launch();
}
