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
#include <cgogn/core/types/cmap/phi.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_emr_modeling/volume_emr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;
using EMR_Map3 = cgogn::EMR_Map3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<MRMesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<MRMesh>::Vertex;
using Edge = typename cgogn::mesh_traits<MRMesh>::Edge;
using Face = typename cgogn::mesh_traits<MRMesh>::Face;
using Volume = typename cgogn::mesh_traits<MRMesh>::Volume;
using Dart = cgogn::Dart;

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
	app.set_window_title("EMR Volume");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> vr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);

	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrmp);
	v1->link_module(&vr);
	v1->link_module(&vs);

	/*cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mp);
	v2->link_module(&mrmp);
	v2->link_module(&vr);
	v2->link_module(&vs);*/

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	MRMesh* mrm = vmrm.create_mrmesh(*m, mp.mesh_name(m));
	// MRMesh* mrm2 = vmrm.create_mrmesh(*m, mp.mesh_name(m));
	MRMesh* mrm2 = mrm->get_copy();
	mrm2->parent = mrm;
	mrmp.register_mesh(mrm2, "copy");
	vs.selected_mesh_ = mrm;
	cgogn::index_cells<Mesh::Face>(*mrm);
	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);

	m->add_resolution();
	mrm->change_resolution_level(1);
	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");

	vmrm.subdivide(*mrm, position.get());
	m->add_resolution();
	mrm->change_resolution_level(2);
	vmrm.subdivide(*mrm, position.get());

	auto md = mrmp.mesh_data(mrm);
	md->template add_cells_set<Edge>();

	mrmp.set_mesh_bb_vertex_position(mrm, position);

	vr.set_vertex_position(*v1, *mrm, position);
	vr.set_vertex_position(*v1, *mrm2, nullptr);

	// std::srand(std::time(nullptr));
	std::srand(2124512438);

	vs.f_keypress = [&](cgogn::ui::View* view, MRMesh* selected_mesh, std::int32_t k,
						cgogn::ui::CellsSet<MRMesh, Vertex>* selected_vertices,
						cgogn::ui::CellsSet<MRMesh, Edge>* selected_edges) {
		switch (k)
		{
		case GLFW_KEY_I: {
			MRMesh tmp(*mrm);
			tmp.change_resolution_level(0);
			if (tmp.check_integrity())
				std::cout << "ok lvl 0 integrity" << std::endl;
			cgogn_message_assert(tmp.check_integrity(), "check_integrity failed");
			tmp.change_resolution_level(1);
			if (tmp.check_integrity())
				std::cout << "ok lvl 1 integrity" << std::endl;
			cgogn_message_assert(tmp.check_integrity(), "check_integrity failed");
			tmp.change_resolution_level(2);
			if (tmp.check_integrity())
				std::cout << "ok lvl 2 integrity" << std::endl;
			cgogn_message_assert(tmp.check_integrity(), "check_integrity failed");
			std::cout << "ok check " << std::endl;
			break;
		}
		case GLFW_KEY_E:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex v) {
					std::vector<Edge> vec_edge;
					cgogn::foreach_incident_edge(*selected_mesh, v, [&](Edge e) -> bool {
						vec_edge.push_back(e);

						return true;
					});
					for (Edge& e : vec_edge)
					{
						if (view->shift_pressed())
						{
							selected_mesh->disable_edge_subdivision(e);
						}
						else
						{
							selected_mesh->activate_edge_subdivision(e);
						}
					}
				});
				cgogn_message_assert(selected_mesh->check_integrity(), "check_integrity failed");
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}

			break;
		case GLFW_KEY_F:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex v) {
					std::vector<Face> vec_face;
					cgogn::foreach_incident_face(*selected_mesh, v, [&](Face f) -> bool {
						vec_face.push_back(f);
						return true;
					});
					for (auto& f : vec_face)
					{
						if (view->shift_pressed())
						{
							if (selected_mesh->disable_face_subdivision(f, true))
								std::cout << "ok pour la subdiv de face " << std::endl;
						}
						else
						{
							selected_mesh->activate_face_subdivision(f);
						}
					}
				});
				cgogn_message_assert(selected_mesh->check_integrity(), "check_integrity failed");
				vmrm.changed_connectivity(*selected_mesh, position.get());
			}

			std::cout << "hello" << std::endl;

			break;
		case GLFW_KEY_V:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex v) {
					std::vector<Volume> vec_volume;
					cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
						vec_volume.push_back(w);
						return true;
					});
					for (auto& w : vec_volume)
					{
						if (view->shift_pressed())
						{
							if (selected_mesh->disable_volume_subdivision(w, true))
								std::cout << "ok pour la subdiv de face " << std::endl;
						}
						else
						{
							selected_mesh->activate_volume_subdivision(w);
						}
					}
				});
				vmrm.changed_connectivity(*mrm, position.get());
				vmrm.changed_connectivity(*mrm2, position.get());
			}
			cgogn_message_assert(mrm->check_integrity(), "check_integrity failed");
			cgogn_message_assert(mrm2->check_integrity(), "check_integrity failed");
			std::cout << "hello" << std::endl;

			break;
		case GLFW_KEY_Z: {
			std::vector<Edge> vec_edge;
			std::clock_t start;
			double duration;

			start = std::clock();

			cgogn::foreach_cell(*selected_mesh, [&](Vertex) -> bool { return true; });
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps parcours vertex : " << duration << std::endl;
			start = std::clock();
			cgogn::foreach_cell(*selected_mesh, [&](Edge) -> bool { return true; });
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps parcours edge : " << duration << std::endl;
			start = std::clock();
			cgogn::foreach_cell(*selected_mesh, [&](Face) -> bool { return true; });
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps parcours face : " << duration << std::endl;
			start = std::clock();
			cgogn::foreach_cell(*selected_mesh, [&](Volume) -> bool { return true; });
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps parcours volume : " << duration << std::endl;
			start = std::clock();
			duration = 0;
			for (cgogn::Dart d = selected_mesh->begin(), e = selected_mesh->end(); d != e; d = selected_mesh->next(d))
			{
				start = std::clock();

				if (selected_mesh->edge_level(d) != 0)
				{
					foreach_dart_of_orbit(*selected_mesh, Face(d), [&](Dart) -> bool { return true; });
					duration += (std::clock() - start) / (double)CLOCKS_PER_SEC;
				}
			}
			std::cout << "temps face foreach dart : " << duration << std::endl;
			start = std::clock();
			duration = 0;
			int cpt = 0;
			for (cgogn::Dart d = selected_mesh->begin(), e = selected_mesh->end(); d != e; d = selected_mesh->next(d))
			{
				start = std::clock();

				if (selected_mesh->edge_level(d) != 0)
				{
					foreach_dart_of_orbit(*selected_mesh, Volume(d), [&](cgogn::Dart) -> bool {
						cpt++;
						return true;
					});
					duration += (std::clock() - start) / (double)CLOCKS_PER_SEC;
				}
			}
			std::cout << "temps volume foreach dart : " << duration << std::endl;
			std::cout << "nb dart : " << cpt << std::endl;
			start = std::clock();
			for (cgogn::Dart d = selected_mesh->begin(), e = selected_mesh->end(); d != e; d = selected_mesh->next(d))
			{
				cgogn::phi1(*selected_mesh, d);
			}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps phi1 : " << duration << std::endl;
			start = std::clock();
			for (cgogn::Dart d = selected_mesh->begin(), e = selected_mesh->end(); d != e; d = selected_mesh->next(d))
			{
				phi2(*selected_mesh, d);
			}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps phi2 : " << duration << std::endl;
			start = std::clock();
			for (cgogn::Dart d = selected_mesh->begin(), e = selected_mesh->end(); d != e; d = selected_mesh->next(d))
			{
				phi3(*selected_mesh, d);
			}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps phi3 : " << duration << std::endl;
		}

		break;
		case GLFW_KEY_T: {
			std::vector<Edge> vec_edge;

			cgogn::foreach_cell(*mrm, [&](Edge e) -> bool {
				if ((rand() / (double)RAND_MAX) * 100 < 10)
				{
					vec_edge.push_back(e);
				}
				return true;
			});
			for (auto e : vec_edge)
			{
				mrm->activate_edge_subdivision(e);
			}
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}

		break;
		case GLFW_KEY_U: {
			std::vector<Face> vec_face;

			cgogn::foreach_cell(*mrm, [&](Face f) -> bool {
				if ((rand() / (double)RAND_MAX) * 100 < 10)
				{
					vec_face.push_back(f);
				}
				return true;
			});
			for (auto f : vec_face)
			{
				if (view->shift_pressed())
				{
					selected_mesh->disable_face_subdivision(f, true);
				}
				else
				{
					selected_mesh->activate_face_subdivision(f);
				}
			}
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}

		break;
		case GLFW_KEY_R:
#define CELL_RANDOM Volume
			if (!view->shift_pressed())
			{
				std::vector<CELL_RANDOM> vec_volume;

				cgogn::foreach_cell(*selected_mesh, [&](CELL_RANDOM v) -> bool {
					if ((rand() / (double)RAND_MAX) * 100 < 10)
					{
						vec_volume.push_back(v);
					}
					return true;
				});
				std::clock_t start;
				double duration;
				auto md = mrmp.mesh_data(selected_mesh);
				double diff_volume = md->nb_cells<CELL_RANDOM>();

				start = std::clock();
				for (auto v : vec_volume)
				{
					selected_mesh->activate_volume_subdivision(v);
				}

				duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				md->update_nb_cells();
				std::cout << "temps activate " << md->nb_cells<CELL_RANDOM>() - diff_volume << " volume : " << duration
						  << std::endl;
				start = std::clock();
				vmrm.changed_connectivity(*selected_mesh, position.get());
				duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				std::cout << "temps update topo : " << duration << std::endl;
			}
			else
			{
				std::vector<CELL_RANDOM> vec_volume;

				cgogn::foreach_cell(*selected_mesh, [&](CELL_RANDOM v) -> bool {
					if ((rand() / (double)RAND_MAX) * 100 < 10)
					{
						vec_volume.push_back(v);
					}
					return true;
				});
				std::clock_t start;
				double duration;
				auto md = mrmp.mesh_data(selected_mesh);
				double diff_volume = md->nb_cells<CELL_RANDOM>();

				start = std::clock();
				for (auto v : vec_volume)
				{
					if (selected_mesh->get_dart_visibility(v.dart) > selected_mesh->current_level_)
						continue;
					selected_mesh->disable_volume_subdivision(v, true);
				}

				duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				md->update_nb_cells();
				std::cout << "temps disable " << diff_volume - md->nb_cells<CELL_RANDOM>() << " volume : " << duration
						  << std::endl;
				start = std::clock();
				vmrm.changed_connectivity(*selected_mesh, position.get());
				duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				std::cout << "temps update topo : " << duration << std::endl;
			}

			break;
		case GLFW_KEY_C:
			std::vector<int> bucket;
			for (uint i = 0; i <= mrm2->maximum_level_; i++)
			{
				bucket.push_back(0);
			}
			cgogn::foreach_cell(*mrm2, [&](Edge f) -> bool {
				bucket[mrm2->edge_level(f.dart)]++;
				return true;
			});
			std::cout << "Edge level : " << std::endl;
			for (auto i : bucket)
			{
				std::cout << i << std::endl;
			}

			for (uint i = 0; i <= mrm2->maximum_level_; i++)
			{
				bucket[i] = 0;
			}
			cgogn::foreach_cell(*mrm2, [&](Face f) -> bool {
				bucket[mrm2->face_level(f.dart)]++;
				return true;
			});
			std::cout << "Face level : " << std::endl;
			for (auto i : bucket)
			{
				std::cout << i << std::endl;
			}

			for (uint i = 0; i <= mrm2->maximum_level_; i++)
			{
				bucket[i] = 0;
			}
			cgogn::foreach_cell(*mrm2, [&](Volume f) -> bool {
				bucket[mrm2->volume_level(f.dart)]++;
				return true;
			});
			std::cout << "Volume level : " << std::endl;
			for (auto i : bucket)
			{
				std::cout << i << std::endl;
			}

			break;
		}
	};

	return app.launch();
}
