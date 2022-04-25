/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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
#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>

#include <GLFW/glfw3.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/core/types/cmap/phi.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/animation_multiresolution/animation_multiresolution.h>
#include <cgogn/ui/modules/linked_volumes/linked_volumes.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_emr_modeling/volume_emr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Face2 = typename cgogn::mesh_traits<Mesh>::Face2;
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
	app.set_window_title("Shape Matching");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);
	cgogn::ui::AnimationMultiresolution<MRMesh> am(app);
	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::LinkedVolumes<MRMesh> lv(app);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mrmp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&am);
	v1->link_module(&lv);

	/*cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mrmp);
	v2->link_module(&mrsr);
	v2->link_module(&vs);
	v2->link_module(&am);
	v2->link_module(&lv);

	cgogn::ui::View* v3 = app.add_view();
	v3->link_module(&mrmp);
	v3->link_module(&mrsr);
	v3->link_module(&vs);
	v3->link_module(&am);
	v3->link_module(&lv);*/

	/*cgogn::ui::View* v4 = app.add_view();
	v4->link_module(&mrmp);
	v4->link_module(&mrsr);
	v4->link_module(&vs);
	v4->link_module(&am);*/

	app.init_modules();

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

	MRMesh* meca_mesh = vmrm.create_mrmesh(*m, "mecanic");
	MRMesh* geometry_mesh = vmrm.create_mrmesh(*m, "geometry");
	MRMesh* Visu_mesh = vmrm.create_mrmesh(*m, "Visu");

	vmrm.selected_vertex_parents_ = cgogn::add_attribute<std::array<Vertex, 4>, Vertex>(*m, "parents");
	vmrm.selected_vertex_relative_position_ = cgogn::add_attribute<Vec3, Vertex>(*m, "relative_position");

	cgogn::index_cells<Mesh::Volume>(*m);
	cgogn::index_cells<Mesh::Edge>(*m);
	cgogn::index_cells<Mesh::Face>(*m);

	vmrm.subdivide(*meca_mesh, position.get());
	vmrm.subdivide(*meca_mesh, position.get());
	std::vector<Volume> list_cut_volumes;

	std::srand(164512792);
	while (std::rand() / ((RAND_MAX + 1u) / 5) > 1)
	{
		cgogn::foreach_cell(*meca_mesh, [&list_cut_volumes](Volume v) -> bool {
			list_cut_volumes.push_back(v);
			return true;
		});
		for (Volume v : list_cut_volumes)
		{

			int tmp = std::rand() / ((RAND_MAX + 1u) / 7);
			if (tmp == 1)
			{
				meca_mesh->activate_volume_subdivision(v);
			}
		}
		list_cut_volumes.clear();
	}

	/*cgogn::foreach_cell(*geometry_mesh, [&](Face f) -> bool {
		if (is_incident_to_boundary(*geometry_mesh, f))
		{

			geometry_mesh->activate_face_subdivision(f);
		}
		return true;
	});*/
	for (cgogn::Dart d = geometry_mesh->begin(); d != geometry_mesh->end(); d = geometry_mesh->next(d))
	{
		if (is_boundary(*geometry_mesh, phi3(*geometry_mesh, d)))
		{
			cgogn::Dart d2 = phi2(*geometry_mesh, d);
			while (!is_boundary(*geometry_mesh, phi3(*geometry_mesh, d2)))
			{
				d2 = cgogn::phi<32>(*geometry_mesh, d2);
			}
			// First attribute is the one watch for adaptive subdivision
			auto edge_angle = cgogn::geometry::angle(*geometry_mesh, Face2(d), Face2(d2), position.get());
			if (std::abs(edge_angle) > 0.5f)
			{
				geometry_mesh->activate_face_subdivision(Face(d));
			}
		}
	}

	vmrm.changed_connectivity(*meca_mesh, position.get());
	vmrm.changed_connectivity(*geometry_mesh, position.get());

	mrsr.set_vertex_position(*v1, *meca_mesh, nullptr);
	// mrsr.set_vertex_position(*v2, *topo_mesh, nullptr);
	// mrsr.set_vertex_position(*v2, *meca_mesh, position);
	// mrsr.set_vertex_position(*v3, *topo_mesh, nullptr);
	// mrsr.set_vertex_position(*v3, *meca_mesh, nullptr);
	/*mrsr.set_vertex_position(*v4, *cph1, position);
	mrsr.set_vertex_position(*v4, *cph2, position);*/

	std::srand(std::time(nullptr));

	vs.f_keypress = [&](cgogn::ui::View*, MRMesh* selected_mesh, std::int32_t k,
						cgogn::ui::CellsSet<MRMesh, Vertex>* selected_vertices, cgogn::ui::CellsSet<MRMesh, Edge>*) {
		switch (k)
		{
		case GLFW_KEY_1:
			std::cout << "ok vs" << std::endl;
			mrmp.foreach_mesh([&](MRMesh* m, const std::string&) { mrsr.set_vertex_position(*v1, *m, nullptr); });
			mrsr.set_vertex_position(*v1, *selected_mesh, position);
			break;
		case GLFW_KEY_2:
			mrmp.foreach_mesh([&](MRMesh* m, const std::string&) { vmrm.changed_connectivity(*m, position.get()); });
			break;
		case GLFW_KEY_R: {
			selected_vertices->foreach_cell([&](Vertex v) {
				cgogn::value<Vec3>(*selected_mesh, position.get(), v) =
					cgogn::value<Vec3>(*selected_mesh, position.get(), v) + Vec3((rand() / (double)RAND_MAX) - 0.5,
																				 (rand() / (double)RAND_MAX) - 0.5,
																				 (rand() / (double)RAND_MAX) - 0.5);
			});
			break;
		}
		case GLFW_KEY_U:
			vmrm.changed_connectivity(*selected_mesh, position.get());
			break;
		case GLFW_KEY_K:
			if (selected_vertices != nullptr)
			{
				std::vector<Vertex> new_vertices;
				selected_vertices->foreach_cell([&](Vertex e) {
					/*selected_vertices->unselect(e);
					selected_vertices->select(Vertex(phi1(*selected_mesh, e.dart)));*/

					cgogn::foreach_adjacent_vertex_through_edge(*selected_mesh, e, [&](Vertex w) -> bool {
						new_vertices.push_back(w);
						return true;
					});
				});
				for (auto v : new_vertices)
				{
					selected_vertices->select(v);
				}

				vs.mesh_provider_->emit_cells_set_changed(selected_mesh, selected_vertices);
			}
			break;
		case GLFW_KEY_B:
			if (selected_vertices != nullptr)
			{
				std::vector<Vertex> new_vertices;
				selected_vertices->foreach_cell([&](Vertex e) {
					selected_vertices->unselect(e);
					selected_vertices->select(Vertex(phi2(*selected_mesh, e.dart)));
				});

				vs.mesh_provider_->emit_cells_set_changed(selected_mesh, selected_vertices);
			}
			break;
		case GLFW_KEY_J:
			if (selected_vertices != nullptr)
			{
				std::vector<Vertex> new_vertices;
				selected_vertices->foreach_cell([&](Vertex e) {
					selected_vertices->unselect(e);
					selected_vertices->select(Vertex(phi3(*selected_mesh, e.dart)));
				});

				vs.mesh_provider_->emit_cells_set_changed(selected_mesh, selected_vertices);
			}
			break;
		case GLFW_KEY_L:
			if (selected_vertices != nullptr)
			{
				selected_vertices->foreach_cell([&](Vertex e) {
					cgogn::foreach_incident_volume(*selected_mesh, e, [&](Volume v) -> bool {
						std::cout << selected_mesh->volume_level(v.dart) << std::endl;
						return true;
					});
				});

				vmrm.changed_connectivity(*selected_mesh, position.get());
			}
			break;
		}
	};

	return app.launch();
}
