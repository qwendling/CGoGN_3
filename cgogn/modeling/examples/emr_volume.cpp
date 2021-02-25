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
#include <cgogn/ui/modules/volume_emr_modeling/volume_emr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;

template <typename T>
using Attribute = typename cgogn::mesh_traits<MRMesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<MRMesh>::Vertex;
using Edge = typename cgogn::mesh_traits<MRMesh>::Edge;
using Face = typename cgogn::mesh_traits<MRMesh>::Face;
using Volume = typename cgogn::mesh_traits<MRMesh>::Volume;

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
	cgogn::ui::VolumeRender<MRMesh> vr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);

	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrmp);
	v1->link_module(&vr);
	v1->link_module(&vs);

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	MRMesh* mrm = vmrm.create_mrmesh(*m, mp.mesh_name(m));
	vs.selected_mesh_ = mrm;

	m->add_resolution();
	mrm->change_resolution_level(1);
	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");

	vmrm.subdivide(*mrm, position.get());
	m->add_resolution();
	mrm->change_resolution_level(2);
	// vmrm.subdivide(*mrm, position.get());

	auto md = mrmp.mesh_data(mrm);
	md->template add_cells_set<Edge>();

	mrmp.set_mesh_bb_vertex_position(mrm, position);

	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);
	cgogn::index_cells<Mesh::Face>(*mrm);

	vr.set_vertex_position(*v1, *mrm, position);

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
			std::cout << "hello" << std::endl;

			break;
		case GLFW_KEY_C:
			std::vector<int> bucket;
			for (uint i = 0; i <= mrm->maximum_level_; i++)
			{
				bucket.push_back(0);
			}
			cgogn::foreach_cell(*mrm, [&](Volume f) -> bool {
				bucket[mrm->volume_level(f.dart)]++;
				return true;
			});
			for (auto i : bucket)
			{
				std::cout << i << std::endl;
			}

			break;
		}
	};

	return app.launch();
}
