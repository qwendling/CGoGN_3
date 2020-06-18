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
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>

#include <GLFW/glfw3.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/animation_multiresolution/animation_multiresolution.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_mr_modeling/volume_mr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

using Mesh = cgogn::CMap3;
using MRMesh = cgogn::CPH3_adaptative;

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
	app.set_window_title("Shape Matching");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);
	cgogn::ui::AnimationMultiresolution<MRMesh> am(app);
	cgogn::ui::VolumeMRModeling vmrm(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mrmp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&am);

	cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mrmp);
	v2->link_module(&mrsr);
	v2->link_module(&vs);
	v2->link_module(&am);

	cgogn::ui::View* v3 = app.add_view();
	v3->link_module(&mrmp);
	v3->link_module(&mrsr);
	v3->link_module(&vs);
	v3->link_module(&am);

	cgogn::ui::View* v4 = app.add_view();
	v4->link_module(&mrmp);
	v4->link_module(&mrsr);
	v4->link_module(&vs);
	v4->link_module(&am);

	app.init_modules();

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

	MRMesh* cph1 = vmrm.create_cph3(*m, mp.mesh_name(m));
	MRMesh* cph2 = vmrm.create_cph3(*m, mp.mesh_name(m));

	vmrm.selected_vertex_parents_ = cgogn::add_attribute<std::array<Vertex, 3>, Vertex>(*m, "parents");
	vmrm.selected_vertex_relative_position_ = cgogn::add_attribute<Vec3, Vertex>(*m, "relative_position");

	cgogn::index_cells<Mesh::Volume>(*m);
	cgogn::index_cells<Mesh::Edge>(*m);
	cgogn::index_cells<Mesh::Face>(*m);

	vmrm.subdivide(*cph2, position.get());
	vmrm.subdivide(*cph2, position.get());
	std::vector<Volume> list_cut_volumes;
	cgogn::foreach_cell(*cph2, [&list_cut_volumes](Volume v) -> bool {
		list_cut_volumes.push_back(v);
		return true;
	});
	for (Volume v : list_cut_volumes)
		cph2->activate_volume_subdivision(v);
	vmrm.changed_connectivity(*cph2, position.get());

	mrsr.set_vertex_position(*v1, *cph1, position);
	mrsr.set_vertex_position(*v1, *cph2, nullptr);
	mrsr.set_vertex_position(*v2, *cph1, nullptr);
	mrsr.set_vertex_position(*v2, *cph2, position);
	/*mrsr.set_vertex_position(*v3, *cph1, position);
	mrsr.set_vertex_position(*v3, *cph2, position);
	mrsr.set_vertex_position(*v4, *cph1, position);
	mrsr.set_vertex_position(*v4, *cph2, position);*/

	std::srand(std::time(nullptr));

	vs.f_keypress = [&](cgogn::ui::View*, MRMesh* selected_mesh, std::int32_t k,
						cgogn::ui::CellsSet<MRMesh, Vertex>* selected_vertices, cgogn::ui::CellsSet<MRMesh, Edge>*) {
		switch (k)
		{
		case GLFW_KEY_R: {
			selected_vertices->foreach_cell([&](Vertex v) {
				cgogn::value<Vec3>(*selected_mesh, position.get(), v) =
					cgogn::value<Vec3>(*selected_mesh, position.get(), v) + Vec3((rand() / (double)RAND_MAX) - 0.5,
																				 (rand() / (double)RAND_MAX) - 0.5,
																				 (rand() / (double)RAND_MAX) - 0.5);
			});
			break;
		}
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
		}
	};

	return app.launch();
}
