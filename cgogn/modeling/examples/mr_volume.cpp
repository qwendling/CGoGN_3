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

using MRMesh = cgogn::CPH3;
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

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	vmrm.selected_vertex_parents_ = cgogn::add_attribute<std::array<Vertex, 4>, Vertex>(*m, "parents");
	vmrm.selected_vertex_relative_position_ = cgogn::add_attribute<Vec3, Vertex>(*m, "relative_position");

	MRMesh* cph2 = vmrm.create_cph3(*m, mp.mesh_name(m));

	vs.selected_mesh_ = cph2;
	auto md = mrmp.mesh_data(cph2);
	md->template add_cells_set<Edge>();

	cgogn::index_cells<Mesh::Volume>(*m);
	cgogn::index_cells<Mesh::Edge>(*m);
	cgogn::index_cells<Mesh::Face>(*m);

	mrsr.set_vertex_position(*v1, *cph2, position);
	// mrsr.set_vertex_position(*v1, *cph2, position);
	// mrsr.set_vertex_position(*v2, *cph1, position);

	std::clock_t start;
	double duration;

	start = std::clock();

	cph2->current_level_ = 1;
	vmrm.subdivide(*cph2, position.get());
	cph2->current_level_ = 2;
	vmrm.subdivide(*cph2, position.get());

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	// std::cout << "temps subdivide  : " << duration << std::endl;

	vs.f_keypress = [&](cgogn::ui::View* view, MRMesh* selected_mesh, std::int32_t k,
						cgogn::ui::CellsSet<MRMesh, Vertex>*, cgogn::ui::CellsSet<MRMesh, Edge>*) {
		switch (k)
		{
		case GLFW_KEY_I: {

			if (selected_mesh->disable_volume_subdivision(Volume(cgogn::Dart(0)), true))
			{
				std::cout << "Yes" << std::endl;
			}
			else
			{
				std::cout << "Nope" << std::endl;
			}
			selected_mesh->check_integrity();
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		case GLFW_KEY_J: {
			auto fn = [](Vertex) {};
			cgogn::modeling::butterflySubdivisionVolume(*selected_mesh, 0.0f, {position.get()}, Volume(cgogn::Dart(0)),
														fn, fn, fn);
			selected_mesh->check_integrity();
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		default:
			break;
		}
	};

	return app.launch();
}
