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
#include <random>

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

	/*cgogn::modeling::butterflySubdivisionVolumeRegular(*m, 0.0f, {position.get()});
	cgogn::modeling::butterflySubdivisionVolumeRegular(*m, 0.0f, {position.get()});*/

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

	// cph2->current_level_ = 1;
	//  vmrm.subdivide(*cph2, position.get());
	/*cph2->current_level_ = 2;
	vmrm.subdivide(*cph2, position.get());*/
	// cph2->current_level_ = 2;

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	// std::cout << "temps subdivide  : " << duration << std::endl;

	std::vector<Volume> volume_to_subdivided;
	std::vector<Volume> volume_to_simplified;
	std::random_device rd;
	std::mt19937 g(rd());

	cgogn::foreach_cell(*cph2, [&](Volume v) -> bool {
		volume_to_subdivided.push_back(Volume(cph2->volume_youngest_dart(v.dart)));
		return true;
	});

	vs.f_keypress = [&](cgogn::ui::View* view, MRMesh* selected_mesh, std::int32_t k,
						cgogn::ui::CellsSet<MRMesh, Vertex>*, cgogn::ui::CellsSet<MRMesh, Edge>*) {
		switch (k)
		{
		case GLFW_KEY_I: {

			std::vector<Volume> choix_volume;
			// std::shuffle(volume_to_simplified.begin(), volume_to_simplified.end(), g);
			int n = 0;
			for (auto it = volume_to_simplified.begin(); it != volume_to_simplified.end();)
			{
				Volume v = *it;
				choix_volume.push_back(v);
				volume_to_subdivided.push_back(v);
				it = volume_to_simplified.erase(it);
				n++;
				if (n >= 100)
					break;
			}
			// std::cout << "temps subdivided 10% volume mr : " << duration << std::endl;

			start = std::clock();
			for (Volume v : choix_volume)
			{
				cph2->disable_volume_subdivision(v, true);
			}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << duration << ";" << std::endl;
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		case GLFW_KEY_J: {
			std::vector<Volume> choix_volume;
			// std::shuffle(volume_to_subdivided.begin(), volume_to_subdivided.end(), g);
			int n = 0;
			for (auto it = volume_to_subdivided.begin(); it != volume_to_subdivided.end();)
			{

				Volume v = *it;
				choix_volume.push_back(v);
				volume_to_simplified.push_back(v);
				it = volume_to_subdivided.erase(it);

				n++;
				if (n >= 100)
					break;
			}
			// std::cout << "temps subdivided 10% volume mr : " << duration << std::endl;

			start = std::clock();
			for (Volume v : choix_volume)
			{
				auto fn = [](Vertex) {};
				cgogn::modeling::butterflySubdivisionVolume(*cph2, 0.0f, {position.get()}, {v}, fn, fn, fn);
			}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << duration << ";" << std::endl;
			vmrm.changed_connectivity(*selected_mesh, position.get());
		}
		break;
		default:
			break;
		}
	};

	return app.launch();
}
