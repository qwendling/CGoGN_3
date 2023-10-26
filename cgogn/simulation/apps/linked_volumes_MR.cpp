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
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/volume_selection.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/ui_modules/linked_volumes.h>
#include <cgogn/modeling/ui_modules/volume_emr_modeling.h>
#include <cgogn/modeling/ui_modules/volume_mr_modeling.h>
#include <cgogn/rendering/ui_modules/topo_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/simulation/ui_modules/shape_matching.h>

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;
using EMR_Map3 = cgogn::EMR_Map3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<MRMesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<MRMesh>::Vertex;
using Edge = typename cgogn::mesh_traits<MRMesh>::Edge;
using Face = typename cgogn::mesh_traits<MRMesh>::Face;
using Volume = typename cgogn::mesh_traits<MRMesh>::Volume;

using Vec3 = cgogn::geometry::Vec3;

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "fr-FR");
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
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);
	cgogn::ui::ShapeMatching<MRMesh> sm(app);
	cgogn::ui::LinkedVolumes<MRMesh> lv(app);
	cgogn::ui::TopoRender<MRMesh> tr(app);

	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&sm);
	v1->link_module(&lv);
	v1->link_module(&tr);

	/*cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mp);
	v2->link_module(&mrsr);
	v2->link_module(&vs);
	v2->link_module(&sm);
	v2->link_module(&lv);*/

	app.init_modules();

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	MRMesh* mrm = vmrm.create_mrmesh(*m, mp.mesh_name(*m));
	MRMesh* mrm2 = vmrm.create_mrmesh(*m, mp.mesh_name(*m));
	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");

	cgogn::index_cells<MRMesh::Volume>(*mrm);
	cgogn::index_cells<MRMesh::Edge>(*mrm);
	cgogn::index_cells<MRMesh::Face>(*mrm);

	vmrm.subdivide(*mrm, position.get());

	vmrm.subdivide(*mrm, position.get());

	mrm2->parent = mrm;

	mrsr.set_vertex_position(*v1, *mrm, position);
	mrsr.set_vertex_position(*v1, *mrm2, nullptr);
	v1->scene_bb_locked_ = true;

	/*mrsr.set_vertex_position(*v2, *mrm, nullptr);
	mrsr.set_vertex_position(*v2, *mrm2, position);
	v2->scene_bb_locked_ = true;*/

	vmrm.changed_connectivity(*mrm2, position.get());
	std::vector<Volume> list_cut_volumes;

	std::srand(164512792);
	/*while (std::rand() / ((RAND_MAX + 1u) / 5) > 1)
	{
		cgogn::foreach_cell(*mrm2, [&list_cut_volumes](Volume v) -> bool {
			list_cut_volumes.push_back(v);
			return true;
		});
		for (Volume v : list_cut_volumes)
		{

			int tmp = std::rand() / ((RAND_MAX + 1u) / 2);
			if (tmp == 1)
			{
				mrm2->activate_volume_subdivision(v);
			}
		}
		list_cut_volumes.clear();
	}*/
	vmrm.changed_connectivity(*mrm2, position.get());

	return app.launch();
}
