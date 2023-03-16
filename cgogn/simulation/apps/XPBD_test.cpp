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
#include <cgogn/core/types/cmap/EMR3_compact.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/volume_selection.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/ui_modules/volume_emr_modeling.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/simulation/ui_modules/XPBD.h>

#include <cgogn/simulation/algos/XPBD/XPBD.h>

// using Mesh = cgogn::CMap3;

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;

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
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::VolumeSelection<Mesh> vs(app);
	cgogn::ui::XPBD_View<Mesh> xp_v(app);
	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&xp_v);

	app.init_modules();

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	MRMesh* mrm = vmrm.create_mrmesh(*m, mp.mesh_name(*m));

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");

	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);
	cgogn::index_cells<Mesh::Face>(*mrm);

	// cgogn::simulation::XPBD<Mesh> xp_sim;
	// xp_sim.init_solver(*m, position);

	mrsr.set_vertex_position(*v1, *mrm, position);
	vmrm.subdivide(*mrm, position.get());

	vmrm.subdivide(*mrm, position.get());

	std::srand(std::time(nullptr));

	/*cgogn::foreach_cell(*m, [&](Vertex v) -> bool {
		cgogn::value<Vec3>(*m, position, v) +=
			Vec3(std::rand() / double(RAND_MAX + 1u), std::rand() / double(RAND_MAX + 1u),
				 std::rand() / double(RAND_MAX + 1u));
		return true;
	});*/
	// xp_sim.solver(*m, 0.01f);
	/*std::clock_t start = std::clock();
	double duration = 0;
	for (int i = 0; i < 10; i++)
		xp_sim.solver(*m, 0.01f);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps solve xpbd : " << duration / 100.0f << std::endl;*/

	// mp.emit_attribute_changed(*m, position.get());

	return app.launch();
}
