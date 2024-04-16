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
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/volume_selection.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/ui_modules/Fit_Volume_To_Surface.h>
#include <cgogn/modeling/ui_modules/multiresolution_editing.h>
#include <cgogn/modeling/ui_modules/volume_emr_modeling.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/simulation/ui_modules/XPBD.h>
#include <cgogn/simulation/ui_modules/XPBD_multiresolution.h>

#include <cgogn/simulation/algos/XPBD/XPBD.h>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

// using Mesh = cgogn::CMap3;

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;
// using Mesh = cgogn::EMR_Map3_Compact;
using Surface = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

using Face2 = typename cgogn::mesh_traits<Surface>::Face;
using Vertex2 = typename cgogn::mesh_traits<Surface>::Vertex;

using Vec3 = cgogn::geometry::Vec3;
using uint32 = cgogn::uint32;


int main(int argc, char** argv)
{

	if (argc > 2)
	{
		std::cout << "Usage: " << argv[0] << " [nb_subdiv]" << std::endl;
		return 1;
	}



	int nb_subdivision = 2;
	if (argc == 2)
	{
		nb_subdivision = std::atoi(argv[1]);
	}

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("XPBD");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::MeshProvider<cgogn::CMap2> mps(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::VolumeSelection<Mesh> vs(app);
	cgogn::ui::XPBD_Multiresolution_View<MRMesh> xp_v(app);
	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);
	cgogn::ui::Multiresolution_editing<MRMesh> mre(app);
	cgogn::ui::SurfaceDifferentialProperties<cgogn::CMap2> sdp(app);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&xp_v);
	v1->link_module(&sr);
	v1->link_module(&mre);

	app.init_modules();

	Mesh* m ;

	Volume vol = cgogn::add_prism(static_cast<cgogn::CMap2&>(*m), 4u, false);



	MRMesh* mrm = vmrm.create_mrmesh(*m, "mecanique");
	MRMesh* topo = vmrm.create_mrmesh(*m, "Topology");
	mrm->topology_ = topo;
	MRMesh* geometry_mesh = vmrm.create_mrmesh(*m, "geometry");
	geometry_mesh->parent = mrm;

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");
	std::shared_ptr<Attribute<Vec3>> normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal__anim_multires");


	cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
		cgogn::value<Vec3>(*mrm, position, v) *= 100;
		return true;
	});

	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);
	cgogn::index_cells<Mesh::Face>(*mrm);



	mrsr.set_vertex_position(*v1, *mrm, position);

	std::srand(124575);

	for (int i = 0; i < nb_subdivision; i++)
	{
		vmrm.subdivide(*mrm, position.get());
	}
	mrm->current_level_ = 0;


	 cgogn::foreach_cell(*geometry_mesh, [&](Face f) -> bool {

		if (is_incident_to_boundary(*geometry_mesh, f))
		{
			geometry_mesh->activate_face_subdivision(f);
		}
		return true;
	});
	mrm->current_level_ = 0;

	vmrm.changed_connectivity(*mrm, position.get());
	vmrm.changed_connectivity(*geometry_mesh, position.get());

	return app.launch();
}
