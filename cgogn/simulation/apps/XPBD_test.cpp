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
#include <cgogn/simulation/ui_modules/XPBD_multiresolution.h>

#include <cgogn/simulation/algos/XPBD/XPBD.h>

// using Mesh = cgogn::CMap3;

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;
// using Mesh = cgogn::EMR_Map3_Compact;

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
	std::string filename2;
	bool have_fine_mesh = false;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " volume_mesh_file" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	if (argc == 3)
	{
		filename2 = std::string(argv[2]);
		have_fine_mesh = true;
	}

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Shape Matching");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> mrsr(app);
	cgogn::ui::VolumeSelection<Mesh> vs(app);
	cgogn::ui::XPBD_Multiresolution_View<MRMesh> xp_v(app);
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

	Mesh* m_fine;
	if (have_fine_mesh)
	{
		m_fine = mp.load_volume_from_file(filename2);
		if (!m_fine)
		{
			std::cout << "File fine could not be loaded" << std::endl;
			return 1;
		}
	}

	MRMesh* mrm = vmrm.create_mrmesh(*m, mp.mesh_name(*m));
	MRMesh* geometry_mesh = vmrm.create_mrmesh(*m, "geometry");
	geometry_mesh->parent = mrm;

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");
	std::shared_ptr<Attribute<Vec3>> normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal__anim_multires");

	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);
	cgogn::index_cells<Mesh::Face>(*mrm);

	// cgogn::simulation::XPBD<Mesh> xp_sim;
	// xp_sim.init_solver(*m, position);

	mrsr.set_vertex_position(*v1, *mrm, position);
	vmrm.subdivide(*mrm, position.get());

	vmrm.subdivide(*mrm, position.get());

	cgogn::foreach_cell(*geometry_mesh, [&](Face f) -> bool {
		if (is_incident_to_boundary(*geometry_mesh, f))
		{

			geometry_mesh->activate_face_subdivision(f);
		}
		return true;
	});

	std::srand(124575);

	if (false && have_fine_mesh)
	{
		mrm->current_level_ = 2;
		std::shared_ptr<Attribute<Vec3>> position3 = cgogn::get_attribute<Vec3, Vertex>(*m_fine, "position");
		std::shared_ptr<Attribute<Vec3>> n_vert = cgogn::add_attribute<Vec3, Vertex>(*m, "n_vert");
		std::vector<std::pair<Vec3, Vec3>> vec_normal_centroid;

		cgogn::foreach_cell(*m_fine, [&](Face f) -> bool {
			if (cgogn::is_incident_to_boundary(*m_fine, f))
			{
				vec_normal_centroid.push_back(
					std::make_pair(cgogn::value<Vec3>(*m_fine, position3.get(), Vertex(f.dart)),
								   cgogn::geometry::normal(*m_fine, f, position3.get())));
			}
			return true;
		});
		for (int i = 0; i < 100; i++)
		{
			cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
				if (!(cgogn::is_incident_to_boundary(*mrm, v)))
					return true;
				Vec3 n{0.0, 0.0, 0.0};
				cgogn::foreach_incident_face(*mrm, v, [&](Face f) -> bool {
					if (!(cgogn::is_incident_to_boundary(*mrm, f)))
						return true;
					if (cgogn::is_boundary(*mrm, f.dart))
						f.dart = cgogn::phi3(*mrm, f.dart);
					n += cgogn::geometry::normal(*mrm, f, position.get());
					return true;
				});
				n.normalize();
				cgogn::value<Vec3>(*mrm, n_vert.get(), v) = n;

				return true;
			});

			cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
				if (!(cgogn::is_incident_to_boundary(*mrm, v)))
					return true;
				double dist_min = 1000000;
				Vec3 n_min;
				double coeff;
				Vec3 p = cgogn::value<Vec3>(*mrm, position.get(), v);
				for (auto&& [c, n] : vec_normal_centroid)
				{
					double d = (c - p).norm();
					if (abs(d) < abs(dist_min))
					{
						dist_min = d;
						n_min = n;
						coeff = (c - p).dot(n);
					}
				}
				cgogn::value<Vec3>(*mrm, position.get(), v) += coeff * cgogn::value<Vec3>(*mrm, n_vert.get(), v);

				return true;
			});
		}
		mrm->current_level_ = 0;
	}

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
	vmrm.changed_connectivity(*mrm, position.get());
	vmrm.changed_connectivity(*geometry_mesh, position.get());

	return app.launch();
}
