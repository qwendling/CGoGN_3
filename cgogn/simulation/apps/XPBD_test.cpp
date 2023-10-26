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

struct BVH_Hit
{
	bool hit = false;
	Face2 face;
	Vec3 bcoords;
	double dist;
	Vec3 pos;
};

BVH_Hit intersect_bvh(Surface* surface, Attribute<Vec3>* position, acc::BVHTree<cgogn::uint32, Vec3>* surface_bvh,
					  std::vector<Face2>& surface_faces_, const acc::Ray<Vec3>& r)
{
	acc::BVHTree<cgogn::uint32, Vec3>::Hit h;
	if (surface_bvh->intersect(r, &h))
	{
		Face2 f = surface_faces_[h.idx];
		std::vector<Vertex2> vertices = incident_vertices(*surface, f);
		Vec3 p = h.bcoords[0] * cgogn::value<Vec3>(*surface, position, vertices[0]) +
				 h.bcoords[1] * cgogn::value<Vec3>(*surface, position, vertices[1]) +
				 h.bcoords[2] * cgogn::value<Vec3>(*surface, position, vertices[2]);
		return {true, f, {h.bcoords[0], h.bcoords[1], h.bcoords[2]}, h.t, p};
	}
	else
		return BVH_Hit();
}

int main(int argc, char** argv)
{
	std::string filename;
	std::string filename2;
	bool have_fine_mesh = false;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " volume_mesh_file [fine mesh] [nb_subdiv]" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	if (argc >= 3)
	{
		filename2 = std::string(argv[2]);
		have_fine_mesh = true;
	}

	int nb_subdivision = 0;
	if (argc >= 4)
	{
		nb_subdivision = std::atoi(argv[3]);
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
	cgogn::ui::FitVolumeSurface<Surface, MRMesh> fvs(app);
	cgogn::ui::Multiresolution_editing<MRMesh> mre(app);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&xp_v);
	v1->link_module(&sr);
	v1->link_module(&fvs);
	v1->link_module(&mre);

	app.init_modules();

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	cgogn::CMap2* m_fine;
	if (have_fine_mesh)
	{
		m_fine = mps.load_surface_from_file(filename2);
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
	if (have_fine_mesh)
	{
		std::shared_ptr<Attribute<Vec3>> position_surface = cgogn::get_attribute<Vec3, Vertex2>(*m_fine, "position");
		cgogn::foreach_cell(*m_fine, [&](Vertex2 v) -> bool {
			cgogn::value<Vec3>(*m_fine, position_surface, v) *= 100;
			return true;
		});
	}

	cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
		cgogn::value<Vec3>(*mrm, position, v) *= 100;
		return true;
	});

	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);
	cgogn::index_cells<Mesh::Face>(*mrm);

	// cgogn::simulation::XPBD<Mesh> xp_sim;
	// xp_sim.init_solver(*m, position);

	mrsr.set_vertex_position(*v1, *mrm, position);
	/*vmrm.subdivide(*mrm, position.get());

	vmrm.subdivide(*mrm, position.get());*/

	std::srand(124575);

	if (have_fine_mesh)
	{

		std::vector<Vertex2> surface_vertices_;
		std::vector<Face2> surface_faces_;
		std::shared_ptr<Attribute<Vec3>> position3 = cgogn::get_attribute<Vec3, Vertex2>(*m_fine, "position");
		std::shared_ptr<Attribute<Vec3>> n_vert = cgogn::add_attribute<Vec3, Vertex>(*m, "n_vert");

		fvs.set_current_surface(m_fine);
		fvs.set_current_surface_vertex_position(position3);
		fvs.set_current_volume(mrm);

		uint32 nb_vertices = mps.mesh_data(*m_fine).template nb_cells<Vertex2>();
		uint32 nb_faces = mps.mesh_data(*m_fine).template nb_cells<Face2>();

		auto vertex_index = cgogn::get_or_add_attribute<uint32, Vertex2>(*m_fine, "__bvh_vertex_index");

		std::vector<Vec3> vertex_position;
		vertex_position.reserve(nb_vertices);
		surface_vertices_.clear();
		surface_vertices_.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(*m_fine, [&](Vertex2 v) -> bool {
			cgogn::value<uint32>(*m_fine, vertex_index, v) = idx++;
			surface_vertices_.push_back(v);
			vertex_position.push_back(cgogn::value<Vec3>(*m_fine, position3, v));
			return true;
		});

		surface_faces_.clear();
		surface_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(*m_fine, [&](Face2 f) -> bool {
			surface_faces_.push_back(f);
			foreach_incident_vertex(*m_fine, f, [&](Vertex2 v) -> bool {
				face_vertex_indices.push_back(cgogn::value<uint32>(*m_fine, vertex_index, v));
				return true;
			});
			return true;
		});

		acc::BVHTree<uint32, Vec3>* surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position);

		std::vector<std::pair<Vec3, Vec3>> vec_normal_centroid;

		cgogn::foreach_cell(*m_fine, [&](Face2 f) -> bool {
			vec_normal_centroid.push_back(std::make_pair(cgogn::value<Vec3>(*m_fine, position3.get(), Vertex2(f.dart)),
														 cgogn::geometry::normal(*m_fine, f, position3.get())));

			return true;
		});
		auto fn = [&]() {
			std::vector<Vertex> vec_vertices;
			cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
				if ((cgogn::is_incident_to_boundary(*mrm, v)))
					vec_vertices.push_back(v);
				return true;
			});
#define IMAX 10
			for (int i = 1; i <= IMAX; i++)
			{
				// fvs.regularize_surface_vertices(20);
				/*fvs.optimize_volume_vertices(10.0, true);
				for (int j = 0; j < 1; j++)
					fvs.relocate_interior_vertices();*/
				std::cout << "progress : " << i << std::endl;
				for (Vertex v : vec_vertices)
				{
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
					Vec3 p = cgogn::value<Vec3>(*mrm, position.get(), v);
					// Vec3 n = cgogn::value<Vec3>(*mrm, n_vert.get(), v);
					BVH_Hit h =
						intersect_bvh(m_fine, position3.get(), surface_bvh_, surface_faces_, {p, n, 0, acc::inf});
					Vec3 pos = surface_bvh_->closest_point(p);
					double nPos = (p - pos).norm();
					intersect_bvh(m_fine, position3.get(), surface_bvh_, surface_faces_, {p, n, 0, acc::inf});
					if (h.hit)
					{
						Vec3 nFace = cgogn::geometry::normal(*m_fine, h.face, position3.get());
						Vec3 Dir = p - h.pos;
						if (Dir.dot(nFace) < 0)
						{
							double nDir = Dir.norm();
							Dir = -Dir;
							Dir.normalize();
							Dir *= std::min(nDir, nPos);
							cgogn::value<Vec3>(*mrm, position.get(), v) += (1.0 / double(IMAX + 1 - i)) * Dir;
						}
						else
						{
							h.hit = false;
						}
					}
					if (!h.hit)
					{
						n = -n;
						h = intersect_bvh(m_fine, position3.get(), surface_bvh_, surface_faces_, {p, n, 0, acc::inf});
						Vec3 pos;
						if (h.hit)
						{
							Vec3 nFace = cgogn::geometry::normal(*m_fine, h.face, position3.get());
							Vec3 Dir = p - h.pos;
							if (Dir.dot(nFace) > 0)
							{
								double nDir = Dir.norm();
								Dir = -Dir;
								Dir.normalize();
								Dir *= std::min(nDir, nPos);
								cgogn::value<Vec3>(*mrm, position.get(), v) += (double(i) / double(IMAX)) * Dir;
							}
						}
					}
				}
				for (int j = 0; j < 1; j++)
					fvs.relocate_interior_vertices();
				//  fvs.optimize_volume_vertices(10.0, true);
			}
			for (int i = 0; i < 10; i++)
			{
				// fvs.regularize_surface_vertices(20);
				// fvs.relocate_interior_vertices();
				//  fvs.optimize_volume_vertices(10.0, true);
			}
		};

		// fn();
		/*vmrm.subdivide(*mrm, position.get());
		vmrm.subdivide(*mrm, position.get());
		vmrm.subdivide(*mrm, position.get());
		mrm->current_level_ = mrm->maximum_level_;*/
		/*for (int i = 0; i < 3; i++)
			fvs.optimize_volume_vertices(10.0, true);*/
		for (int j = 0; j < nb_subdivision; j++)
		{
			/*for (int i = 0; i < 10; i++)
			{
				fvs.optimize_volume_vertices(10.0, true);
				fvs.relocate_interior_vertices();
			}*/
			vmrm.subdivide(*mrm, position.get());

			mrm->current_level_ = mrm->maximum_level_;
			vmrm.changed_connectivity(*mrm, position.get());
			fvs.update_topo();
			fn();
		}
		/*for (int i = 0; i < 10; i++)
		{
			fvs.optimize_volume_vertices(10.0, true);
			fvs.relocate_interior_vertices();
		}*/
		// fn();
		/*for (int j = 0; j < nb_subdivision; j++)
		{
			for (int i = 0; i < 10; i++)
			{
				fvs.optimize_volume_vertices(10.0, true);
				fvs.relocate_interior_vertices();
			}
			vmrm.subdivide(*mrm, position.get());

			mrm->current_level_ = mrm->maximum_level_;
			vmrm.changed_connectivity(*mrm, position.get());
			fvs.update_topo();
			for (int i = 0; i < 10; i++)
			{
				fvs.optimize_volume_vertices(10.0, true);
				fvs.relocate_interior_vertices();
			}
			fn();


		}*/

		mrm->current_level_ = 0;
	}
	else
	{
		vmrm.subdivide(*mrm, position.get());
		vmrm.subdivide(*mrm, position.get());
		// vmrm.subdivide(*mrm, position.get());
		mrm->current_level_ = 0;
	}

	cgogn::foreach_cell(*geometry_mesh, [&](Face f) -> bool {
		if (is_incident_to_boundary(*geometry_mesh, f))
		{
			geometry_mesh->activate_face_subdivision(f);
		}
		return true;
	});

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
	/*fvs.set_current_volume(geometry_mesh);
	fvs.update_topo();
	fvs.refresh_volume_skin();
	mrm->current_level_ = 0;*/

	vmrm.changed_connectivity(*mrm, position.get());
	vmrm.changed_connectivity(*geometry_mesh, position.get());

	return app.launch();
}
