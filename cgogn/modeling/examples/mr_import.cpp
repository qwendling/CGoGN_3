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
#include <cgogn/core/types/cmap/phi.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/animation_multiresolution/animation_multiresolution.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_emr_modeling/volume_emr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;
using EMR_Map3 = cgogn::EMR_Map3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<MRMesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<MRMesh>::Vertex;
using Edge = typename cgogn::mesh_traits<MRMesh>::Edge;
using Face = typename cgogn::mesh_traits<MRMesh>::Face;
using Volume = typename cgogn::mesh_traits<MRMesh>::Volume;
using Dart = cgogn::Dart;

using Vec3 = cgogn::geometry::Vec3;

int main(int argc, char** argv)
{
	std::string filename, filename2, filename3;
	if (argc < 4)
	{
		std::cout << "Usage: " << argv[0] << " volume_mesh_file0 volume_mesh_file1 volume_mesh_file2" << std::endl;
		return 1;
	}
	else
	{
		filename = std::string(argv[1]);
		filename2 = std::string(argv[2]);
		filename3 = std::string(argv[3]);
	}

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("MR Import");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	// cgogn::ui::VolumeRender<Mesh> vr(app);
	cgogn::ui::VolumeRender<MRMesh> vrmr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);

	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrmp);
	// v1->link_module(&vr);
	v1->link_module(&vrmr);
	v1->link_module(&vs);

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	Mesh* m2 = mp.load_volume_from_file(filename2);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	Mesh* m3 = mp.load_volume_from_file(filename3);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	std::shared_ptr<Attribute<Vec3>> position1 = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> position2 = cgogn::get_attribute<Vec3, Vertex>(*m2, "position");
	std::shared_ptr<Attribute<Vec3>> position3 = cgogn::get_attribute<Vec3, Vertex>(*m3, "position");

	MRMesh* mrm = vmrm.create_mrmesh(*m, "meca");
	MRMesh* geometry_mesh = vmrm.create_mrmesh(*m, "geometry");
	MRMesh* Visu_mesh = vmrm.create_mrmesh(*m, "Visu");

	vmrm.selected_vertex_parents_ = cgogn::add_attribute<std::array<Vertex, 4>, Vertex>(*m, "parents");
	vmrm.selected_vertex_relative_position_ = cgogn::add_attribute<Vec3, Vertex>(*m, "relative_position");

	vs.selected_mesh_ = mrm;
	cgogn::index_cells<Mesh::Face>(*mrm);
	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);
	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");
	std::shared_ptr<Attribute<Vec3>> n_vert = cgogn::add_attribute<Vec3, Vertex>(*mrm, "n_vert");

	std::vector<std::pair<Vec3, Vec3>> vec_normal_centroid;

	cgogn::foreach_cell(*m3, [&](Face f) -> bool {
		if (cgogn::is_incident_to_boundary(*m3, f))
		{
			vec_normal_centroid.push_back(std::make_pair(cgogn::value<Vec3>(*m3, position3.get(), Vertex(f.dart)),
														 cgogn::geometry::normal(*m3, f, position3.get())));
		}
		return true;
	});

	auto fn = [&]() {
		for (int i = 0; i < 100; i++)
		{
			cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
				if (mrm->dart_level(v.dart) != mrm->maximum_level_)
				{
					return true;
				}
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
				if (mrm->dart_level(v.dart) != mrm->maximum_level_)
				{
					return true;
				}
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
		cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
			if (mrm->dart_level(v.dart) == 0)
			{
				return true;
			}
			std::array<Vertex, 4> p = cgogn::value<std::array<Vertex, 4>>(*mrm, vmrm.selected_vertex_parents_, v);
			Vec3 A = cgogn::value<Vec3>(*mrm, position, p[0]);
			Vec3 B = cgogn::value<Vec3>(*mrm, position, p[1]);
			Vec3 C = cgogn::value<Vec3>(*mrm, position, p[2]);
			Vec3 D = cgogn::value<Vec3>(*mrm, position, p[3]);

			/*Vec3 X = (B - A).normalized();
			Vec3 Y = (C - A).normalized();
			Vec3 Z = (D - A).normalized();*/

			Vec3 X = (B - A);
			Vec3 Y = (C - A);
			Vec3 Z = (D - A);
			Eigen::Matrix3d mb;

			mb.col(0) = X;
			mb.col(1) = Y;
			mb.col(2) = Z;
			Eigen::Matrix3d inv;
			bool is_inversible;
			mb.computeInverseWithCheck(inv, is_inversible);
			assert(is_inversible);
			if (!is_inversible)
				std::cout << "pb" << std::endl;
			cgogn::value<Vec3>(*mrm, vmrm.selected_vertex_relative_position_, v) =
				inv * (cgogn::value<Vec3>(*mrm, position, v) - A);
			return true;
		});
	};

	fn();

	m->add_resolution();
	mrm->change_resolution_level(1);
	vmrm.subdivide(*mrm, position.get());

	fn();

	m->add_resolution();
	mrm->change_resolution_level(2);
	vmrm.subdivide(*mrm, position.get());

	fn();

	auto md = mrmp.mesh_data(mrm);
	md->template add_cells_set<Edge>();

	mrmp.set_mesh_bb_vertex_position(mrm, position);

	/*vr.set_vertex_position(*v1, *m, nullptr);
	vr.set_vertex_position(*v1, *m2, nullptr);
	vr.set_vertex_position(*v1, *m3, nullptr);*/

	// std::srand(std::time(nullptr));
	std::srand(2124512438);

	cgogn::foreach_cell(*geometry_mesh, [&](Face f) -> bool {
		if (is_incident_to_boundary(*geometry_mesh, f))
		{

			geometry_mesh->activate_face_subdivision(f);
		}
		return true;
	});

	Visu_mesh->current_level_ = Visu_mesh->maximum_level_;

	mrm->current_level_ = 0;
	std::vector<Volume> vol_vec;
	std::vector<Volume> vol_vec_simpl;
	cgogn::CellMarker<MRMesh, Volume> vm(*mrm);
	vmrm.changed_connectivity(*mrm, position.get());
	vmrm.changed_connectivity(*Visu_mesh, position.get());
	vmrm.changed_connectivity(*geometry_mesh, position.get());

	vs.f_keypress = [&](cgogn::ui::View*, MRMesh*, std::int32_t k, cgogn::ui::CellsSet<MRMesh, Vertex>*,
						cgogn::ui::CellsSet<MRMesh, Edge>*) {
		switch (k)
		{
		case GLFW_KEY_I: {
			mrm->current_level_ = 0;

			vmrm.changed_connectivity(*mrm, position.get());
		}
		break;
		case GLFW_KEY_U:
			vmrm.changed_connectivity(*mrm, position.get());
			break;
		case GLFW_KEY_T: {
			auto fn = [&]() {
				if (vol_vec.empty())
				{
					vm.unmark_all();
					Volume selected_volume;
					while (!selected_volume.is_valid())
					{
						cgogn::foreach_cell(*mrm, [&](Volume v) {
							if (rand() % 100 < 3)
							{
								selected_volume = v;
								return false;
							}
							return true;
						});
					}
					vol_vec.push_back(selected_volume);
					vm.mark(selected_volume);
				}
				for (Volume v : vol_vec_simpl)
				{
					mrm->disable_volume_subdivision(v, true);
				}
				vol_vec_simpl.clear();
				std::vector<Volume> vol_vec_tmp;
				for (Volume v : vol_vec)
				{
					cgogn::foreach_incident_vertex(*mrm, v, [&](Vertex w) -> bool {
						cgogn::foreach_incident_volume(*mrm, w, [&](Volume v2) -> bool {
							cgogn::foreach_incident_vertex(*mrm, v2, [&](Vertex w2) -> bool {
								cgogn::foreach_incident_volume(*mrm, w2, [&](Volume v3) -> bool {
									if (vm.is_marked(v3))
									{
										return true;
									}
									vm.mark(v3);
									vol_vec_tmp.push_back(v3);
									return true;
								});
								return true;
							});
							return true;
						});
						return true;
					});
					vol_vec_simpl.push_back(v);
				}

				for (Volume v : vol_vec)
				{
					mrm->activate_volume_subdivision(v);
				}
				vol_vec.swap(vol_vec_tmp);
				vol_vec_tmp.clear();
				vmrm.changed_connectivity(*mrm, position.get());
			};
			fn();
		}

		break;
		}
	};

	return app.launch();
}
