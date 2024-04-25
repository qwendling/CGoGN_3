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
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/core/types/cmap/EMR3_compact.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/algos/centroid.h>
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

#include <chrono>
#include <random>

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

std::random_device rd;

class LocalInterface : public cgogn::ui::ViewModule
{

public:
	LocalInterface(const cgogn::ui::App& app)
		: cgogn::ui::ViewModule(app, "LocalInterface"), mesh_(nullptr), vertex_position_(nullptr),
		  mesh_provider_(nullptr), vol_render_(nullptr), moving_color_(1.0f, 0.0f, 1.0f, 1.0f), cm_not_cut(nullptr)
	{
		view_ = app.current_view();
	}

	~LocalInterface()
	{
	}
	void force_update()
	{
		for (cgogn::ui::View* v : linked_views_)
			v->request_update();
	}

	double dist_plan(Vec3 n, Vec3 p)
	{
		return n.dot(p);
	}

	void apply_random_cut()
	{
		bool is_running = xmv->running_;
		xmv->stop();
		Vertex v_select;
		int nb_boucle = 0;

		std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

		std::vector<Vertex> shuffle_vertices;
		cgogn::foreach_cell(*mesh_, [&](Vertex v) -> bool {
			if (mesh_->dart_level(v.dart) != 0)
				return true;
			if (cm_not_cut->is_marked(v))
				return true;
			shuffle_vertices.push_back(v);

			return true;
		});

		if (!shuffle_vertices.empty())
		{
			std::shuffle(shuffle_vertices.begin(), shuffle_vertices.end(), gen);
			v_select = shuffle_vertices.front();
		}

		if (!v_select.dart.is_nil())
		{

			Vec3 pos = cgogn::value<Vec3>(*mesh_, vertex_position_.get(), v_select);
			Vec3 pos_voisin = cgogn::value<Vec3>(*mesh_, vertex_position_.get(), Vertex(phi2(*mesh_, v_select.dart)));

			/*double taille_sphere = (pos - pos_voisin).norm() * 2.;
			Vec3 centre_sphere =
				Vec3(double(rand()) / double(RAND_MAX) - 0.5, double(rand()) / double(RAND_MAX) - 0.5,
					 double(rand()) / double(RAND_MAX) - 0.5)
						.normalized() *
					taille_sphere +
				pos;*/

			std::uniform_real_distribution<> dis(-1.0, 1.0);
			std::uniform_real_distribution<> dis2(0.2, 1.0);

			Vec3 a(dis(gen), dis(gen), 0.);
			a.normalize();
			double d = pos.dot(a);
			mesh_->start_writer();
			std::vector<std::shared_ptr<Attribute<Vec3>>> list_update_attribute;
			std::cout << "Début découpe" << std::endl;

			cgogn::CellMarker<MRMesh, Volume> cm_cut(*mesh_);

			std::vector<Volume> list_volume;
			std::vector<Volume> list_volume_topo;
			std::vector<Volume> list_volume_not_cut;

			cgogn::foreach_incident_volume(*mesh_->topology_, v_select, [&](Volume v) -> bool {
				list_volume_topo.push_back(v);
				std::vector<Volume> tmp = simu_solver->get_all_current_child(*mesh_->topology_, v);
				for (Volume w : tmp)
				{
					list_volume.push_back(w);
					cm_cut.mark(w);
				}
				cgogn::foreach_adjacent_volume_through_vertex(*mesh_->topology_, v, [&](Volume w) -> bool {
					list_volume_not_cut.push_back(w);
					return true;
				});
				return true;
			});

			auto fn_activate_volume = [&](MRMesh* m, std::vector<Volume> list_volume_, bool simu_activate) {
				std::vector<Volume> list_volume_activate;
				do
				{
					list_volume_activate.clear();
					for (Volume v : list_volume_)
					{
						cm_cut.mark(v);
						if (m->volume_level(v.dart) == m->maximum_level_)
							continue;
						Vec3 cm = cgogn::geometry::centroid<Vec3>(*m, v, vertex_position_.get());
						double dist_plan_volume = cm.dot(a) - d;
						// double dist_sphere_volume = (centre_sphere - cm).norm() - taille_sphere;
						cgogn::foreach_incident_vertex(*m, v, [&](Vertex w) -> bool {
							if ((cgogn::value<Vec3>(*m, vertex_position_.get(), w).dot(a) - d) * dist_plan_volume < 0)
							{
								list_volume_activate.push_back(v);
								return false;
							}
							/*if (((cgogn::value<Vec3>(*mesh_, vertex_position_.get(), w) - centre_sphere).norm() -
								 taille_sphere) *
									dist_sphere_volume <
								0)
							{
								list_volume_activate.push_back(v);
								return false;
							}*/
							return true;
						});
					}

					list_volume_.clear();
					for (Volume v : list_volume_activate)
					{
						cgogn::foreach_incident_vertex(*m, v, [&](Vertex w) -> bool {
							list_volume_.push_back(Volume(w.dart));
							return true;
						});
					}
					if (simu_activate)
						simu_solver->activate_volume_tree(*m, list_volume_activate);
					else
						for (Volume v : list_volume_activate)
						{
							m->activate_volume_subdivision(v);
						}
				} while (!list_volume_activate.empty());
			};

			fn_activate_volume(mesh_, list_volume, true);
			// fn_activate_volume(mesh_->topology_, list_volume_topo, false);

			/*std::vector<Volume> list_volume_cut;

			cgogn::foreach_cell(*mesh_, [&](Volume v) -> bool {
				if (cm_cut.is_marked(v))
					list_volume_cut.push_back(v);
				return true;
			});*/

			std::vector<Face> face_unsew;

			parallel_foreach_cell(*mesh_, [&](Face f) -> bool {
				if (cgogn::is_incident_to_boundary(*mesh_, f) || !cm_cut.is_marked(Volume(f.dart)) ||
					!cm_cut.is_marked(Volume(cgogn::phi3(*mesh_, f.dart))))
					return true;

				Volume v1(f.dart), v2(cgogn::phi3(*mesh_, f.dart));

				Vec3 cm1 = cgogn::geometry::centroid<Vec3>(*mesh_, v1, vertex_position_.get());
				Vec3 cm2 = cgogn::geometry::centroid<Vec3>(*mesh_, v2, vertex_position_.get());

				double d1 = cm1.dot(a) - d;
				double d2 = cm2.dot(a) - d;

				/*double d1_sphere = (centre_sphere - cm1).norm() - taille_sphere;
				double d2_sphere = (centre_sphere - cm2).norm() - taille_sphere;*/

				if (d1 * d2 < 0)
					face_unsew.push_back(f);
				/*if (d1_sphere * d2_sphere < 0)
					face_unsew.push_back(f);*/

				return true;
			});

			if (!face_unsew.empty())
			{

				for (Face f : face_unsew)
				{
					unsew_volume(
						*mesh_, f,
						[&](std::pair<Vertex, Vertex> p) -> bool {
							list_update_attribute.push_back(simu_solver->pos_);
							list_update_attribute.push_back(simu_solver->init_pos_);
							list_update_attribute.push_back(simu_solver->speed_);
							for (auto attr : list_update_attribute)
							{
								cgogn::value<Vec3>(*mesh_, attr, p.second) = cgogn::value<Vec3>(*mesh_, attr, p.first);
							}
							return true;
						},
						true);
				}
				simu_solver->update_topo(*mesh_);
			}

			std::cout << "Fin découpe" << std::endl;

			for (Volume v : list_volume_not_cut)
			{

				cgogn::foreach_incident_vertex(*mesh_->topology_, v, [&](Vertex v2) -> bool {
					cm_not_cut->mark(v2);
					return true;
				});
			}

			/*list_update_attribute.push_back(vertex_position_);

			for (auto attr : list_update_attribute)
			{
				mesh_provider_->emit_attribute_changed(*mesh_, attr.get());
			}
			mesh_provider_->emit_connectivity_changed(*mesh_);
			// mesh_provider_->emit_connectivity_changed(*selected_mesh_->topology_);
			xmv->refresh_volume_skin();
			xmv->surface_provider_->emit_attribute_changed(*xmv->volume_skin_,
			xmv->volume_skin_vertex_position_.get());*/
			// sdp->update_normal();
			mesh_->end_writer();

			if (is_running)
				xmv->start();
		}
	}

	void adapt_random()
	{
		bool is_running = xmv->running_;
		xmv->stop();
		mesh_->start_writer();

		std::vector<std::shared_ptr<Attribute<Vec3>>> list_update_attribute;
		std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(0.0, 1.0);

		xmv->simu_solver.Update_error_quotat(*mesh_, 7., [&](MRMesh&, Volume) -> double { return dis(gen); });

		list_update_attribute.push_back(vertex_position_);
		for (auto attr : list_update_attribute)
		{
			mesh_provider_->emit_attribute_changed(*mesh_, attr.get());
		}
		mesh_provider_->emit_connectivity_changed(*mesh_);
		mesh_->end_writer();
		if (is_running)
			xmv->start();
	}

	void init() override
	{
		mesh_provider_ = static_cast<cgogn::ui::MeshProvider<MRMesh>*>(
			app_.module("MeshProvider (" + std::string{cgogn::mesh_traits<MRMesh>::name} + ")"));

		vol_render_ = static_cast<cgogn::ui::VolumeRender<MRMesh>*>(
			app_.module("VolumeRender (" + std::string{cgogn::mesh_traits<MRMesh>::name} + ")"));

		shape_ = cgogn::rendering::ShapeDrawer::instance();
		shape_->color(cgogn::rendering::ShapeDrawer::CYLINDER) = cgogn::rendering::GLColor(0.5294, 0.6078, 0.6078, 1);
		shape_->color(cgogn::rendering::ShapeDrawer::CUBE) = cgogn::rendering::GLColor(0.5294, 0.6078, 0.6078, 1);
		shape_->color(cgogn::rendering::ShapeDrawer::SPHERE) = cgogn::rendering::GLColor(0.88, 0.55, 0.17, 1);
	}

	void draw(cgogn::ui::View* view) override
	{
		using namespace cgogn;
		const rendering::GLMat4& proj_matrix = view->projection_matrix();
		const rendering::GLMat4& view_matrix = view->modelview_matrix();

		// Eigen::Affine3f transfo = Eigen::Translation3f(Eigen::Vector3f(0., 0., 100.)) * Eigen::Scaling(500.f);
		// shape_->draw(rendering::ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());
	}

	void left_panel() override
	{
		if (ImGui::Button("Perform random cut"))
		{
			apply_random_cut();
		}

		if (ImGui::Button("Adapt max random"))
		{
			adapt_random();
		}

		if (ImGui::Button("Start anim"))
		{
			xmv->take_screenshot_ = true;
			cgogn::launch_thread([this]() {
				xmv->frame_number_ = 0;
				for (int i = 0; i < 10; i++)
				{
					for (int j = 0; j < 25; j++)
					{
						while (xmv->need_update_)
						{
							std::this_thread::sleep_for(std::chrono::milliseconds(5));
						}
						mesh_->start_writer();
						xmv->step();
						mesh_->end_writer();
						// xmv->selected_view_->request_update();
					}
					if (i % 2 == 0)
					{
						mesh_->start_writer();
						std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
						std::uniform_real_distribution<> dis(0.0, 1.0);
						xmv->simu_solver.Update_error_quotat(*mesh_, 0.5,
															 [&](MRMesh&, Volume) -> double { return dis(gen); });
						// xmv->step();
						xmv->simu_solver.Update_error_quotat(*mesh_, 0.5,
															 [&](MRMesh&, Volume) -> double { return dis(gen); });
						// xmv->step();
						mesh_->end_writer();

						apply_random_cut();
						mesh_->start_writer();
						xmv->step();
						xmv->simu_solver.Update_error_quotat(*mesh_, 7.,
															 [&](MRMesh&, Volume) -> double { return dis(gen); });
						mesh_->end_writer();
					}
					else
					{
						mesh_->start_writer();
						// xmv->step();

						std::vector<std::shared_ptr<Attribute<Vec3>>> list_update_attribute;
						std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
						std::uniform_real_distribution<> dis(0.0, 1.0);

						xmv->simu_solver.Update_error_quotat(*mesh_, 7.,
															 [&](MRMesh&, Volume) -> double { return dis(gen); });
						mesh_->end_writer();
					}
				}
			});
		}

		if (ImGui::Button("Fix border"))
		{
			uint32 cur = mesh_->current_level_;
			mesh_->current_level_ = mesh_->maximum_level_;

			double x_min = DBL_MAX, x_max = -1e10, y_min = DBL_MAX, y_max = -1e10;
			cgogn::foreach_cell(*mesh_, [&](Vertex v) -> bool {
				const Vec3& p = cgogn::value<Vec3>(*mesh_, vertex_position_.get(), v);
				if (p.x() < x_min)
					x_min = p.x();
				if (p.y() < y_min)
					y_min = p.y();

				if (p.x() > x_max)
					x_max = p.x();
				if (p.y() > y_max)
					y_max = p.y();
				return true;
			});

			auto fixed_vertex = cgogn::get_attribute<bool, Vertex>(*mesh_, "fixed_vertex");
			double delta = 10.1;
			cgogn::foreach_cell(*mesh_, [&](Vertex v) -> bool {
				const Vec3& p = cgogn::value<Vec3>(*mesh_, vertex_position_.get(), v);
				if (p.x() < x_min + delta || p.y() < y_min + delta || p.x() > x_max - delta || p.y() > y_max - delta)
				{
					cgogn::value<bool>(*mesh_, fixed_vertex.get(), v) = true;
					cgogn::foreach_adjacent_vertex_through_edge(*mesh_, v, [&](Vertex w) -> bool {
						cm_not_cut->mark(w);
						return true;
					});
				}

				return true;
			});
			mesh_->current_level_ = cur;

			cgogn::foreach_cell(*mesh_, [&](Vertex v) -> bool {
				const Vec3& p = cgogn::value<Vec3>(*mesh_, vertex_position_.get(), v);
				if (p.x() < x_min + delta || p.y() < y_min + delta || p.x() > x_max - delta || p.y() > y_max - delta)
				{
					cgogn::foreach_adjacent_vertex_through_edge(*mesh_, v, [&](Vertex w) -> bool {
						cm_not_cut->mark(w);
						return true;
					});
				}

				return true;
			});
		}

		if (ImGui::Button("Load mesh"))
		{
			Vec3 cm = cgogn::geometry::centroid<Vec3>(*mesh_, vertex_position_.get());
			cgogn::foreach_cell(*mesh_, [&](Vertex v) -> bool {
				Vec3& p = cgogn::value<Vec3>(*mesh_, vertex_position_.get(), v);
				p = (p - cm) * 1.1 + cm;
				return true;
			});
		}

		if (ImGui::Button("move mesh"))
		{
			cgogn::CellMarkerStore<MRMesh, Vertex> vm(*mesh_);
			std::vector<Vertex> CC_0;
			CC_0.push_back(Vertex(cgogn::Dart(0)));
			vm.mark(Vertex(cgogn::Dart(0)));
			cgogn::value<Vec3>(*mesh_, vertex_position_.get(), Vertex(cgogn::Dart(0))) += Vec3(3, 0, 0);
			int cur = mesh_->current_level_;
			mesh_->current_level_ = mesh_->maximum_level_;
			while (!CC_0.empty())
			{
				Vertex v = CC_0.back();
				CC_0.pop_back();
				cgogn::foreach_adjacent_vertex_through_edge(*mesh_, v, [&](Vertex w) -> bool {
					if (!vm.is_marked(w))
					{
						vm.mark(w);
						CC_0.push_back(w);
						cgogn::value<Vec3>(*mesh_, vertex_position_.get(), w) += Vec3(3, 0, 0);
					}
					return true;
				});
			}
			mesh_->current_level_ = cur;

			mesh_provider_->emit_attribute_changed(*mesh_, vertex_position_.get());
			mesh_provider_->emit_connectivity_changed(*geom_);

			force_update();
		}
		if (!moving_dart_.is_nil())
			ImGui::Text("Dart index : %d", moving_dart_.index);
	}
	MRMesh* mesh_;
	MRMesh* mesh_meca_;
	MRMesh* geom_;
	cgogn::ui::View* view_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	cgogn::ui::MeshProvider<MRMesh>* mesh_provider_;
	cgogn::ui::VolumeRender<MRMesh>* vol_render_;
	cgogn::Dart d_hexa_;
	cgogn::Dart d_pyra_;
	cgogn::Dart moving_dart_;
	Eigen::Vector4f moving_color_;
	cgogn::CellMarker<MRMesh, Vertex>* cm_not_cut;
	cgogn::simulation::XPBD_Multiresolution* simu_solver;
	cgogn::ui::XPBD_Multiresolution_View<MRMesh>* xmv;
	cgogn::ui::SurfaceDifferentialProperties<cgogn::CMap2>* sdp;
	float expl_vol_;
	cgogn::rendering::ShapeDrawer* shape_;
};

int main(int argc, char** argv)
{

	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [nb_subdiv]" << std::endl;
		return 1;
	}
	std::string filename = argv[1];

	int nb_subdivision = 2;
	if (argc == 3)
	{
		nb_subdivision = std::atoi(argv[2]);
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
	cgogn::ui::FitVolumeSurface<Surface, MRMesh> fvs(app);
	LocalInterface interf(app);

	xp_v.sdp = &sdp;

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrsr);
	v1->link_module(&vs);
	v1->link_module(&xp_v);
	v1->link_module(&sr);
	v1->link_module(&mre);
	v1->link_module(&interf);
	v1->link_module(&fvs);

	app.init_modules();

	Mesh* m = mp.load_volume_from_file(filename);

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

	interf.mesh_ = mrm;
	interf.mesh_meca_ = mrm;
	interf.vertex_position_ = position;
	interf.geom_ = geometry_mesh;
	interf.cm_not_cut = new cgogn::CellMarker<MRMesh, Vertex>(*mrm);
	interf.simu_solver = &xp_v.simu_solver;
	interf.xmv = &xp_v;
	interf.sdp = &sdp;

	mrsr.set_vertex_position(*v1, *mrm, position);

	std::srand(17101995);

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

	fvs.set_current_volume(geometry_mesh);
	fvs.update_topo();
	fvs.refresh_volume_skin();

	mrsr.set_vertex_position(*v1, *geometry_mesh, nullptr);
	mrsr.set_vertex_position(*v1, *geometry_mesh->topology_, nullptr);
	mrsr.parameters_[mrsr.selected_view_][mrm].render_volumes_ = false;
	mrsr.parameters_[mrsr.selected_view_][mrm].render_vertices_ = true;
	mrsr.parameters_[mrsr.selected_view_][mrm].vertex_scale_factor_ = 0.3;

	v1->update_scene_bb();
	v1->request_update();

	return app.launch();
}
