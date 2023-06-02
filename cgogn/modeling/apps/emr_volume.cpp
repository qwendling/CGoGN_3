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
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/volume_selection.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/ui_modules/volume_emr_modeling.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/topo_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/simulation/ui_modules/animation_multiresolution.h>

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

class LocalInterface : public cgogn::ui::ViewModule
{

public:
	LocalInterface(const cgogn::ui::App& app)
		: cgogn::ui::ViewModule(app, "LocalInterface"), mesh_(nullptr), vertex_position_(nullptr),
		  mesh_provider_(nullptr), vol_render_(nullptr), topo_render_(nullptr), moving_color_(1.0f, 0.0f, 1.0f, 1.0f)
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

	Eigen::Vector4f get_color(Dart d)
	{
		Eigen::Vector4f result;
		switch (mesh_->dart_level(d))
		{
		case 0:
			result = Eigen::Vector4f(1.0f, 0.0, 0.0f, 1.0f);
			break;
		case 1:
			result = Eigen::Vector4f(0.0f, 1.0, 0.0f, 1.0f);
			break;
		case 2:
			result = Eigen::Vector4f(0.0f, 0.0, 1.0f, 1.0f);
			break;
		case 3:
			result = Eigen::Vector4f(1.0f, 1.0, 0.0f, 1.0f);
			break;
		default:
			result = Eigen::Vector4f(1.0f, 0.0, 1.0f, 1.0f);
		}
		return result;
	}

	void init() override
	{
		mesh_provider_ = static_cast<cgogn::ui::MeshProvider<MRMesh>*>(
			app_.module("MeshProvider (" + std::string{cgogn::mesh_traits<MRMesh>::name} + ")"));

		vol_render_ = static_cast<cgogn::ui::VolumeRender<MRMesh>*>(
			app_.module("VolumeRender (" + std::string{cgogn::mesh_traits<MRMesh>::name} + ")"));
		topo_render_ = static_cast<cgogn::ui::TopoRender<MRMesh>*>(
			app_.module("TopoRender (" + std::string{cgogn::mesh_traits<MRMesh>::name} + ")"));
	}

	void left_panel() override
	{
		if (ImGui::SliderFloat("Explode", &expl_vol_, 0.01f, 1.0f))
		{
			vol_render_->set_volume_explode(*view_, *mesh_, expl_vol_);
			topo_render_->set_volume_explode(expl_vol_ + 0.02f);
			force_update();
		}

		if (ImGui::Button("init moving"))
		{

			for (Dart d : moving_darts)
			{
				topo_render_->reset_dart_color(d);
			}
			moving_darts.clear();
			cgogn::foreach_cell(*mesh_, [&](Face f) -> bool {
				if (cgogn::is_boundary(*mesh_, f.dart))
					f.dart = phi3(*mesh_, f.dart);
				if (!cgogn::is_boundary(*mesh_, phi3(*mesh_, f.dart)))
				{
					moving_darts.push_back(phi3(*mesh_, f.dart));
					topo_render_->set_dart_color(phi3(*mesh_, f.dart), get_color(phi3(*mesh_, f.dart)));
				}
				moving_darts.push_back(f.dart);
				topo_render_->set_dart_color(f.dart, get_color(f.dart));
				return true;
			});
			force_update();
		}

		if (ImGui::Button("phi 1"))
		{
			for (Dart& d : moving_darts)
			{
				topo_render_->reset_dart_color(d);
				d = phi1(*mesh_, d);
				topo_render_->set_dart_color(d, get_color(d));
			}
			force_update();
		}
		if (ImGui::Button("phi -1"))
		{
			for (Dart& d : moving_darts)
			{
				topo_render_->reset_dart_color(d);
				d = phi_1(*mesh_, d);
				topo_render_->set_dart_color(d, get_color(d));
			}
			force_update();
		}

		if (ImGui::Button("phi 2"))
		{
			for (Dart& d : moving_darts)
			{
				topo_render_->reset_dart_color(d);
				d = phi2(*mesh_, d);
				topo_render_->set_dart_color(d, get_color(d));
			}
			force_update();
		}

		if (ImGui::Button("phi 3"))
		{
			cgogn::Dart new_moving_dart_ = cgogn::phi3(*mesh_, moving_dart_);
			topo_render_->set_dart_color(new_moving_dart_, moving_color_);
			topo_render_->reset_dart_color(moving_dart_);
			moving_dart_ = new_moving_dart_;
			force_update();
		}
	}
	MRMesh* mesh_;
	cgogn::ui::View* view_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	cgogn::ui::MeshProvider<MRMesh>* mesh_provider_;
	cgogn::ui::VolumeRender<MRMesh>* vol_render_;
	cgogn::ui::TopoRender<MRMesh>* topo_render_;
	cgogn::Dart moving_dart_;
	std::vector<Dart> moving_darts;
	Eigen::Vector4f moving_color_;
	float expl_vol_;
};

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
	app.set_window_title("EMR Volume");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::MeshProvider<MRMesh> mrmp(app);
	cgogn::ui::VolumeRender<MRMesh> vr(app);
	cgogn::ui::VolumeSelection<MRMesh> vs(app);
	cgogn::ui::TopoRender<MRMesh> tr(app);
	LocalInterface interf(app);

	cgogn::ui::VolumeEMRModeling<MRMesh> vmrm(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&mrmp);
	v1->link_module(&vr);
	v1->link_module(&vs);
	v1->link_module(&vmrm);
	v1->link_module(&tr);
	v1->link_module(&interf);

	/*cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mp);
	v2->link_module(&mrmp);
	v2->link_module(&vr);
	v2->link_module(&vs);*/

	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	MRMesh* mrm = vmrm.create_mrmesh(*m, mp.mesh_name(*m));
	// MRMesh* mrm2 = vmrm.create_mrmesh(*m, mp.mesh_name(m));
	cgogn::index_cells<Mesh::Face>(*mrm);
	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);

	m->add_resolution();
	mrm->change_resolution_level(1);
	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");

	vmrm.subdivide(*mrm, position.get());
	m->add_resolution();
	mrm->change_resolution_level(2);
	vmrm.subdivide(*mrm, position.get());

	auto& md = mrmp.mesh_data(*mrm);
	md.template add_cells_set<Edge>();

	mrmp.set_mesh_bb_vertex_position(*mrm, position);

	vr.set_vertex_position(*v1, *mrm, position);

	// std::srand(std::time(nullptr));
	std::srand(2124512438);
	mrm->current_level_ = 0;
	std::vector<Volume> vol_vec;
	std::vector<Volume> vol_vec_simpl;
	cgogn::CellMarker<MRMesh, Volume> vm(*mrm);

	vmrm.changed_connectivity(*mrm, position.get());
	interf.mesh_ = mrm;
	interf.vertex_position_ = position;

	return app.launch();
}
