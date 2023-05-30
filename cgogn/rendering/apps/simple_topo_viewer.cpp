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


#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/rendering/ui_modules/topo_render.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/functions/distance.h>
#include <cgogn/modeling/algos/volume_utils.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/types/cmap/dart.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Mesh = cgogn::CMap3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;
using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;
using namespace cgogn::numerics;

class LocalInterface : public cgogn::ui::ViewModule
{

public:
	LocalInterface(const cgogn::ui::App& app)
		: cgogn::ui::ViewModule(app, "LocalInterface"), mesh_(nullptr), vertex_position_(nullptr), mesh_provider_(nullptr),
		  vol_render_(nullptr), topo_render_(nullptr), moving_color_(1.0f, 0.0f, 1.0f,1.0f)
	{
		view_ = app.current_view();
	}

	~LocalInterface()
	{
	}

	void create()
	{
		mesh_ = mesh_provider_->add_mesh("pyra_and_hexa");
		vertex_position_ = cgogn::add_attribute<Vec3, Vertex>(*mesh_, "position");

		d_pyra_ = cgogn::add_pyramid(*mesh_, 4, false).dart;
		d_hexa_ = cgogn::add_prism(*mesh_, 4, false).dart;

		auto setPosV = [&](cgogn::Dart d, const Vec3& P) { cgogn::value<Vec3>(*mesh_, vertex_position_, Vertex(d)) = P; };
		cgogn::Dart dp = d_pyra_;
		cgogn::Dart dh = d_hexa_;
		for (int i = 0; i < 4; ++i)
		{
			cgogn::phi3_sew(*mesh_, dp, dh);
			dp = cgogn::phi1(*mesh_, dp);
			dh = cgogn::phi_1(*mesh_, dh);
		}

		cgogn::close(*mesh_, false);

		cgogn::index_cells<Vertex>(*mesh_);
		setPosV(dp, Vec3(-1, -1, -1));
		dp = cgogn::phi1(*mesh_, dp);
		setPosV(dp, Vec3(-1.4, 1.4, -1));
		dp = cgogn::phi1(*mesh_, dp);
		setPosV(dp, Vec3(1, 1, -1));
		dp = cgogn::phi1(*mesh_, dp);
		setPosV(dp, Vec3(1.4, -1.4, -1));
		setPosV(cgogn::phi<2, -1>(*mesh_, dp), Vec3(0, 0, 1));

		dh = cgogn::phi<2, 1, 1, 2>(*mesh_, dh);
		setPosV(dh, Vec3(-1.2, -1.2, -3));
		dh = cgogn::phi1(*mesh_, dh);
		setPosV(dh, Vec3(-1, 1, -3));
		dh = cgogn::phi1(*mesh_, dh);
		setPosV(dh, Vec3(1.2, 1.2, -3));
		dh = cgogn::phi1(*mesh_, dh);
		setPosV(dh, Vec3(1, -1, -3));

		cgogn::index_cells<Volume>(*mesh_);

		vol_render_->set_vertex_position(*app_.current_view(), *mesh_, vertex_position_);
		
		topo_render_->set_selected_mesh(*mesh_);
		topo_render_->set_vertex_position(*app_.current_view(), *mesh_, vertex_position_);
		topo_render_->set_dart_color(d_pyra_, {1.0f, 0.0f, 1.0f, 1.0f});
		mesh_provider_->set_mesh_bb_vertex_position(*mesh_, vertex_position_);


		mesh_provider_->emit_connectivity_changed(*mesh_);
	}

protected:

	void force_update()
	{
		for (cgogn::ui::View* v : linked_views_)
			v->request_update();
	}


	void init() override
	{
		mesh_provider_ = static_cast<cgogn::ui::MeshProvider<Mesh>*>(
			app_.module("MeshProvider (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));

		vol_render_ = static_cast<cgogn::ui::VolumeRender<Mesh>*>(
			app_.module("VolumeRender (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));
		topo_render_ = static_cast<cgogn::ui::TopoRender<Mesh>*>(
			app_.module("TopoRender (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));

	}

	void left_panel() override
	{
		if (ImGui::SliderFloat("Explode", &expl_vol_, 0.01f, 1.0f))
		{
			vol_render_->set_volume_explode(*view_,*mesh_,expl_vol_);
			topo_render_->set_volume_explode(expl_vol_+0.02f);
			force_update();
		}

		if (ImGui::Button("init moving"))
		{
			if (!moving_dart_.is_nil())
				topo_render_->set_dart_color(moving_dart_, moving_color_);
			moving_dart_ = d_hexa_;
			topo_render_->set_dart_color(moving_dart_, moving_color_);
			force_update();
		}

		if (ImGui::Button("ph1"))
		{
			cgogn::Dart new_moving_dart_ = cgogn::phi1(*mesh_, moving_dart_);
			topo_render_->set_dart_color(new_moving_dart_, moving_color_);
			topo_render_->reset_dart_color(moving_dart_);
			moving_dart_ = new_moving_dart_;
			force_update();
		}

		if (ImGui::Button("phi2"))
		{
			cgogn::Dart new_moving_dart_ = cgogn::phi2(*mesh_, moving_dart_);
			topo_render_->set_dart_color(new_moving_dart_, moving_color_);
			topo_render_->reset_dart_color(moving_dart_);
			moving_dart_ = new_moving_dart_;
			force_update();
		}

		if (ImGui::Button("phi3"))
		{
			cgogn::Dart new_moving_dart_ = cgogn::phi3(*mesh_, moving_dart_);
			topo_render_->set_dart_color(new_moving_dart_, moving_color_);
			topo_render_->reset_dart_color(moving_dart_);
			moving_dart_ = new_moving_dart_;
			force_update();
		}
	}

private:
	Mesh* mesh_;
	cgogn::ui::View* view_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	cgogn::ui::MeshProvider<Mesh>* mesh_provider_;
	cgogn::ui::VolumeRender<Mesh>* vol_render_;
	cgogn::ui::TopoRender<Mesh>* topo_render_;
	cgogn::Dart d_hexa_;
	cgogn::Dart d_pyra_;
	cgogn::Dart moving_dart_;
	Eigen::Vector4f moving_color_;
	float expl_vol_;
};


int main(int argc, char** argv)
{
	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple volume viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::VolumeRender<Mesh> vr(app);
	cgogn::ui::TopoRender<Mesh> tr(app);
	LocalInterface interf(app);
	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&vr);
	v1->link_module(&tr);
	v1->link_module(&interf);

	interf.create();
	
	return app.launch();
}
