/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_TOPO_RENDER_H_
#define CGOGN_MODULE_TOPO_RENDER_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/rendering/topo_drawer.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class TopoRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "Topo_Render can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Volume = typename mesh_traits<MESH>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters() : vertex_position_(nullptr), render_topo_(true)
		{
			topo_drawer_ = std::make_unique<rendering::TopoDrawer>();
			topo_renderer_ = topo_drawer_->generate_renderer();
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		inline void update_topo(const MESH& m)
		{
			if constexpr (mesh_traits<MESH>::dimension == 2)
				topo_drawer_->update2D(m, vertex_position_.get());
			else if constexpr (mesh_traits<MESH>::dimension == 3)
				topo_drawer_->update3D(m, vertex_position_.get());
		}

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::unique_ptr<rendering::TopoDrawer> topo_drawer_;
		std::unique_ptr<rendering::TopoDrawer::Renderer> topo_renderer_;
		bool render_topo_;
	};

public:
	TopoRender(const App& app)
		: ViewModule(app, "TopoRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~TopoRender()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		for (View* v : linked_views_)
		{
			auto& p = parameters_[v][m];

			p.topo_drawer_->reset_all_colors(rendering::GLVec3(1.0f,1.0f,1.0f));

			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);

			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, v, m]() {
					Parameters& p = parameters_[v][m];
					p.update_topo(*m);
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
							p.update_topo(*m);
						v->request_update();
					}));
		}
	}

public:
	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
			p.update_topo(m);
		v.request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &TopoRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_topo_)
			{
				p.topo_renderer_->draw(proj_matrix, view_matrix);
			}
		}
	}

	void left_panel() override
	{
		bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Mesh", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
								[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
									set_vertex_position(*selected_view_, *selected_mesh_, attribute);
									});
			
			ImGui::Separator();
			need_update |= ImGui::Checkbox("Topo", &p.render_topo_);

			if (p.render_topo_)
			{
				ImGui::Separator();
				need_update |= (ImGui::Checkbox("Render darts", &p.topo_renderer_->render_darts_));
				if (p.topo_renderer_->render_darts_)
				{
					if (ImGui::ColorEdit3("colorDarts", p.topo_drawer_->dart_color_.data(),
											ImGuiColorEditFlags_NoInputs))
					{
						const auto& col = p.topo_drawer_->dart_color_;
						p.topo_drawer_->reset_all_colors(rendering::GLVec3(col[0], col[1], col[2]));
						need_update = true;
					}
				}
				//ImGui::TextUnformatted("Volume parameters");
				if (mesh_traits<MESH>::dimension >= 2)
				{
					need_update |= ImGui::Checkbox("Render phi2", &p.topo_renderer_->render_phi2_);
					if (p.topo_renderer_->render_phi2_)
					{
						need_update |= ImGui::ColorEdit3("colorPhi2", p.topo_drawer_->phi2_color_.data(),
														 ImGuiColorEditFlags_NoInputs);
					}
				}

				if (mesh_traits<MESH>::dimension >= 3)
				{
					need_update |= ImGui::Checkbox("Render phi3", &p.topo_renderer_->render_phi3_);
					if (p.topo_renderer_->render_phi3_)
					{
						need_update |= ImGui::ColorEdit3("colorPhi3", p.topo_drawer_->phi3_color_.data(),
														 ImGuiColorEditFlags_NoInputs);
					}
				}

				need_update |= ImGui::SliderFloat("Width", &p.topo_renderer_->width_, 1.0f, 5.0f);

				need_update |= ImGui::SliderFloat("explodeEdges", &(p.topo_drawer_->shrink_e_), 0.01f, 1.0f);
				if (mesh_traits<MESH>::dimension >= 2)
					need_update |= ImGui::SliderFloat("explodeFaces", &(p.topo_drawer_->shrink_f_), 0.01f, 1.0f);
				if (mesh_traits<MESH>::dimension >= 3)
					need_update |= ImGui::SliderFloat("explodeVolumes", &(p.topo_drawer_->shrink_v_), 0.01f, 1.0f);
			}
			if (need_update)
				p.update_topo(*selected_mesh_);
		}

		float64 remain = mesh_provider_->mesh_data(*selected_mesh_).outlined_until_ - App::frame_time_;
		if (remain > 0)
			need_update = true;

		if (need_update)
		{
			for (View* v : linked_views_)
				v->request_update();
		}
	}

public:

	void set_dart_color(Dart d, const Eigen::Vector4f& color)
	{
		Parameters& p = parameters_[selected_view_][selected_mesh_];
		p.topo_drawer_->update_color(d, color);

	}

	void reset_dart_color(Dart d)
	{
		Parameters& p = parameters_[selected_view_][selected_mesh_];
		p.topo_drawer_->update_color(d, p.topo_drawer_->dart_color_);
	}

	void set_selected_mesh(const MESH& m)
	{
		selected_mesh_ = &m;
	}

	void set_volume_explode(float expl)
	{
		Parameters& p = parameters_[selected_view_][selected_mesh_];
		p.topo_drawer_->shrink_v_ = expl;
		p.update_topo(*selected_mesh_);
	}

private:
	View* selected_view_;
	const MESH* selected_mesh_;
	std::unordered_map<View*, std::unordered_map<const MESH*, Parameters>> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif
