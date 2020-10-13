/*******************************************************************************
 * CGoGN                                                                        *
 * Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_MODULE_LINKED_VOLUMES_H_
#define CGOGN_MODULE_LINKED_VOLUMES_H_

#include <GLFW/glfw3.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/frame_manipulator.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>
#include <cgogn/simulation/algos/Simulation_solver.h>
#include <cgogn/simulation/algos/Simulation_solver_multiresolution.h>
#include <cgogn/simulation/algos/linked_volumes/linked_volumes.h>

#include <boost/synapse/connect.hpp>
#include <imgui/imgui.h>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class LinkedVolumes : public ViewModule
{

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  show_frame_manipulator_(false), manipulating_frame_(false)
		{

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 0, 1, 0.65);
			param_edge_->width_ = 2.0f;
			param_edge_->set_vbos({&edges_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO move_vertex_vbo_;
		rendering::VBO edges_vbo_;
		Vertex selected_vertex_;

		rendering::FrameManipulator frame_manipulator_;
		bool show_frame_manipulator_;
		bool manipulating_frame_;
	};

public:
	LinkedVolumes(const App& app)
		: ViewModule(app, "LinkedVolumes (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_view_(app.current_view())
	{
		f_keypress = [](View*, MESH*, int32, CellsSet<MESH, Vertex>*, CellsSet<MESH, Edge>*) {};
	}

	~LinkedVolumes()
	{
	}

	std::function<void(View*, MESH*, int32, CellsSet<MESH, Vertex>*, CellsSet<MESH, Edge>*)> f_keypress;

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6.0;
					}
				}));
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = geometry::mean_edge_length(m, p.vertex_position_.get()) / 6.0;
		}
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &LinkedVolumes<MESH>::init_mesh));
	}

	void mouse_press_event(View* view, int32, int32 x, int32 y) override
	{
		Parameters& p = parameters_[selected_mesh_];
		if (p.manipulating_frame_)
		{
			auto [P, Q] = view->pixel_ray(x, y);
			p.frame_manipulator_.pick(x, y, P, Q);
			view->request_update();
		}
	}

	void key_press_event(View*, int32 key_code)
	{
		if (key_code == GLFW_KEY_C)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				if (p.show_frame_manipulator_)
					p.manipulating_frame_ = true;
			}
		}
		if (key_code == GLFW_KEY_T)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				Vec3 pos;
				p.frame_manipulator_.get_position(pos);
				Vec3 a;
				p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, a);
				double d = pos.dot(a);
				lv_.compute_cut_plan(a, d, p.vertex_position_.get(), [&](std::pair<Vertex, Vertex> p) -> bool {
					for (auto attr : list_update_attribute)
					{
						value<Vec3>(*selected_mesh_, attr, p.second) = value<Vec3>(*selected_mesh_, attr, p.first);
					}
					return true;
				});
				for (auto attr : list_update_attribute)
				{
					mesh_provider_->emit_attribute_changed(selected_mesh_, attr.get());
					std::cout << "hello" << std::endl;
				}
				mesh_provider_->emit_connectivity_changed(selected_mesh_);
			}
		}
	}

	void key_release_event(View*, int32 key_code)
	{
		if (key_code == GLFW_KEY_C)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				p.manipulating_frame_ = false;
			}
		}
	}

	void mouse_release_event(View* view, int32, int32, int32) override
	{
		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];
			p.frame_manipulator_.release();
			view->request_update();
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y)
	{
		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];
			bool leftpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_LEFT);
			bool rightpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_RIGHT);
			if (p.manipulating_frame_ && (rightpress || leftpress))
			{
				p.frame_manipulator_.drag(leftpress, x, y);
				view->stop_event();
				view->request_update();
			}
		}
	}

	void draw(View* view) override
	{
		if (selected_mesh_)
		{
			auto& m = selected_mesh_;
			auto& p = parameters_[selected_mesh_];
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.show_frame_manipulator_)
			{
				Scalar size = (md->bb_max_ - md->bb_min_).norm() / 10;
				p.frame_manipulator_.set_size(size);
				p.frame_manipulator_.draw(true, true, proj_matrix, view_matrix);
			}
		}
	}

	void interface() override
	{

		std::stringstream ss;
		ss << std::setw(6) << std::fixed << std::setprecision(2) << App::fps();
		std::string str_fps = ss.str() + " fps";
		ImGui::Text(str_fps.c_str());

		if (ImGui::BeginCombo("View", selected_view_->name().c_str()))
		{
			for (View* v : linked_views_)
			{
				bool is_selected = v == selected_view_;
				if (ImGui::Selectable(v->name().c_str(), is_selected))
					selected_view_ = v;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}

		if (ImGui::ListBoxHeader("Mesh"))
		{
			mesh_provider_->foreach_mesh([this](MESH* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
				{
					selected_mesh_ = m;
					lv_.init_mesh(m);
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			Parameters& p = parameters_[selected_mesh_];

			need_update_ |= ImGui::Checkbox("Show plane", &p.show_frame_manipulator_);

			if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_,
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == p.vertex_position_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														set_vertex_position(*selected_mesh_, attribute);
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (p.vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					set_vertex_position(*selected_mesh_, nullptr);
			}
			foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
				bool is_selected = false;
				auto it = list_update_attribute.begin();
				for (; it != list_update_attribute.end(); ++it)
				{
					if (*it == attribute)
					{
						is_selected = true;
						break;
					}
				}
				if (ImGui::Checkbox(attribute->name().c_str(), &is_selected))
				{
					if (is_selected)
					{
						list_update_attribute.push_back(attribute);
					}
					else
					{
						list_update_attribute.erase(it);
					}
				}
			});
		}
	}

public:
	MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
	simulation::Linked_volumes<MESH> lv_;
	bool need_update_;
	View* selected_view_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> list_update_attribute;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_LINKED_VOLUMES_H_
