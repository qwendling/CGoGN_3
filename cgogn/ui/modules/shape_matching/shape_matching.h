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

#ifndef CGOGN_MODULE_SHAPE_MATCHING_H_
#define CGOGN_MODULE_SHAPE_MATCHING_H_

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

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>
#include <cgogn/simulation/algos/Simulation_solver.h>
#include <cgogn/simulation/algos/shape_matching/shape_matching.h>

#include <boost/synapse/connect.hpp>
#include <imgui/imgui.h>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class ShapeMatching : public ViewModule
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
			: vertex_position_(nullptr), init_vertex_position_(nullptr), vertex_forces_(nullptr),
			  vertex_masse_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  have_selected_vertex_(false), move_vertex_(0, 0, 0)
		{
			param_move_vertex_ = rendering::ShaderPointSprite::generate_param();
			param_move_vertex_->color_ = rendering::GLColor(1, 1, 0, 0.65);
			param_move_vertex_->set_vbos({&move_vertex_vbo_});

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 0, 1, 0.65);
			param_edge_->width_ = 2.0f;
			param_edge_->set_vbos({&edges_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		void update_move_vertex_vbo()
		{
			std::vector<Vec3> vertices_position;
			if (have_selected_vertex_)
			{
				vertices_position.push_back(move_vertex_);
				vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_.get(), selected_vertex_));
			}
			rendering::update_vbo(vertices_position, &move_vertex_vbo_);
			rendering::update_vbo(vertices_position, &edges_vbo_);
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> init_vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_forces_;
		std::shared_ptr<Attribute<double>> vertex_masse_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_move_vertex_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO move_vertex_vbo_;
		rendering::VBO edges_vbo_;

		Vec3 move_vertex_;
		bool have_selected_vertex_;
		Vertex selected_vertex_;
	};

public:
	ShapeMatching(const App& app)
		: ViewModule(app, "ShapeMatching (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_view_(app.current_view()), sm_solver_(0.5f), running_(false)
	{
		f_keypress = [](View*, MESH*, int32, CellsSet<MESH, Vertex>*, CellsSet<MESH, Edge>*) {};
	}

	~ShapeMatching()
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

	void set_init_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& init_vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.init_vertex_position_ = init_vertex_position;
	}

	void set_vertex_force(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_forces)
	{
		Parameters& p = parameters_[&m];

		p.vertex_forces_ = vertex_forces;
		simu_solver.forces_ext_ = vertex_forces;
	}
	void set_vertex_masse(const MESH& m, const std::shared_ptr<Attribute<double>>& vertex_masse)
	{
		Parameters& p = parameters_[&m];

		p.vertex_masse_ = vertex_masse;
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &ShapeMatching<MESH>::init_mesh));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		Parameters& p = parameters_[selected_mesh_];
		if (button == 1 && p.have_selected_vertex_)
		{
			p.move_vertex_ =
				view->pixel_scene_(x, y, value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_));
			p.update_move_vertex_vbo();
			view->request_update();
		}
		if (selected_mesh_ && view->shift_pressed())
		{
			if (p.vertex_position_)
			{

				rendering::GLVec3d near = view->unproject(x, y, 0.0);
				rendering::GLVec3d far = view->unproject(x, y, 1.0);
				Vec3 A{near.x(), near.y(), near.z()};
				Vec3 B{far.x(), far.y(), far.z()};
				std::vector<Vertex> picked;
				cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
				if (!picked.empty())
				{
					p.selected_vertex_ = picked[0];
					p.have_selected_vertex_ = true;
					p.move_vertex_ = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), picked[0]);
					p.update_move_vertex_vbo();
					view->request_update();
				}
			}
		}
	}

	void key_press_event(View* v, int32 key_code)
	{
		if (key_code == GLFW_KEY_LEFT_CONTROL)
		{
			v->lock_rotation_ = true;
			can_move_vertex_ = true;
		}
	}

	void key_release_event(View* v, int32 key_code)
	{
		if (key_code == GLFW_KEY_LEFT_CONTROL)
		{
			v->lock_rotation_ = false;
			can_move_vertex_ = false;
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y)
	{
		Parameters& p = parameters_[selected_mesh_];
		if (p.have_selected_vertex_ && can_move_vertex_)
		{
			p.move_vertex_ =
				view->pixel_scene_(x, y, value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_));
			p.update_move_vertex_vbo();
			view->request_update();
		}
	}

	void start()
	{
		running_ = true;

		launch_thread([this]() {
			while (this->running_)
			{
				Parameters& p = parameters_[selected_mesh_];
				simu_solver.compute_time_step(*selected_mesh_, p.vertex_position_.get(), p.vertex_masse_.get(), 0.005);
				need_update_ = true;
				std::this_thread::sleep_for(std::chrono::milliseconds(5));
			}
		});

		// app_.start_timer(5000, [this]() -> bool { return !running_; });
	}

	void stop()
	{
		running_ = false;
	}

	void step()
	{
		Parameters& p = parameters_[selected_mesh_];
		simu_solver.compute_time_step(*selected_mesh_, p.vertex_position_.get(), p.vertex_masse_.get(), 0.005);
		need_update_ = true;
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.have_selected_vertex_ && p.param_move_vertex_->vao_initialized())
			{
				p.param_move_vertex_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_move_vertex_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, 2);
				p.param_move_vertex_->release();
			}

			if (p.have_selected_vertex_ && p.param_edge_->vao_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_LINES, 0, 2);
				p.param_edge_->release();
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
					simu_solver.init_solver(*selected_mesh_, &sm_solver_);
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			Parameters& p = parameters_[selected_mesh_];

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
			if (ImGui::BeginCombo("Forces", p.vertex_forces_ ? p.vertex_forces_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_,
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == p.vertex_forces_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														set_vertex_force(*selected_mesh_, attribute);
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (p.vertex_forces_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##force"))
					set_vertex_force(*selected_mesh_, nullptr);
			}
			if (ImGui::BeginCombo("Init vertex position",
								  p.init_vertex_position_ ? p.init_vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(
					*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
						bool is_selected = attribute == p.init_vertex_position_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						{
							set_init_vertex_position(*selected_mesh_, attribute);
							if (p.vertex_masse_)
								sm_solver_.init_solver(*selected_mesh_, p.init_vertex_position_, p.vertex_masse_);
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
				ImGui::EndCombo();
			}
			if (p.init_vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##init"))
					set_init_vertex_position(*selected_mesh_, nullptr);
			}
			if (ImGui::BeginCombo("vertex masse", p.vertex_masse_ ? p.vertex_masse_->name().c_str() : "-- select --"))
			{
				foreach_attribute<double, Vertex>(
					*selected_mesh_, [&](const std::shared_ptr<Attribute<double>>& attribute) {
						bool is_selected = attribute == p.vertex_masse_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						{
							set_vertex_masse(*selected_mesh_, attribute);
							if (p.init_vertex_position_)
								sm_solver_.init_solver(*selected_mesh_, p.init_vertex_position_, p.vertex_masse_);
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
				ImGui::EndCombo();
			}
			if (p.vertex_masse_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##masse"))
					set_vertex_masse(*selected_mesh_, nullptr);
			}
			if (p.init_vertex_position_ && p.vertex_position_ && p.vertex_masse_)
			{
				if (ImGui::Button("init initial pos"))
				{
					foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						value<Vec3>(*selected_mesh_, p.init_vertex_position_.get(), v) =
							value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
						return true;
					});
					sm_solver_.update_topo(*selected_mesh_, {});
				}
			}
			if (p.vertex_masse_)
			{
				if (ImGui::Button("init masse"))
				{
					foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						value<double>(*selected_mesh_, p.vertex_masse_.get(), v) = 1.0f;
						return true;
					});
					sm_solver_.update_topo(*selected_mesh_, {});
				}
			}
			if (ImGui::Button("new attribute Vec3"))
			{
				static uint32 nb_new_attribute_Vec3 = 0;
				add_attribute<Vec3, Vertex>(*selected_mesh_,
											"SM_attributte_vec3_" + std::to_string(nb_new_attribute_Vec3++));
			}
			if (ImGui::Button("new attribute double"))
			{
				static uint32 nb_new_attribute_double = 0;
				add_attribute<double, Vertex>(*selected_mesh_,
											  "SM_attributte_double_" + std::to_string(nb_new_attribute_double++));
			}
			if (p.vertex_position_ && p.init_vertex_position_ && p.vertex_forces_ && p.vertex_masse_)
			{
				ImGui::Separator();

				MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);
				Parameters& p = parameters_[selected_mesh_];

				if (!running_)
				{
					if (ImGui::Button("play"))
					{
						start();
					}
					if (ImGui::Button("step"))
					{
						step();
					}
				}
				else
				{
					if (ImGui::Button("stop"))
					{
						stop();
					}
				}
				if (need_update_)
				{
					mesh_provider_->emit_attribute_changed(selected_mesh_, p.vertex_position_.get());
					need_update_ = false;
				}
				double min = 0, max = 1;
				ImGui::SliderScalar("Level", ImGuiDataType_Double, &sm_solver_.stiffness_, &min, &max);
			}
		}

		ImGui::End();
	}

public:
	MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
	simulation::shape_matching_constraint_solver<MESH> sm_solver_;
	simulation::Simulation_solver<MESH> simu_solver;
	bool running_;
	bool need_update_;
	bool can_move_vertex_;
	View* selected_view_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHAPE_MATCHING_H_
