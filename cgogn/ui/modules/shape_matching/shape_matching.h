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
class ShapeMatching : public Module
{

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	enum SelectingCell
	{
		VertexSelect = 0,
		EdgeSelect
	};

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), init_vertex_position_(nullptr), vertex_forces_(nullptr),
			  vertex_masse_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  selected_vertices_set_(nullptr), selected_edges_set_(nullptr), selecting_cell_(VertexSelect)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 0, 0, 0.65);
			param_edge_->width_ = 2.0f;
			param_edge_->set_vbos({&selected_edges_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell(
					[&](Vertex v) { selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v)); });
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		void update_selected_edges_vbo()
		{
			if (selected_edges_set_)
			{
				std::vector<Vec3> selected_edges_position;
				selected_edges_position.reserve(selected_edges_set_->size() * 2);
				selected_edges_set_->foreach_cell([&](Edge e) {
					std::vector<Vertex> vertices = incident_vertices(*mesh_, e);
					selected_edges_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[0]));
					selected_edges_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[1]));
				});
				rendering::update_vbo(selected_edges_position, &selected_edges_vbo_);
			}
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> init_vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_forces_;
		std::shared_ptr<Attribute<double>> vertex_masse_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO selected_vertices_vbo_;
		rendering::VBO selected_edges_vbo_;

		CellsSet<MESH, Vertex>* selected_vertices_set_;
		CellsSet<MESH, Edge>* selected_edges_set_;

		SelectingCell selecting_cell_;
	};

public:
	ShapeMatching(const App& app)
		: Module(app, "ShapeMatching (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  sm_solver_(0.5f), running_(false)
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
						p.update_selected_vertices_vbo();
						p.update_selected_edges_vbo();
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
			p.update_selected_vertices_vbo();
			p.update_selected_edges_vbo();
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

	void key_press_event(View* v, int32 key_code)
	{
		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];
			this->f_keypress(v, selected_mesh_, key_code, p.selected_vertices_set_, p.selected_edges_set_);
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

	void interface() override
	{

		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

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
					sm_solver_.update_topo(*selected_mesh_);
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
					sm_solver_.update_topo(*selected_mesh_);
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
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHAPE_MATCHING_H_
