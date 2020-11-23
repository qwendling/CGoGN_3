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

#ifndef CGOGN_MODULE_ANIMATION_MULTIRESOLUTION_H_
#define CGOGN_MODULE_ANIMATION_MULTIRESOLUTION_H_

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
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_forces.h>
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_plastique.h>
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_spring.h>
#include <cgogn/simulation/algos/shape_matching/shape_matching.h>

#include <boost/synapse/connect.hpp>
#include <imgui/imgui.h>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MR_MESH>
class AnimationMultiresolution : public ViewModule
{

	template <typename T>
	using Attribute = typename mesh_traits<MR_MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MR_MESH>::Vertex;
	using Edge = typename mesh_traits<MR_MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_relative_position_(nullptr), init_vertex_position_(nullptr),
			  vertex_forces_(nullptr), vertex_masse_(nullptr), vertex_parents_(nullptr), vertex_scale_factor_(1.0),
			  sphere_scale_factor_(10.0), have_selected_vertex_(false), move_vertex_(0, 0, 0),
			  show_frame_manipulator_(false), manipulating_frame_(false)
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
			if (have_selected_vertex_)
			{
				std::vector<Vec3> vertices_position;
				vertices_position.push_back(move_vertex_);
				vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_.get(), selected_vertex_));

				rendering::update_vbo(vertices_position, &move_vertex_vbo_);
				rendering::update_vbo(vertices_position, &edges_vbo_);
			}
		}

		MR_MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_relative_position_;
		std::shared_ptr<Attribute<Vec3>> init_vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_forces_;
		std::shared_ptr<Attribute<double>> vertex_masse_;
		std::shared_ptr<Attribute<std::array<Vertex, 3>>> vertex_parents_;

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
		rendering::FrameManipulator frame_manipulator_;
		bool show_frame_manipulator_;
		bool manipulating_frame_;
	};

public:
	AnimationMultiresolution(const App& app)
		: ViewModule(app, "Animation_multiresolution (" + std::string{mesh_traits<MR_MESH>::name} + ")"),
		  mecanical_mesh_(nullptr), selected_view_(app.current_view()), sm_solver_(0.9f), running_(false), ps_(0.9f),
		  geometric_mesh_(nullptr), modif_topo_(false)
	{
		f_keypress = [](View*, MR_MESH*, int32, CellsSet<MR_MESH, Vertex>*, CellsSet<MR_MESH, Edge>*) {};
	}

	~AnimationMultiresolution()
	{
	}

	std::function<void(View*, MR_MESH*, int32, CellsSet<MR_MESH, Vertex>*, CellsSet<MR_MESH, Edge>*)> f_keypress;

private:
	void init_mesh(MR_MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MR_MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6.0;
					}
				}));
	}

public:
	void set_vertex_position(const MR_MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = geometry::mean_edge_length(m, p.vertex_position_.get()) / 6.0;
		}
	}

	void set_init_vertex_position(const MR_MESH& m, const std::shared_ptr<Attribute<Vec3>>& init_vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.init_vertex_position_ = init_vertex_position;
	}

	void set_vertex_force(const MR_MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_forces)
	{
		Parameters& p = parameters_[&m];

		p.vertex_forces_ = vertex_forces;
		simu_solver.forces_ext_ = vertex_forces;
	}
	void set_vertex_masse(const MR_MESH& m, const std::shared_ptr<Attribute<double>>& vertex_masse)
	{
		Parameters& p = parameters_[&m];

		p.vertex_masse_ = vertex_masse;
	}

	void set_vertex_parents(const MR_MESH& m, const std::shared_ptr<Attribute<std::array<Vertex, 3>>>& vertex_parents)
	{
		Parameters& p = parameters_[&m];

		p.vertex_parents_ = vertex_parents;
		simu_solver.parents_ = vertex_parents;
	}

	void set_vertex_relative_position(const MR_MESH& m,
									  const std::shared_ptr<Attribute<Vec3>>& vertex_relative_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_relative_position_ = vertex_relative_position;
		simu_solver.relative_pos_ = vertex_relative_position;
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MR_MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MR_MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MR_MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MR_MESH>::mesh_added>(
			mesh_provider_, this, &AnimationMultiresolution<MR_MESH>::init_mesh));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		Parameters& p = parameters_[mecanical_mesh_];
		if (button == 1 && p.have_selected_vertex_)
		{
			p.move_vertex_ =
				view->pixel_scene_(x, y, value<Vec3>(*mecanical_mesh_, p.vertex_position_.get(), p.selected_vertex_));
			p.update_move_vertex_vbo();
			view->request_update();
		}
		if (p.manipulating_frame_)
		{
			auto [P, Q] = view->pixel_ray(x, y);
			p.frame_manipulator_.pick(x, y, P, Q);
			view->request_update();
		}
		view->pixel_ray(x, y);
		if (mecanical_mesh_ && view->shift_pressed())
		{
			if (p.vertex_position_)
			{
				rendering::GLVec3d near = view->unproject(x, y, 0.0);
				rendering::GLVec3d far = view->unproject(x, y, 1.0);
				Vec3 A{near.x(), near.y(), near.z()};
				Vec3 B{far.x(), far.y(), far.z()};
				std::vector<Vertex> picked;
				cgogn::geometry::picking(*mecanical_mesh_, p.vertex_position_.get(), A, B, picked);
				if (!picked.empty())
				{
					p.selected_vertex_ = picked[0];
					p.have_selected_vertex_ = true;
					p.move_vertex_ = value<Vec3>(*mecanical_mesh_, p.vertex_position_.get(), picked[0]);
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
		if (key_code == GLFW_KEY_C)
		{
			if (mecanical_mesh_)
			{
				Parameters& p = parameters_[mecanical_mesh_];
				if (p.show_frame_manipulator_)
					p.manipulating_frame_ = true;
			}
		}
	}

	void key_release_event(View* v, int32 key_code)
	{
		if (key_code == GLFW_KEY_LEFT_CONTROL)
		{
			v->lock_rotation_ = false;
			can_move_vertex_ = false;
		}
		if (key_code == GLFW_KEY_C)
		{
			if (mecanical_mesh_)
			{
				Parameters& p = parameters_[mecanical_mesh_];
				p.manipulating_frame_ = false;
			}
		}
	}

	void mouse_release_event(View* view, int32, int32, int32) override
	{
		if (mecanical_mesh_)
		{
			Parameters& p = parameters_[mecanical_mesh_];
			p.frame_manipulator_.release();
			view->request_update();
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y)
	{
		Parameters& p = parameters_[mecanical_mesh_];
		if (p.have_selected_vertex_ && can_move_vertex_)
		{
			p.move_vertex_ =
				view->pixel_scene_(x, y, value<Vec3>(*mecanical_mesh_, p.vertex_position_.get(), p.selected_vertex_));
			p.update_move_vertex_vbo();
			view->request_update();
		}
		bool leftpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_LEFT);
		bool rightpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_RIGHT);
		if (p.manipulating_frame_ && (rightpress || leftpress))
		{
			p.frame_manipulator_.drag(leftpress, x, y);
			view->stop_event();
			view->request_update();
		}
	}
#define TIME_STEP 0.005f
	void start()
	{
		running_ = true;

		Parameters& p = parameters_[mecanical_mesh_];

		launch_thread([this, &p]() {
			typename MR_MESH::CMAP& map = static_cast<typename MR_MESH::CMAP&>(*mecanical_mesh_);
			while (this->running_)
			{

				map.start_writer();
				if (p.have_selected_vertex_)
				{
					Vec3 pos = value<Vec3>(*mecanical_mesh_, p.vertex_position_.get(), p.selected_vertex_);
					double m = value<double>(*mecanical_mesh_, p.vertex_masse_.get(), p.selected_vertex_);
					value<Vec3>(*mecanical_mesh_, p.vertex_forces_.get(), p.selected_vertex_) =
						m * (p.move_vertex_ - pos) / TIME_STEP;
				}

				/*parallel_foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool {
					double m = value<double>(*mecanical_mesh_, p.vertex_masse_.get(), v);
					value<Vec3>(*mecanical_mesh_, p.vertex_forces_.get(), v) += m * Vec3(0, 0, -9.81);
					return true;
				});*/

				simu_solver.compute_time_step(*geometric_mesh_, p.vertex_position_.get(), p.vertex_masse_.get(),
											  TIME_STEP, modif_topo_);
				if (true || p.show_frame_manipulator_)
				{
					Vec3 position;
					Vec3 axis_z;
					p.frame_manipulator_.get_position(position);
					p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, axis_z);
					double d = position.dot(axis_z);
					parallel_foreach_cell(mecanical_mesh_->m_, [&](Vertex v) -> bool {
						double tmp = value<Vec3>(*geometric_mesh_, p.vertex_position_.get(), v).dot(axis_z);
						if (value<Vec3>(*geometric_mesh_, p.vertex_position_.get(), v).dot(axis_z) < d)
						{
							value<Vec3>(*geometric_mesh_, p.vertex_position_.get(), v) += (d - tmp) * axis_z;
							value<Vec3>(*geometric_mesh_, simu_solver.speed_.get(), v) -=
								axis_z.dot(value<Vec3>(*geometric_mesh_, simu_solver.speed_.get(), v)) * axis_z;
						}
						return true;
					});
				}
				map.end_writer();
				// std::this_thread::sleep_for(std::chrono::milliseconds(5));
			}
		});

		launch_thread([this, &p]() {
			cv_m.try_lock();
			while (this->running_)
			{
				meca_update_ = true;
				// cv_m.try_lock();
				std::this_thread::sleep_for(std::chrono::milliseconds(16));
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
		Parameters& p = parameters_[mecanical_mesh_];
		if (p.have_selected_vertex_)
		{
			Vec3 pos = value<Vec3>(*mecanical_mesh_, p.vertex_position_.get(), p.selected_vertex_);
			double m = value<double>(*mecanical_mesh_, p.vertex_masse_.get(), p.selected_vertex_);
			value<Vec3>(*mecanical_mesh_, p.vertex_forces_.get(), p.selected_vertex_) =
				m * (p.move_vertex_ - pos) / TIME_STEP;
		}

		/*parallel_foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool {
			double m = value<double>(*mecanical_mesh_, p.vertex_masse_.get(), v);
			value<Vec3>(*mecanical_mesh_, p.vertex_forces_.get(), v) += m * Vec3(0, 0, -9.81);
			return true;
		});*/

		simu_solver.compute_time_step(*geometric_mesh_, p.vertex_position_.get(), p.vertex_masse_.get(), TIME_STEP,
									  modif_topo_);
		if (p.show_frame_manipulator_)
		{
			Vec3 position;
			Vec3 axis_z;
			p.frame_manipulator_.get_position(position);
			p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, axis_z);
			double d = position.dot(axis_z);
			parallel_foreach_cell(mecanical_mesh_->m_, [&](Vertex v) -> bool {
				double tmp = value<Vec3>(*geometric_mesh_, p.vertex_position_.get(), v).dot(axis_z);
				if (value<Vec3>(*geometric_mesh_, p.vertex_position_.get(), v).dot(axis_z) < d)
				{
					value<Vec3>(*geometric_mesh_, p.vertex_position_.get(), v) += (d - tmp) * axis_z;
					value<Vec3>(*geometric_mesh_, simu_solver.speed_.get(), v) -=
						axis_z.dot(value<Vec3>(*geometric_mesh_, simu_solver.speed_.get(), v)) * axis_z;
				}
				return true;
			});
		}
		need_update_ = true;
	}

	void reset_forces()
	{
		simu_solver.reset_forces(*mecanical_mesh_);
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			MeshData<MR_MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.have_selected_vertex_ && p.param_move_vertex_->vao_initialized())
			{
				p.param_move_vertex_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
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

		if (ImGui::ListBoxHeader("Mesh_meca"))
		{
			mesh_provider_->foreach_mesh([this](MR_MESH* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == mecanical_mesh_))
				{
					mecanical_mesh_ = m;
				}
			});
			ImGui::ListBoxFooter();
		}
		if (ImGui::ListBoxHeader("Mesh_geom"))
		{
			mesh_provider_->foreach_mesh([this](MR_MESH* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == geometric_mesh_))
				{
					geometric_mesh_ = m;
					// simu_solver.init_solver(*geometric_mesh_, &sm_solver_, &ps_);
				}
			});
			ImGui::ListBoxFooter();
		}

		if (mecanical_mesh_ && geometric_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			Parameters& p = parameters_[mecanical_mesh_];

			need_update_ |= ImGui::Checkbox("Show ground", &p.show_frame_manipulator_);
			if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(
					*mecanical_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
						bool is_selected = attribute == p.vertex_position_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						{
							set_vertex_position(*mecanical_mesh_, attribute);
							simu_solver.init_solver(*mecanical_mesh_, &sm_solver_, p.vertex_position_.get(), &ps_);
							mesh_provider_->register_mesh(simu_solver.coarse_meca_mesh_, "coarse_mesh");
							mesh_provider_->register_mesh(simu_solver.fine_meca_mesh_, "fine_mesh");
							p.vertex_masse_ = sm_solver_.masse_;
							p.init_vertex_position_ = sm_solver_.vertex_init_position_;
							p.vertex_forces_ = simu_solver.forces_ext_;
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
				ImGui::EndCombo();
			}
			if (p.vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					set_vertex_position(*mecanical_mesh_, nullptr);
			}

			if (ImGui::BeginCombo("vertex relative position", p.vertex_relative_position_
																  ? p.vertex_relative_position_->name().c_str()
																  : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*mecanical_mesh_,
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == p.vertex_relative_position_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
													{
														set_vertex_relative_position(*mecanical_mesh_, attribute);
													}
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (p.vertex_relative_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##relative"))
					set_vertex_relative_position(*mecanical_mesh_, nullptr);
			}
			if (ImGui::BeginCombo("vertex parent",
								  p.vertex_parents_ ? p.vertex_parents_->name().c_str() : "-- select --"))
			{
				foreach_attribute<std::array<Vertex, 3>, Vertex>(
					*mecanical_mesh_, [&](const std::shared_ptr<Attribute<std::array<Vertex, 3>>>& attribute) {
						bool is_selected = attribute == p.vertex_parents_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
							set_vertex_parents(*mecanical_mesh_, attribute);
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
				ImGui::EndCombo();
			}
			if (p.vertex_parents_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##parents"))
					set_vertex_parents(*mecanical_mesh_, nullptr);
			}

			if (p.vertex_position_ && p.init_vertex_position_ && p.vertex_forces_)
			{
				ImGui::Separator();

				MeshData<MR_MESH>* md = mesh_provider_->mesh_data(mecanical_mesh_);
				Parameters& p = parameters_[mecanical_mesh_];

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
					p.update_move_vertex_vbo();
					mesh_provider_->foreach_mesh([&](MR_MESH* m, const std::string&) {
						mesh_provider_->emit_attribute_changed(m, p.vertex_position_.get());
						mesh_provider_->emit_attribute_changed(m, simu_solver.diff_volume_current_fine_.get());
						mesh_provider_->emit_attribute_changed(m, simu_solver.diff_volume_coarse_current_.get());
						mesh_provider_->emit_attribute_changed(m, simu_solver.pos_current_.get());
						mesh_provider_->emit_attribute_changed(m, simu_solver.pos_coarse_.get());
						if (modif_topo_)
						{
							mesh_provider_->emit_connectivity_changed(m);
						}
					});
					need_update_ = false;
					meca_update_ = false;
					modif_topo_ = false;
					// cv.notify_all();
					cv_m.unlock();
				}
				simulation::shape_matching_constraint_solver<MR_MESH>* sm1 =
					static_cast<simulation::shape_matching_constraint_solver<MR_MESH>*>(simu_solver.sc_);
				simulation::shape_matching_constraint_solver<MR_MESH>* sm2 =
					static_cast<simulation::shape_matching_constraint_solver<MR_MESH>*>(simu_solver.sc_fine_.get());
				simulation::shape_matching_constraint_solver<MR_MESH>* sm3 =
					static_cast<simulation::shape_matching_constraint_solver<MR_MESH>*>(simu_solver.sc_coarse_.get());

				double min = 0, max = 1, new_alpha = sm1->stiffness_;
				ImGui::SliderScalar("stiffness", ImGuiDataType_Double, &new_alpha, &min, &max);

				sm1->stiffness_ = new_alpha;
				sm2->stiffness_ = new_alpha;
				sm3->stiffness_ = new_alpha;
			}
		}
	}

public:
	MR_MESH* mecanical_mesh_;
	MR_MESH* geometric_mesh_;
	std::unordered_map<const MR_MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MR_MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MR_MESH>* mesh_provider_;
	simulation::shape_matching_constraint_solver<MR_MESH> sm_solver_;
	simulation::Simulation_solver_multiresolution<MR_MESH> simu_solver;
	simulation::Propagation_Plastique<MR_MESH> ps_;
	std::condition_variable cv;
	std::mutex cv_m;
	bool running_;
	bool need_update_;
	bool meca_update_;
	bool can_move_vertex_;
	bool modif_topo_;
	View* selected_view_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_ANIMATION_MULTIRESOLUTION_H_
