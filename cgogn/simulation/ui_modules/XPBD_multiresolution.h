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

#ifndef CGOGN_MODULE_XPBD_MULTIRESOLUTION_H_
#define CGOGN_MODULE_XPBD_MULTIRESOLUTION_H_

#include <GLFW/glfw3.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
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
#include <cgogn/simulation/algos/XPBD/XPBD_multiresolution.h>
#include <cgogn/simulation/algos/linked_volumes/linked_volumes.h>

#include <cgogn/modeling/algos/volume_utils.h>

#include <boost/synapse/connect.hpp>
#include <imgui/imgui.h>

#include <cgogn/rendering/shape_drawer.h>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class XPBD_Multiresolution_View : public ViewModule
{

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using SURFACE = CMap2;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Volume = typename mesh_traits<MESH>::Volume;
	using Face = typename mesh_traits<MESH>::Face;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), init_vertex_position_(nullptr), vertex_forces_(nullptr),
			  vertex_masse_(nullptr), fixed_vertex(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  have_selected_vertex_(false), move_vertex_(0, 0, 0), show_frame_manipulator_(false),
			  manipulating_frame_(false)
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

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> init_vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_forces_;
		std::shared_ptr<Attribute<double>> vertex_masse_;
		std::shared_ptr<Attribute<bool>> fixed_vertex;

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
	XPBD_Multiresolution_View(const App& app)
		: ViewModule(app, "XPBD (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  geom_mesh_(nullptr), selected_view_(app.current_view()), running_(false), apply_gravity(false),
		  take_screenshot_(false), ground_(false), inverse_control_(nullptr), draw_cylinder(false),
		  radius_cylinder(50.0f), pos_cylinder1(-272, 31, -79), Zaxis_cylinder1(1. / sqrt(2.), 1. / sqrt(2.), 0),
		  pos_cylinder2(-800, -1400, 5), Zaxis_cylinder2(0, 0, 1), pos_cylinder3(170, -2100, 5),
		  Zaxis_cylinder3(0, 0, 1), pos_sphere(700, 0, 100), pos_sphere2(600, 400, 400), shape_(nullptr),
		  show_sphere_(false), sphere_radius_(100.0f), gravity_intensity_(1.)
	{
		f_keypress = [](View*, MESH*, int32, CellsSet<MESH, Vertex>*, CellsSet<MESH, Edge>*) {};
	}

	~XPBD_Multiresolution_View()
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
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, m]() {
				Parameters& p = parameters_[m];
				if (p.vertex_position_ && p.init_vertex_position_ && p.vertex_forces_ && p.vertex_masse_)
				{
					// sm_solver_.update_topo(*m, {});
				}
			}));
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		simu_solver.init_solver(*selected_mesh_, vertex_position);
		p.vertex_position_ = vertex_position;
		p.vertex_forces_ = simu_solver.f_ext_;
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
		simu_solver.f_ext_ = vertex_forces;
	}
	void set_vertex_masse(const MESH& m, const std::shared_ptr<Attribute<double>>& vertex_masse)
	{
		Parameters& p = parameters_[&m];

		p.vertex_masse_ = vertex_masse;
	}

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &XPBD_Multiresolution_View<MESH>::init_mesh));
		shape_ = rendering::ShapeDrawer::instance();
		shape_->color(rendering::ShapeDrawer::CYLINDER) = rendering::GLColor(0.5294, 0.6078, 0.6078, 1);
		shape_->color(rendering::ShapeDrawer::CUBE) = rendering::GLColor(0.5294, 0.6078, 0.6078, 1);
		shape_->color(rendering::ShapeDrawer::SPHERE) = rendering::GLColor(0.88, 0.55, 0.17, 1);
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
		if (p.manipulating_frame_)
		{
			auto [P, Q] = view->pixel_ray(x, y);
			p.frame_manipulator_.pick(x, y, P, Q);
			view->request_update();
		}
	}
	void key_press_event(View* v, int32 key_code) override
	{
		if (key_code == GLFW_KEY_LEFT_CONTROL)
		{
			v->lock_rotation_ = true;
			can_move_vertex_ = true;
		}
		if (key_code == GLFW_KEY_LEFT_SHIFT)
		{
			inverse_control_ = true;
		}
		if (key_code == GLFW_KEY_C)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				if (p.show_frame_manipulator_)
					p.manipulating_frame_ = true;
			}
		}
		if (key_code == GLFW_KEY_G)
		{
			apply_gravity = !apply_gravity;
		}
		if (key_code == GLFW_KEY_P)
		{
			ground_ = !ground_;
			std::cout << ground_ << std::endl;
		}
		if (key_code == GLFW_KEY_S)
		{
			take_screenshot_ = !take_screenshot_;
			frame_number_ = 0;
		}
		if (key_code == GLFW_KEY_F)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				Vec3 pos;
				p.frame_manipulator_.get_position(pos);
				Vec3 a;
				p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, a);
				std::cout << "pos : " << pos << std::endl;
				std::cout << " normale : " << a << std::endl;
				double d = pos.dot(a);
				parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
					if (value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a) < d)
					{
						value<bool>(*selected_mesh_, p.fixed_vertex.get(), v) = true;
					}
					return true;
				});
			}
		}
		if (key_code == GLFW_KEY_M)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				Vec3 pos;
				p.frame_manipulator_.get_position(pos);
				Vec3 a;
				p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, a);
				std::cout << "pos : " << pos << std::endl;
				std::cout << " normale : " << a << std::endl;
				double d = pos.dot(a);
				moving_vertices.clear();
				parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
					if (value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a) < d)
					{
						value<bool>(*selected_mesh_, p.fixed_vertex.get(), v) = true;
						moving_vertices.push_back(v);
					}
					return true;
				});
			}
		}
		if (key_code == GLFW_KEY_B)
		{
			draw_cylinder = !draw_cylinder;
			v->request_update();
		}
		if (key_code == GLFW_KEY_K)
		{
			draw_sphere2 = !draw_sphere2;
			v->request_update();
		}
		if (key_code == GLFW_KEY_E)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				Vec3 pos1(110, 0, 0);
				Vec3 a1(1, 0, 0);
				double d1 = pos1.dot(a1);
				Vec3 pos2(1210, 0, 0);
				Vec3 a2(-1, 0, 0);
				double d2 = pos2.dot(a2);
				Vec3 pos3(0, 60, 0);
				Vec3 a3(0, 1, 0);
				double d3 = pos3.dot(a3);
				Vec3 pos4(0, 680, 0);
				Vec3 a4(0, -1, 0);
				double d4 = pos4.dot(a4);
				parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
					if (value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a1) < d1 ||
						value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a2) < d2 ||
						value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a3) < d3 ||
						value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a4) < d4)
					{
						value<bool>(*selected_mesh_, p.fixed_vertex.get(), v) = true;
					}
					/* if (value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(a2) < d2)
						moving_vertices.push_back(v);*/
					return true;
				});
			}
		}
		if (key_code == GLFW_KEY_U)
		{
			if (selected_mesh_)
			{
				// std::clock_t start = std::clock();
				// double duration = 0;
				Parameters& p = parameters_[selected_mesh_];
				for (int i = 0; i < 50; i++)
				{
					simu_solver.solver(*selected_mesh_, geom_mesh_, 0.01f, false);
					simu_solver.compute_contact(*selected_mesh_, *geom_mesh_, [&](Vertex v) -> bool {
						Vec3& pos = value<Vec3>(*geom_mesh_, p.vertex_position_.get(), v);
						Vec3 axis_z = Zaxis_cylinder1;

						Vec3 pos2 = pos - pos_cylinder1.cast<double>();

						double dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							return true;
						}
						return false;
					});
					parallel_foreach_cell(*geom_mesh_, [&](Vertex v) -> bool {
						Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
						Vec3& speed = value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v);
						Vec3 axis_z = Zaxis_cylinder1;

						Vec3 pos2 = pos - pos_cylinder1.cast<double>();

						double dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							Vec3 dir_collision = (axis_z * (axis_z.dot(pos2)) - pos2).normalized() * dist;
							pos += dir_collision;

							Vec3 dir_col_norm = dir_collision.normalized();
							double tmp = speed.dot(-dir_col_norm);
							if (tmp > 0)
							{
								speed += dir_col_norm * tmp;
							}
							return true;
						}
						return true;
					});
					pos_cylinder1 = Eigen::Vector3f(260, 370 + (cos(2 * M_PI / 100 * i) + 1) / 2.0 * 300.0 - 300.0, 0);
				}
				/*duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				std::cout << "temps solve xpbd : " << duration / 100.0f << std::endl;*/
				need_update_ = true;
			}
		}
		if (key_code == GLFW_KEY_V)
		{
			Parameters& p = parameters_[selected_mesh_];
			for (int i = 0; i < 500; i++)
			{
				static Vec3 cm = geometry::centroid<Vec3>(*selected_mesh_, p.vertex_position_.get());
				Eigen::Affine3f transfo = Eigen::Translation3f(cm.cast<float>()) *
										  Eigen::AngleAxisf(0.01, Eigen::Vector3f::UnitZ()) *
										  Eigen::Translation3f(-cm.cast<float>());
				pos_sphere = transfo * pos_sphere;
				simu_solver.solver(*selected_mesh_, geom_mesh_, 0.01f, false);
				simu_solver.compute_error_point(*selected_mesh_, pos_sphere.cast<double>(),
												4 + 2 * cos(float(i) / 10.0f));
			}
			/*duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "temps solve xpbd : " << duration / 100.0f << std::endl;*/
			need_update_ = true;
		}
		if (key_code == GLFW_KEY_X)
		{
			if (inverse_control_)
			{
				pos_cylinder1 -= Eigen::Vector3f(0.1, 0, 0);
			}
			else
			{
				pos_cylinder1 += Eigen::Vector3f(0.1, 0, 0);
			}
			std::cout << pos_cylinder1 << std::endl;
			v->request_update();
		}
		if (key_code == GLFW_KEY_Z)
		{
			if (inverse_control_)
			{
				pos_cylinder1 -= Eigen::Vector3f(0, 0, 0.1);
			}
			else
			{
				pos_cylinder1 += Eigen::Vector3f(0, 0, 0.1);
			}

			v->request_update();
		}
		if (key_code == GLFW_KEY_N)
		{
			show_sphere_ = !show_sphere_;

			v->request_update();
		}
	}

	void key_release_event(View* v, int32 key_code)
	{
		if (key_code == GLFW_KEY_LEFT_CONTROL)
		{
			v->lock_rotation_ = false;
			can_move_vertex_ = false;
		}
		if (key_code == GLFW_KEY_LEFT_SHIFT)
		{
			inverse_control_ = false;
		}
		if (key_code == GLFW_KEY_C)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				p.manipulating_frame_ = false;
			}
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y)
	{
		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];
			if (p.have_selected_vertex_ && can_move_vertex_)
			{
				p.move_vertex_ = view->pixel_scene_(
					x, y, value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_));
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

#define TIME_STEP 0.001666f
	void start()
	{
		running_ = true;

		launch_thread([this]() {
			while (this->running_)
			{
				Parameters& p = parameters_[selected_mesh_];
				if (p.have_selected_vertex_)
				{
					Vec3 pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_);
					double m = value<double>(*selected_mesh_, p.vertex_masse_.get(), p.selected_vertex_);
					value<Vec3>(*selected_mesh_, p.vertex_forces_.get(), p.selected_vertex_) =
						m * (p.move_vertex_ - pos) / TIME_STEP;
					std::cout << value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_)
							  << std::endl;
				}
				for (Vertex v : moving_vertices)
				{
					value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v) += Vec3(0.1, 0, 0);
				}
				selected_mesh_->start_writer();
				
				if (show_sphere_)
				{
					static Vec3 cm = geometry::centroid<Vec3>(*selected_mesh_, p.vertex_position_.get());
					Eigen::Affine3f transfo = Eigen::Translation3f(cm.cast<float>()) *
											  Eigen::AngleAxisf(0.01, Eigen::Vector3f::UnitZ()) *
											  Eigen::Translation3f(-cm.cast<float>());
					pos_sphere = transfo * pos_sphere;
				}
				for (int i = 0; i < 1; i++)
				{

					if (apply_gravity)
					{
						parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
							value<Vec3>(*selected_mesh_, p.vertex_forces_, v) +=
								value<double>(*selected_mesh_, simu_solver.masse_, v) * Vec3(0, -98.1, 0) *
								gravity_intensity_;
							return true;
						});
					}
					simu_solver.solver(*selected_mesh_, geom_mesh_, TIME_STEP, false);
					parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						value<Vec3>(*selected_mesh_, p.vertex_forces_, v) = Vec3(0, 0, 0);
						return true;
					});
				}
				if (!moving_vertices.empty())
				{
					parallel_foreach_cell(*selected_mesh_, [&](Volume v) -> bool {
						std::shared_ptr<Attribute<Eigen::Matrix3d>> F =
							get_attribute<Eigen::Matrix3d, Volume>(*selected_mesh_, "XPBD_F_volume");
						Eigen::Matrix3d deformation_gradient = value<Eigen::Matrix3d>(*selected_mesh_, F, v);
						Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(deformation_gradient);
						auto eigenvalue = eigensolver.eigenvalues();
						int max_id = 0;
						for (int i = 1; i < 3; i++)
						{
							if (eigenvalue(i) > eigenvalue(max_id))
							{
								max_id = i;
							}
						}
						std::shared_ptr<Attribute<double>> vol_eigen =
							get_or_add_attribute<double, Volume>(*selected_mesh_, "Volume_eigen_value");
						std::shared_ptr<Attribute<Vec3>> vol_eigen_vector =
							get_or_add_attribute<Vec3, Volume>(*selected_mesh_, "Volume_eigen_vector");
						value<double>(*selected_mesh_, vol_eigen, v) = eigenvalue(max_id);
						value<Vec3>(*selected_mesh_, vol_eigen_vector, v) = eigensolver.eigenvectors().col(max_id);
						std::shared_ptr<Attribute<double>> separation_volume =
							get_or_add_attribute<double, Volume>(*selected_mesh_, "Separation_volume");
						value<double>(*selected_mesh_, separation_volume, v) = 0;
						return true;
					});
					std::vector<Face> face_unsew;
					parallel_foreach_cell(*selected_mesh_, [&](Face f) -> bool {
						if (is_incident_to_boundary(*selected_mesh_, f))
							return true;
						std::shared_ptr<Attribute<double>> vol_eigen =
							get_attribute<double, Volume>(*selected_mesh_, "Volume_eigen_value");
						std::shared_ptr<Attribute<Vec3>> vol_eigen_vector =
							get_attribute<Vec3, Volume>(*selected_mesh_, "Volume_eigen_vector");
						std::shared_ptr<Attribute<double>> face_tensor =
							get_or_add_attribute<double, Face>(*selected_mesh_, "Face_tensor");
						double evalue_1 = value<double>(*selected_mesh_, vol_eigen, Volume(f.dart));
						double evalue_2 =
							value<double>(*selected_mesh_, vol_eigen, Volume(phi3(*selected_mesh_, f.dart)));
						Vec3 evector_1 = value<Vec3>(*selected_mesh_, vol_eigen_vector, Volume(f.dart));
						Vec3 evector_2 =
							value<Vec3>(*selected_mesh_, vol_eigen_vector, Volume(phi3(*selected_mesh_, f.dart)));
						Vec3 face_normale = geometry::normal(*selected_mesh_, f, p.vertex_position_.get());
						double tension = (evalue_1 - 1.) * abs(evector_1.normalized().dot(face_normale)) +
										 (evalue_2 - 1.) * abs(evector_2.normalized().dot(face_normale));
						value<double>(*selected_mesh_, face_tensor, f) = tension;
						if (tension > 0.3f)
						{
							std::shared_ptr<Attribute<double>> separation_volume =
								get_attribute<double, Volume>(*selected_mesh_, "Separation_volume");
							value<double>(*selected_mesh_, separation_volume, Volume(f.dart)) += 1;
							value<double>(*selected_mesh_, separation_volume, Volume(phi3(*selected_mesh_, f.dart))) +=
								1;
							face_unsew.push_back(f);
						}
						return true;
					});
					if (!face_unsew.empty())
					{

						for (Face f : face_unsew)
						{
							unsew_volume(
								*selected_mesh_, f,
								[&](std::pair<Vertex, Vertex> pv) -> bool {
									std::vector<std::shared_ptr<Attribute<Vec3>>> list_update_attribute;
									list_update_attribute.push_back(simu_solver.pos_);
									list_update_attribute.push_back(simu_solver.init_pos_);
									list_update_attribute.push_back(simu_solver.speed_);
									list_update_attribute.push_back(simu_solver.f_ext_);
									for (auto attr : list_update_attribute)
									{
										value<Vec3>(*selected_mesh_, attr, pv.second) =
											value<Vec3>(*selected_mesh_, attr, pv.first);
									}
									return true;
								},
								true);
						}
						
						foreach_cell(*selected_mesh_, [&](Volume vol) -> bool {
							Vec3 cm = value<Vec3>(*selected_mesh_, simu_solver.centroid_, vol);
							geometry::Mat3d F = value<geometry::Mat3d>(*selected_mesh_, simu_solver.F_, vol);

							foreach_incident_vertex(*selected_mesh_, vol, [&](Vertex w) -> bool {
								Vec3 init_r_i = value<Vec3>(*selected_mesh_, simu_solver.init_pos_.get(), w) -
												value<Vec3>(*selected_mesh_, simu_solver.init_cm_, vol);
								value<Vec3>(*selected_mesh_, simu_solver.pos_.get(), w) = cm + F * init_r_i;
								return true;
							});
							return true;
						});
						simu_solver.update_topo(*selected_mesh_);
					}
				}

				if (show_sphere_)
				{

					static Vec3 cm = geometry::centroid<Vec3>(*selected_mesh_, p.vertex_position_.get());
					Eigen::Affine3f transfo = Eigen::Translation3f(cm.cast<float>()) *
											  Eigen::AngleAxisf(0.01, Eigen::Vector3f::UnitZ()) *
											  Eigen::Translation3f(-cm.cast<float>());
					pos_sphere = transfo * pos_sphere;
				}
				for (int i = 0; i < 1; i++)
				{

					if (apply_gravity)
					{
						parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
							value<Vec3>(*selected_mesh_, p.vertex_forces_, v) +=
								value<double>(*selected_mesh_, simu_solver.masse_, v) * Vec3(0, 0, -98.1) *
								gravity_intensity_;
							return true;
						});
					}
					simu_solver.solver(*selected_mesh_, geom_mesh_, TIME_STEP, false);
					parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						value<Vec3>(*selected_mesh_, p.vertex_forces_, v) = Vec3(0, 0, 0);
						return true;
					});
				}
				if (show_sphere_)
				{
					static double it_sphere = 0;
					it_sphere += 0.1;
					simu_solver.compute_error_point(*selected_mesh_, pos_sphere.cast<double>(), 4 + 2 * cos(it_sphere));
				}
				if (draw_sphere2)
				{
					simu_solver.compute_contact(*selected_mesh_, *geom_mesh_, [&](Vertex v) -> bool {
						Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);

						Vec3 pos2 = pos - pos_sphere2.cast<double>();

						double dist = pos2.norm() - sphere_radius2_;

						if (dist < 0)
						{
							return true;
						}
						return false;
					});
					parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
						Vec3& speed = value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v);

						Vec3 pos2 = pos - pos_sphere2.cast<double>();

						double dist = pos2.norm() - sphere_radius2_;

						if (dist < 0)
						{
							Vec3 dir_collision = -pos2.normalized() * dist;
							pos += dir_collision;

							Vec3 dir_col_norm = dir_collision.normalized();
							double tmp = speed.dot(-dir_col_norm);
							if (tmp > 0)
							{
								speed += dir_col_norm * tmp;
							}
							return true;
						}
						return true;
					});
					pos_sphere2 -= Eigen::Vector3f(0, 0, 1);
					parallel_foreach_cell(*selected_mesh_, [&](Volume v) -> bool {
						std::shared_ptr<Attribute<Eigen::Matrix3d>> F =
							get_attribute<Eigen::Matrix3d, Volume>(*selected_mesh_, "XPBD_F_volume");
						Eigen::Matrix3d deformation_gradient = value<Eigen::Matrix3d>(*selected_mesh_, F, v);
						Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(deformation_gradient);
						auto eigenvalue = eigensolver.eigenvalues();
						int max_id = 0;
						for (int i = 1; i < 3; i++)
						{
							if (eigenvalue(i) > eigenvalue(max_id))
							{
								max_id = i;
							}
						}
						std::shared_ptr<Attribute<double>> vol_eigen =
							get_or_add_attribute<double, Volume>(*selected_mesh_, "Volume_eigen_value");
						std::shared_ptr<Attribute<Vec3>> vol_eigen_vector =
							get_or_add_attribute<Vec3, Volume>(*selected_mesh_, "Volume_eigen_vector");
						value<double>(*selected_mesh_, vol_eigen, v) = eigenvalue(max_id);
						value<Vec3>(*selected_mesh_, vol_eigen_vector, v) = eigensolver.eigenvectors().col(max_id);
						std::shared_ptr<Attribute<double>> separation_volume =
							get_or_add_attribute<double, Volume>(*selected_mesh_, "Separation_volume");
						value<double>(*selected_mesh_, separation_volume, v) = 0;
						return true;
					});
					std::vector<Face> face_unsew;
					parallel_foreach_cell(*selected_mesh_, [&](Face f) -> bool {
						if (is_incident_to_boundary(*selected_mesh_, f))
							return true;
						std::shared_ptr<Attribute<double>> vol_eigen =
							get_attribute<double, Volume>(*selected_mesh_, "Volume_eigen_value");
						std::shared_ptr<Attribute<Vec3>> vol_eigen_vector =
							get_attribute<Vec3, Volume>(*selected_mesh_, "Volume_eigen_vector");
						std::shared_ptr<Attribute<double>> face_tensor =
							get_or_add_attribute<double, Face>(*selected_mesh_, "Face_tensor");
						double evalue_1 = value<double>(*selected_mesh_, vol_eigen, Volume(f.dart));
						double evalue_2 =
							value<double>(*selected_mesh_, vol_eigen, Volume(phi3(*selected_mesh_, f.dart)));
						Vec3 evector_1 = value<Vec3>(*selected_mesh_, vol_eigen_vector, Volume(f.dart));
						Vec3 evector_2 =
							value<Vec3>(*selected_mesh_, vol_eigen_vector, Volume(phi3(*selected_mesh_, f.dart)));
						Vec3 face_normale = geometry::normal(*selected_mesh_, f, p.vertex_position_.get());
						double tension = (evalue_1 - 1.) * abs(evector_1.normalized().dot(face_normale)) +
										 (evalue_2 - 1.) * abs(evector_2.normalized().dot(face_normale));
						value<double>(*selected_mesh_, face_tensor, f) = tension;
						if (tension > 0.5f)
						{
							std::shared_ptr<Attribute<double>> separation_volume =
								get_attribute<double, Volume>(*selected_mesh_, "Separation_volume");
							value<double>(*selected_mesh_, separation_volume, Volume(f.dart)) += 1;
							value<double>(*selected_mesh_, separation_volume, Volume(phi3(*selected_mesh_, f.dart))) +=
								1;
							face_unsew.push_back(f);
						}
						return true;
					});
					if (!face_unsew.empty())
					{

						for (Face f : face_unsew)
						{
							unsew_volume(
								*selected_mesh_, f,
								[&](std::pair<Vertex, Vertex> p) -> bool {
									std::vector<std::shared_ptr<Attribute<Vec3>>> list_update_attribute;
									list_update_attribute.push_back(simu_solver.pos_);
									list_update_attribute.push_back(simu_solver.init_pos_);
									list_update_attribute.push_back(simu_solver.speed_);
									for (auto attr : list_update_attribute)
									{
										value<Vec3>(*selected_mesh_, attr, p.second) =
											value<Vec3>(*selected_mesh_, attr, p.first);
									}
									return true;
								},
								true);
						}
						simu_solver.update_topo(*selected_mesh_);
					}
				}
				if (draw_cylinder)
				{
#if 1
					simu_solver.compute_contact(*selected_mesh_, *geom_mesh_, [&](Vertex v) -> bool {
						Vec3& pos = value<Vec3>(*geom_mesh_, p.vertex_position_.get(), v);
						Vec3 axis_z = Zaxis_cylinder1;

						Vec3 pos2 = pos - pos_cylinder1.cast<double>();

						double dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							return true;
						}
						return false;
						axis_z = Zaxis_cylinder2;

						pos2 = pos - pos_cylinder2.cast<double>();

						dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							return true;
						}
						axis_z = Zaxis_cylinder3;

						pos2 = pos - pos_cylinder3.cast<double>();

						dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							return true;
						}
						return false;
					});
#endif
					parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
						Vec3& speed = value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v);
						Vec3 axis_z = Zaxis_cylinder1;

						Vec3 pos2 = pos - pos_cylinder1.cast<double>();

						double dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							Vec3 dir_collision = (axis_z * (axis_z.dot(pos2)) - pos2).normalized() * dist;
							pos += dir_collision;

							Vec3 dir_col_norm = dir_collision.normalized();
							double tmp = speed.dot(-dir_col_norm);
							if (tmp > 0)
							{
								speed += dir_col_norm * tmp;
							}
							return true;
						}
						return false;
						axis_z = Zaxis_cylinder2;

						pos2 = pos - pos_cylinder2.cast<double>();

						dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							Vec3 dir_collision = (axis_z * (axis_z.dot(pos2)) - pos2).normalized() * dist;
							pos += dir_collision;

							Vec3 dir_col_norm = dir_collision.normalized();
							double tmp = speed.dot(-dir_col_norm);
							if (tmp > 0)
							{
								speed += dir_col_norm * tmp;
							}
							return true;
						}
						axis_z = Zaxis_cylinder3;

						pos2 = pos - pos_cylinder3.cast<double>();

						dist = axis_z.cross(pos2).norm() - radius_cylinder;

						if (dist < 0)
						{
							Vec3 dir_collision = (axis_z * (axis_z.dot(pos2)) - pos2).normalized() * dist;
							pos += dir_collision;

							Vec3 dir_col_norm = dir_collision.normalized();
							double tmp = speed.dot(-dir_col_norm);
							if (tmp > 0)
							{
								speed += dir_col_norm * tmp;
							}
							return true;
						}
						return true;
					});
					/*static int nb_iter = 0;
					pos_cylinder1 =
						Eigen::Vector3f(260, 370 + (cos(2 * M_PI / 1000 * nb_iter) + 1) / 2.0 * 200.0 - 200.0, 0);
					nb_iter++;*/
					// pos_cylinder1 += Eigen::Vector3f(0, 1., 0);
				}
				if (ground_)
				{
					Vec3 position;
					Vec3 axis_z;
					p.frame_manipulator_.get_position(position);
					p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, axis_z);
					double d = position.dot(axis_z);
					simu_solver.compute_contact(*selected_mesh_, *geom_mesh_, [&](Vertex v) -> bool {
						double tmp = value<Vec3>(*geom_mesh_, p.vertex_position_.get(), v).dot(axis_z);
						if (tmp < d)
						{
							return true;
						}
						return false;
					});

					parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
						double tmp = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(axis_z);
						if (value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v).dot(axis_z) < d)
						{
							value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v) += (d - tmp) * axis_z;
							value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v) -=
								axis_z.dot(value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v)) * axis_z;
						}
						return true;
					});
				}
				need_update_ = true;
				selected_mesh_->end_writer();
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
		if (p.have_selected_vertex_)
		{
			Vec3 pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_);
			double m = value<double>(*selected_mesh_, p.vertex_masse_.get(), p.selected_vertex_);
			value<Vec3>(*selected_mesh_, p.vertex_forces_.get(), p.selected_vertex_) =
				m * (p.move_vertex_ - pos) / TIME_STEP;
			std::cout << value<Vec3>(*selected_mesh_, p.vertex_position_.get(), p.selected_vertex_) << std::endl;
		}
		for (int i = 0; i < 1; i++)
		{
			if (apply_gravity)
			{
				parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
					value<Vec3>(*selected_mesh_, p.vertex_forces_, v) +=
						value<double>(*selected_mesh_, simu_solver.masse_, v) * Vec3(0, -9.81, 0);
					return true;
				});
			}
			simu_solver.solver(*selected_mesh_, geom_mesh_, TIME_STEP, false);
			parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
				value<Vec3>(*selected_mesh_, p.vertex_forces_, v) = Vec3(0, 0, 0);
				return true;
			});
		}
		need_update_ = true;
	}

	void refresh_volume_skin()
	{
		if (!geom_mesh_)
			return;
		if (!volume_skin_)
			volume_skin_ = surface_provider_->add_mesh("volume_skin");
		Parameters& p = parameters_[selected_mesh_];

		surface_provider_->clear_mesh(*volume_skin_);

		volume_skin_vertex_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "position");
		volume_skin_vertex_normal_ = get_or_add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "normal");
		volume_skin_vertex_index_ = get_or_add_attribute<uint32, SurfaceVertex>(*volume_skin_, "vertex_index");
		volume_skin_vertex_volume_vertex_ = get_or_add_attribute<Vertex, SurfaceVertex>(*volume_skin_, "hex_vertex");
		modeling::extract_volume_surface(*geom_mesh_, p.vertex_position_.get(), *volume_skin_,
										 volume_skin_vertex_position_.get(), volume_skin_vertex_volume_vertex_.get());

		uint32 nb_vertices = 0;
		foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			value<uint32>(*volume_skin_, volume_skin_vertex_index_, v) = nb_vertices++;
			return true;
		});

		surface_provider_->emit_connectivity_changed(*volume_skin_);
	}

	void draw(View* view) override
	{
		const rendering::GLMat4& proj_matrix = view->projection_matrix();
		const rendering::GLMat4& view_matrix = view->modelview_matrix();
		if (selected_mesh_)
		{
			auto& m = selected_mesh_;
			auto& p = parameters_[selected_mesh_];

			MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

			if (p.have_selected_vertex_ && p.param_move_vertex_->attributes_initialized())
			{
				p.param_move_vertex_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_move_vertex_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, 2);
				p.param_move_vertex_->release();
			}

			if (p.have_selected_vertex_ && p.param_edge_->attributes_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_LINES, 0, 2);
				p.param_edge_->release();
			}

			if (p.show_frame_manipulator_)
			{
				double size = (md.bb_max_ - md.bb_min_).norm() / 10;
				p.frame_manipulator_.set_size(size);
				p.frame_manipulator_.draw(true, true, proj_matrix, view_matrix);
			}
		}
		if (draw_cylinder)
		{
			Eigen::Affine3f transfo = Eigen::Translation3f(pos_cylinder1) *
									  Eigen::AngleAxisf(std::acos(Zaxis_cylinder1.x()), Eigen::Vector3f::UnitZ()) *
									  Eigen::AngleAxisf(std::acos(Zaxis_cylinder1.z()), Eigen::Vector3f::UnitY()) *
									  Eigen::Scaling(radius_cylinder, radius_cylinder, 1000.0f);
			shape_->draw(rendering::ShapeDrawer::CYLINDER, proj_matrix, view_matrix * transfo.matrix());
			transfo = Eigen::Translation3f(pos_cylinder2) *
					  Eigen::AngleAxisf(std::acos(Zaxis_cylinder2.x()), Eigen::Vector3f::UnitZ()) *
					  Eigen::AngleAxisf(std::acos(Zaxis_cylinder2.z()), Eigen::Vector3f::UnitY()) *
					  Eigen::Scaling(radius_cylinder, radius_cylinder, 1000.0f);
			// shape_->draw(rendering::ShapeDrawer::CYLINDER, proj_matrix, view_matrix * transfo.matrix());
			transfo = Eigen::Translation3f(pos_cylinder3) *
					  Eigen::AngleAxisf(std::acos(Zaxis_cylinder3.x()), Eigen::Vector3f::UnitZ()) *
					  Eigen::AngleAxisf(std::acos(Zaxis_cylinder3.z()), Eigen::Vector3f::UnitY()) *
					  Eigen::Scaling(radius_cylinder, radius_cylinder, 1000.0f);
			// shape_->draw(rendering::ShapeDrawer::CYLINDER, proj_matrix, view_matrix * transfo.matrix());
		}
		if (show_sphere_)
		{
			Eigen::Affine3f transfo = Eigen::Translation3f(pos_sphere) * Eigen::Scaling(sphere_radius_);
			shape_->draw(rendering::ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());
		}
		if (draw_sphere2)
		{
			Eigen::Affine3f transfo = Eigen::Translation3f(pos_sphere2) * Eigen::Scaling(sphere_radius2_);
			shape_->draw(rendering::ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());
		}
	}

	void left_panel() override
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
			mesh_provider_->foreach_mesh([this](MESH& m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), &m == selected_mesh_))
				{
					selected_mesh_ = &m;
					Parameters& p = parameters_[selected_mesh_];
					p.fixed_vertex = get_attribute<bool, Vertex>(m, "fixed_vertex");
					if (p.fixed_vertex == nullptr)
						p.fixed_vertex = add_attribute<bool, Vertex>(m, "fixed_vertex");
					simu_solver.fixed_vertex = p.fixed_vertex;
				}
			});
			ImGui::ListBoxFooter();
		}
		if (ImGui::ListBoxHeader("Geometry Mesh"))
		{
			mesh_provider_->foreach_mesh([this](MESH& m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), &m == geom_mesh_))
				{
					geom_mesh_ = &m;
				}
			});
			ImGui::ListBoxFooter();
		}

		if (ImGui::ListBoxHeader("Surface Mesh"))
		{
			surface_provider_->foreach_mesh([this](SURFACE& m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), &m == volume_skin_))
				{
					volume_skin_ = &m;
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
				std::shared_ptr<Attribute<Vec3>> attr = nullptr;
				foreach_attribute<Vec3, Vertex>(*selected_mesh_,
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == p.vertex_position_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														attr = attribute;

													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
				if (attr != nullptr)
				{
					set_vertex_position(*selected_mesh_, attr);
				}
			}
			if (p.vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					set_vertex_position(*selected_mesh_, nullptr);
			}
			if (p.vertex_position_)
			{
				ImGui::Separator();

				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
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
				double min = 0, max = 100;
				ImGui::SliderScalar("gravity", ImGuiDataType_Double, &gravity_intensity_, &min, &max);
				if (need_update_)
				{
					p.update_move_vertex_vbo();

					selected_mesh_->start_reader();
					std::shared_ptr<Attribute<double>> separation_volume =
						get_or_add_attribute<double, Volume>(*selected_mesh_, "Separation_volume");
					mesh_provider_->emit_attribute_changed(*selected_mesh_, separation_volume.get());
					mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
					mesh_provider_->emit_attribute_changed(*selected_mesh_, simu_solver.Det_F_Volume_.get());
					mesh_provider_->emit_connectivity_changed(*selected_mesh_);
					mesh_provider_->emit_connectivity_changed(*selected_mesh_->topology_);
					if (geom_mesh_)
						mesh_provider_->emit_attribute_changed(*geom_mesh_, p.vertex_position_.get());

					if (volume_skin_)
					{
						refresh_volume_skin();
						surface_provider_->emit_attribute_changed(*volume_skin_, volume_skin_vertex_position_.get());
					}

					if (take_screenshot_)
					{
						if (frame_number_ % 1 == 0)
						{
							// mesh_provider_->emit_connectivity_changed(*selected_mesh_);
							selected_view_->save_screenshot("../../../screen_video/screen" +
															std::to_string(frame_number_ / 1) + ".jpg");
						}
						frame_number_++;
					}
					selected_mesh_->end_reader();
					need_update_ = false;
				}
			}
		}
	}

public:
	MESH* selected_mesh_;
	MESH* geom_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
	cgogn::simulation::XPBD_Multiresolution simu_solver;
	bool running_;
	bool need_update_;
	bool can_move_vertex_;
	bool apply_gravity;
	bool take_screenshot_;
	int frame_number_;
	bool ground_;
	View* selected_view_;
	bool inverse_control_;
	bool draw_cylinder;
	Eigen::Vector3f pos_sphere2;
	bool draw_sphere2 = false;
	float sphere_radius2_ = 150.0f;
	float radius_cylinder;
	Eigen::Vector3f pos_cylinder1;
	Vec3 Zaxis_cylinder1;
	Eigen::Vector3f pos_cylinder2;
	Vec3 Zaxis_cylinder2;
	Eigen::Vector3f pos_cylinder3;
	Vec3 Zaxis_cylinder3;
	Eigen::Vector3f pos_sphere;
	bool show_sphere_;
	float sphere_radius_;
	rendering::ShapeDrawer* shape_;
	double gravity_intensity_;
	std::vector<Vertex> moving_vertices;

	SURFACE* volume_skin_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_position_ = nullptr;
	std::shared_ptr<SurfaceAttribute<uint32>> volume_skin_vertex_index_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_normal_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vertex>> volume_skin_vertex_volume_vertex_ = nullptr;
	ui::MeshProvider<SURFACE>* surface_provider_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_XPBD_MULTIRESOLUTION_H_
