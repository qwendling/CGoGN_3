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

#ifndef CGOGN_MODULE_SURFACE_RENDER_H_
#define CGOGN_MODULE_SURFACE_RENDER_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

class App;
class View;
template <typename MESH>
class MeshProvider;

template <typename MESH>
class SurfaceRender : public Module
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters() :
			initialized_(false),
			render_vertices_(false),
			render_edges_(false),
			render_faces_(true),
			phong_shading_(false),
			vertex_scale_factor_(1.0)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 1);

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 1, 0, 1);
			param_edge_->width_= 2.5f;

			param_flat_ =  rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(0, 0.8f, 0, 1);
			param_flat_->back_color_ = rendering::GLColor(0, 0, 0.8f, 1);
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_phong_ = rendering::ShaderPhong::generate_param();
			param_phong_->front_color_ = rendering::GLColor(0, 0.8f, 0, 1);
			param_phong_->back_color_ = rendering::GLColor(0, 0, 0.8f, 1);
			param_phong_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_phong_->specular_coef_ = 250.0f;
		}

		bool initialized_;

		std::string vertex_position_name_;
		std::string vertex_normal_name_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderPhong::Param> param_phong_;

		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool phong_shading_;
		
		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

public:

	SurfaceRender(const App& app) :
		ui::Module(app, "SurfaceRender (" + mesh_traits<MESH>::name + ")"),
		selected_mesh_(nullptr)
	{}
	~SurfaceRender()
	{}

	void init()
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + mesh_traits<MESH>::name + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
		{
			parameters_.emplace(m, Parameters());
		});
	}

	void update_data(const MESH& m)
	{
		Parameters& p = parameters_[&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		md->update_vbo(p.vertex_position_name_);
		if (!p.vertex_normal_name_.empty())
			md->update_vbo(p.vertex_normal_name_);

		std::shared_ptr<Attribute<Vec3>> vertex_position = get_attribute<Vec3, Vertex>(m, p.vertex_position_name_);
		p.vertex_base_size_ = geometry::mean_edge_length(m, vertex_position.get()) / 7.0;

		md->update_bb(vertex_position.get());

		p.param_point_sprite_->set_vbos(md->vbo(p.vertex_position_name_));
		p.param_edge_->set_vbos(md->vbo(p.vertex_position_name_));
		p.param_flat_->set_vbos(md->vbo(p.vertex_position_name_));
		p.param_phong_->set_vbos(md->vbo(p.vertex_position_name_), md->vbo(p.vertex_normal_name_));

		p.initialized_ = true;
	}

	void set_vertex_position(const MESH& m, const std::string& vertex_position_name)
	{
		parameters_[&m].vertex_position_name_ = vertex_position_name;
	}

	void set_vertex_normal(const MESH& m, const std::string& vertex_normal_name)
	{
		parameters_[&m].vertex_normal_name_ = vertex_normal_name;
	}

protected:

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			if (!p.initialized_)
				continue;
						
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_faces_)
			{
				if (p.render_edges_)
				{
					glEnable(GL_POLYGON_OFFSET_FILL);
					glPolygonOffset(1.0f, 2.0f);
				}

				if (p.phong_shading_)
				{
					p.param_phong_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES);
					p.param_phong_->release();
				}
				else
				{
					p.param_flat_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES);
					p.param_flat_->release();
				}
			
				if (p.render_edges_)
					glDisable(GL_POLYGON_OFFSET_FILL);
			}

			if (p.render_vertices_)
			{
				p.param_point_sprite_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				md->draw(rendering::POINTS);
				p.param_point_sprite_->release();
			}

			if (p.render_edges_)
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				md->draw(rendering::LINES);
				glDisable(GL_BLEND);
				p.param_edge_->release();
			}
		}
	}

    void interface() override
	{
		bool need_update = false;

		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Select mesh"))
		{
			mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
			{
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
					selected_mesh_ = m;
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];

			std::vector<Attribute<Vec3>*> vec3_attributes;
			foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (Attribute<Vec3>* attribute)
			{
				vec3_attributes.push_back(attribute);
			});

			std::string selected_vertex_position_name_ = p.vertex_position_name_.empty() ? "-- select --" : p.vertex_position_name_;
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				for (Attribute<Vec3>* attribute : vec3_attributes)
				{
					bool is_selected = attribute->name() == p.vertex_position_name_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						p.vertex_position_name_ = attribute->name();
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
			std::string selected_vertex_normal_name_ = p.vertex_normal_name_.empty() ? "-- select --" : p.vertex_normal_name_;
			if (ImGui::BeginCombo("Normal", selected_vertex_normal_name_.c_str()))
			{
				for (Attribute<Vec3>* attribute : vec3_attributes)
				{
					bool is_selected = attribute->name() == p.vertex_normal_name_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						p.vertex_normal_name_ = attribute->name();
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}

			if (!p.vertex_position_name_.empty())
			{
				if (ImGui::Button("Update data"))
					update_data(*selected_mesh_);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			need_update |= ImGui::Checkbox("Faces", &p.render_faces_);
			if (p.render_faces_)
				need_update |= ImGui::Checkbox("Phong shading", &p.phong_shading_);

			if (p.render_faces_)
			{
				if (p.phong_shading_)
				{
					ImGui::Separator();
					ImGui::Text("Phong parameters");
					need_update |= ImGui::ColorEdit3("front color##phong", p.param_phong_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
					ImGui::SameLine();
					need_update |= ImGui::ColorEdit3("back color##phong", p.param_phong_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("spec##phong", &(p.param_phong_->specular_coef_), 10.0f, 1000.0f);
					need_update |= ImGui::Checkbox("double side##phong", &(p.param_phong_->double_side_));
				}
				else
				{
					ImGui::Separator();
					ImGui::Text("Flat parameters");
					need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
					ImGui::SameLine();
					need_update |= ImGui::ColorEdit3("back color##flat", p.param_flat_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_->double_side_));
				}
			}

			if (p.render_edges_)
			{
				ImGui::Separator();
				ImGui::Text("Edges parameters");
				need_update |= ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
			}

			if (p.render_vertices_)
			{
				ImGui::Separator();
				ImGui::Text("Vertices parameters");
				need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1, 2.0);
			}
		}

		ImGui::End();

		if (need_update)
			for (ui::View* v : linked_views_)
				v->request_update();
	}

private:

	const MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_