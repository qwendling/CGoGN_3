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

#ifndef CGOGN_MODULE_MULTIRESOLUTION_MODELING_H_
#define CGOGN_MODULE_MULTIRESOLUTION_MODELING_H_

#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/algos/decimation/decimation.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class MultiresolutionModeling : public Module
{
	static_assert(mesh_traits<MESH>::dimension >= 3, "MultiresolutionModeling can only be used with meshes of dimension >= 3");

    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	MultiresolutionModeling(const App& app) :
		Module(app, "MultiresolutionModeling (" + std::string{mesh_traits<MESH>::name} + ")"),
		selected_mesh_(nullptr),
		selected_vertex_position_(nullptr)
	{}
	~MultiresolutionModeling()
	{}

	void subdivide_mesh(MESH& m, std::shared_ptr<Attribute<Vec3>> vertex_position)
	{
		MESH m2(m,m.cph().maximum_level());
		modeling::butterflySubdivisionVolumeAdaptative(m2,0.5f,vertex_position);
		mesh_provider_->emit_connectivity_changed(&m);
		mesh_provider_->emit_attribute_changed(&m, vertex_position.get());
	}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

    void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Mesh"))
		{
			mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
			{
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
				{
					selected_mesh_ = m;
					selected_vertex_position_.reset();
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			std::string selected_vertex_position_name_ = selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [this] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_position_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_position_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##attribute"))
					selected_vertex_position_.reset();
			}

			if (selected_vertex_position_)
			{
				if (ImGui::Button("Subdivide"))
					subdivide_mesh(*selected_mesh_, selected_vertex_position_);
			}
			auto resolution_update = ImGui::SliderInt("resolution level", &mesh_resolution, 0, selected_mesh_->cph().maximum_level());
			if(resolution_update){
				selected_mesh_->current_level(mesh_resolution);
				mesh_provider_->emit_connectivity_changed(selected_mesh_);
				mesh_provider_->emit_attribute_changed(selected_mesh_, selected_vertex_position_.get());
			}
		}
		
		ImGui::End();
	}

private:

	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	MeshProvider<MESH>* mesh_provider_;
	int32 mesh_resolution;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MULTIRESOLUTION_MODELING_H_
