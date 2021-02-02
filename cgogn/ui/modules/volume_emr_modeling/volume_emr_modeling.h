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

#ifndef CGOGN_MODULE_VOLUME_EMR_MODELING_H_
#define CGOGN_MODULE_VOLUME_EMR_MODELING_H_

#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/decimation/decimation.h>
#include <cgogn/modeling/algos/subdivision.h>

#include <imgui/imgui.h>

namespace cgogn
{

namespace ui
{

class VolumeEMRModeling : public Module
{
	using MRMesh = EMR_Map3;

	template <typename T>
	using Attribute = typename mesh_traits<MRMesh>::template Attribute<T>;

	using Vertex = typename mesh_traits<MRMesh>::Vertex;
	using Edge = typename mesh_traits<MRMesh>::Edge;

	using Vec3 = geometry::Vec3;

public:
	VolumeEMRModeling(const App& app)
		: Module(app, "VolumeEMRModeling (" + std::string{mesh_traits<MRMesh>::name} + ")"),
		  selected_vertex_relative_position_(nullptr), selected_vertex_parents_(nullptr), selected_cph3_(nullptr),
		  selected_cmap3_(nullptr), selected_vertex_position_(nullptr), selected_vertex_attr2_(nullptr),
		  selected_vertex_attr3_(nullptr)
	{
	}
	~VolumeEMRModeling()
	{
	}

	MRMesh* create_mrmesh(MRMesh::MAP& m, const std::string& name)
	{
		MRMesh* result = new MRMesh(m);
		std::string emr_name;
		uint32 count = 1;
		do
		{
			emr_name = name + "_" + std::to_string(count);
			++count;
		} while (emr_provider_->has_mesh(emr_name));
		emr_provider_->register_mesh(result, emr_name);
		return result;
	}

	void changed_connectivity(MRMesh& m, Attribute<Vec3>* vertex_position)
	{

		emr_provider_->emit_connectivity_changed(&m);
		emr_provider_->emit_attribute_changed(&m, vertex_position);

		cmap3_provider_->emit_connectivity_changed(&static_cast<MRMesh::MAP&>(m));
		cmap3_provider_->emit_attribute_changed(&static_cast<MRMesh::MAP&>(m), vertex_position);
	}

protected:
	void init() override
	{
		emr_provider_ = static_cast<ui::MeshProvider<MRMesh>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MRMesh>::name} + ")"));

		cmap3_provider_ = static_cast<ui::MeshProvider<MRMesh::MAP>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MRMesh::MAP>::name} + ")"));
	}

	void interface() override
	{

		if (ImGui::ListBoxHeader("CMap3"))
		{
			cmap3_provider_->foreach_mesh([this](MRMesh::MAP* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_cmap3_))
				{
					selected_cmap3_ = m;
					selected_cmap3_name_ = name;
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_cmap3_)
		{
			if (ImGui::Button("Create MRMesh"))
				create_mrmesh(*selected_cmap3_, selected_cmap3_name_);
		}

		if (ImGui::ListBoxHeader("MRMesh"))
		{
			emr_provider_->foreach_mesh([this](MRMesh* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_cph3_))
				{
					selected_cph3_ = m;
					selected_vertex_position_.reset();
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_cph3_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			uint32 min = 0;

			std::string selected_vertex_position_name_ =
				selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
			std::string selected_vertex_attr2_name_ =
				selected_vertex_attr2_ ? selected_vertex_attr2_->name() : "-- select --";
			std::string selected_vertex_attr3_name_ =
				selected_vertex_attr3_ ? selected_vertex_attr3_->name() : "-- select --";
			std::string selected_vertex_parents_name_ =
				selected_vertex_parents_ ? selected_vertex_parents_->name() : "-- select --";
			std::string selected_vertex_relative_name_ =
				selected_vertex_relative_position_ ? selected_vertex_relative_position_->name() : "-- select --";
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_cph3_,
												[this](const std::shared_ptr<Attribute<Vec3>>& attribute) {
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
			if (ImGui::BeginCombo("Attr2", selected_vertex_attr2_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_cph3_,
												[this](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == selected_vertex_attr2_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														selected_vertex_attr2_ = attribute;
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (selected_vertex_attr2_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##attr2"))
					selected_vertex_attr2_.reset();
			}
			if (ImGui::BeginCombo("Attr3", selected_vertex_attr3_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_cph3_,
												[this](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == selected_vertex_attr3_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														selected_vertex_attr3_ = attribute;
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (selected_vertex_attr3_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##attr3"))
					selected_vertex_attr3_.reset();
			}

			if (ImGui::BeginCombo("Parents", selected_vertex_parents_name_.c_str()))
			{
				foreach_attribute<std::array<Vertex, 3>, Vertex>(
					*selected_cph3_, [this](const std::shared_ptr<Attribute<std::array<Vertex, 3>>>& attribute) {
						bool is_selected = attribute == selected_vertex_parents_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
							selected_vertex_parents_ = attribute;
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
				ImGui::EndCombo();
			}
			if (selected_vertex_parents_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##parents"))
					selected_vertex_parents_.reset();
			}

			if (ImGui::BeginCombo("Relative position", selected_vertex_relative_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_cph3_,
												[this](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == selected_vertex_relative_position_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														selected_vertex_relative_position_ = attribute;
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (selected_vertex_relative_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##relative"))
					selected_vertex_relative_position_.reset();
			}
		}
	}

public:
	std::shared_ptr<Attribute<Vec3>> selected_vertex_relative_position_;
	std::shared_ptr<Attribute<std::array<Vertex, 3>>> selected_vertex_parents_;

private:
	MRMesh* selected_cph3_;
	MRMesh::MAP* selected_cmap3_;
	std::string selected_cmap3_name_;

	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_attr2_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_attr3_;

	MeshProvider<MRMesh>* emr_provider_;
	MeshProvider<MRMesh::MAP>* cmap3_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_VOLUME_EMR_MODELING_H_
