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

#include <cgogn/core/types/cmap/cmap_base.h>

namespace cgogn
{

CMapBase::CMapBase() : nb_reader(0), nb_writer(0), is_modify(false)
{
	darts_ = std::shared_ptr<AttributeContainer>(new AttributeContainer());
	attributes_ =
		std::shared_ptr<std::unordered_map<std::string, std::any>>(new std::unordered_map<std::string, std::any>());
	boundary_marker_ = darts_->get_mark_attribute();
	relations_ = std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>(
		new std::vector<std::shared_ptr<Attribute<Dart>>>());
	cells_indices_ = std::shared_ptr<std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>>(
		new std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>());
	attribute_containers_ =
		std::shared_ptr<std::array<AttributeContainer, NB_ORBITS>>(new std::array<AttributeContainer, NB_ORBITS>());
}

CMapBase::CMapBase(std::shared_ptr<std::unordered_map<std::string, std::any>>& attributes,
				   std::shared_ptr<AttributeContainer>& darts,
				   std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>& relations,
				   std::shared_ptr<std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>>& cells_indices,
				   MarkAttribute* boundary_marker,
				   std::shared_ptr<std::array<AttributeContainer, NB_ORBITS>>& attribute_containers)
	: attributes_(attributes), darts_(darts), relations_(relations), cells_indices_(cells_indices),
	  boundary_marker_(boundary_marker), attribute_containers_(attribute_containers), nb_reader(0), nb_writer(0),
	  is_modify(false)
{
}

CMapBase::~CMapBase()
{
}

} // namespace cgogn
