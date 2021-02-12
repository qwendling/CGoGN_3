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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP0_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP0_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cell.h>
#include <cgogn/core/types/cmap/cmap_base.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap0 : public CMapBase
{
	static const uint8 dimension = 0;

	using Vertex = Cell<DART>;

	using Cells = std::tuple<Vertex>;

	CMap0()
	{
	}
	CMap0(std::shared_ptr<std::unordered_map<std::string, std::any>>& attributes,
		  std::shared_ptr<AttributeContainer>& darts,
		  std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>& relations,
		  std::shared_ptr<std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>>& cells_indices,
		  std::shared_ptr<std::array<AttributeContainer, NB_ORBITS>>& attribute_containers)
		: CMapBase(attributes, darts, relations, cells_indices, attribute_containers){};
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP0_H_
