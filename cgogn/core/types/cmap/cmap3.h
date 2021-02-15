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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP3_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP3_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap2.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap3 : public CMap2
{
	static const uint8 dimension = 3;

	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Cell<PHI2>;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;

	using Cells = std::tuple<Vertex, Vertex2, HalfEdge, Edge, Edge2, Face, Face2, Volume>;

	std::shared_ptr<Attribute<Dart>> phi3_;

	CMap3() : CMap2()
	{
		phi3_ = add_relation("phi3");
	}
	CMap3(std::shared_ptr<Attribute<Dart>>& phi3) : CMap2(), phi3_(phi3)
	{
	}
	CMap3(std::shared_ptr<std::unordered_map<std::string, std::any>>& attributes,
		  std::shared_ptr<AttributeContainer>& darts,
		  std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>& relations,
		  std::shared_ptr<std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>>& cells_indices,
		  MarkAttribute* boundary_marker,
		  std::shared_ptr<std::array<AttributeContainer, NB_ORBITS>>& attribute_containers,
		  std::shared_ptr<Attribute<Dart>>& phi1, std::shared_ptr<Attribute<Dart>>& phi_1,
		  std::shared_ptr<Attribute<Dart>>& phi2, std::shared_ptr<Attribute<Dart>>& phi3)
		: CMap2(attributes, darts, relations, cells_indices, boundary_marker, attribute_containers, phi1, phi_1, phi2),
		  phi3_(phi3){};

	CMap3* get_copy()
	{
		CMap3* result = new CMap3(attributes_, darts_, relations_, cells_indices_, boundary_marker_,
								  attribute_containers_, phi1_, phi_1_, phi2_, phi3_);
		return result;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP3_H_
