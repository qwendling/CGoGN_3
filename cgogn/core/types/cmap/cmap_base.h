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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_H_
#define CGOGN_CORE_CMAP_CMAP_BASE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/container/vector.h>

#include <cgogn/core/types/cmap/cell.h>

#include <any>
#include <array>
#include <condition_variable>
#include <unordered_map>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMapBase
{
	// using AttributeContainer = AttributeContainerT<Vector>;
	using AttributeContainer = AttributeContainerT<ChunkArray>;

	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	using AttributeGen = AttributeContainer::AttributeGen;
	using MarkAttribute = AttributeContainer::MarkAttribute;

	/*************************************************************************/
	// Map attributes container
	/*************************************************************************/
	std::shared_ptr<std::unordered_map<std::string, std::any>> attributes_;

	/*************************************************************************/
	// Dart attributes container
	/*************************************************************************/
	mutable std::shared_ptr<AttributeContainer> darts_;

	// shortcuts to topological relations attributes
	std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>> relations_;
	// shortcuts to cells indices attributes
	std::shared_ptr<std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>> cells_indices_;

	// shortcut to boundary marker attribute
	MarkAttribute* boundary_marker_;

	/*************************************************************************/
	// Cells attributes containers
	/*************************************************************************/
	mutable std::shared_ptr<std::array<AttributeContainer, NB_ORBITS>> attribute_containers_;

	/*************************************************************************/
	// Threads management (Readers/writers)
	/*************************************************************************/

	mutable std::condition_variable cv;
	mutable std::mutex m_;
	mutable int nb_reader;
	int nb_writer;
	bool is_modify;

	CMapBase();
	CMapBase(std::shared_ptr<std::unordered_map<std::string, std::any>>& attributes,
			 std::shared_ptr<AttributeContainer>& darts,
			 std::shared_ptr<std::vector<std::shared_ptr<Attribute<Dart>>>>& relations,
			 std::shared_ptr<std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS>>& cells_indices,
			 std::shared_ptr<std::array<AttributeContainer, NB_ORBITS>>& attribute_containers);
	~CMapBase();

	template <typename T>
	T& get_attribute(const std::string& name)
	{
		auto [it, inserted] = attributes_->try_emplace(name, T());
		return std::any_cast<T&>(it->second);
	}

	CMapBase& get_copy()
	{
		CMapBase* result = new CMapBase(attributes_, darts_, relations_, cells_indices_, attribute_containers_);
		return *result;
	}

	inline std::shared_ptr<Attribute<Dart>> add_relation(const std::string& name)
	{
		return relations_->emplace_back(darts_->add_attribute<Dart>(name));
	}

	inline Dart begin() const
	{
		return Dart(darts_->first_index());
	}
	inline Dart end() const
	{
		return Dart(darts_->last_index());
	}
	inline Dart next(Dart d) const
	{
		return Dart(darts_->next_index(d.index));
	}

	void start_reader() const
	{
		std::unique_lock<std::mutex> lk(m_);
		cv.wait(lk, [&] { return nb_writer <= 0; });
		nb_reader++;
	}

	void end_reader() const
	{
		std::unique_lock<std::mutex> lk(m_);
		nb_reader--;
		cv.notify_all();
	}

	void start_writer()
	{
		std::unique_lock<std::mutex> lk(m_);
		nb_writer++;
		cv.wait(lk, [&] { return nb_reader <= 0 && !is_modify; });
		is_modify = true;
	}

	void end_writer()
	{
		std::unique_lock<std::mutex> lk(m_);
		nb_writer--;
		is_modify = false;
		cv.notify_all();
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
