#include "EMR_Map3_Adaptative.h"
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

uint32 EMR_Map3_Adaptative::nb_views = 0;

bool EMR_Map3_Adaptative::check_integrity() const
{
	for (Dart d = this->begin(), end = this->end(); d != end; d = this->next(d))
	{

		int i = 0;
		Dart it = phi1(*this, d);
		bool is_permutation = false;
		while (i < 1e6 && !is_permutation)
		{
			is_permutation = it == d;
			it = phi1(*this, it);
			i++;
		}
		if (i == 1)
		{
			std::cerr << "phi1 have a fix point : " << d.index << std::endl;
			return false;
		}
		if (!is_permutation)
		{
			std::cerr << "phi1 must be a permutation" << std::endl;
			return false;
		}
		if (phi2(*this, phi2(*this, d)) != d)
		{
			std::cerr << d.index << " phi2 must be an involution" << std::endl;
			return false;
		}
		if (phi2(*this, d) == d)
		{
			std::cerr << "phi2 have a fixpoint" << std::endl;
			return false;
		}
		if (phi3(*this, phi3(*this, d)) != d)
		{
			std::cerr << "phi3 must be an involution" << std::endl;
			return false;
		}
		if (phi3(*this, d) == d)
		{
			std::cerr << "phi3 have a fixpoint" << std::endl;
			return false;
		}
		if (phi1(*this, phi3(*this, phi1(*this, phi3(*this, d)))) != d)
		{
			std::cerr << "phi1(phi3(d)) must be an involution" << std::endl;
			return false;
		}
		if (phi1(*this, phi_1(*this, d)) != d || phi1(*this, phi_1(*this, d)) != phi_1(*this, phi1(*this, d)))
		{
			std::cerr << "phi1(phi_1(d)) must be an involution" << std::endl;
			return false;
		}
		if (is_boundary(*this, d) != is_boundary(*this, phi1(*this, d)))
		{
			std::cerr << "Face must be only boundary or not" << std::endl;
			return false;
		}
	}
	return true;
}

Dart EMR_Map3_Adaptative::get_phi1_buffer(Dart d) const
{

	if (current_level_ == maximum_level_)
		return (*((*m_.MR_phi1_)[current_level_]))[d.index];
	auto& buffer = (*phi1_buffer_)[d.index];
	if (get_parent() && clock_parent_ != get_parent()->clock_views_)
	{
		clock_views_++;
		clock_parent_ = get_parent()->clock_views_;
	}

	if (std::get<0>(buffer) != m_.clock_ || std::get<1>(buffer) != clock_views_ ||
		std::get<2>(buffer) != current_level_)
	{
		auto fn = [&]() -> Dart {
			if (is_boundary(*this, d))
			{
				Dart d3 = phi3(*this, d);
				Dart it = d3;
				Dart it2 = phi1(*this, d3);
				while (it2 != d3)
				{
					it = it2;
					it2 = phi1(*this, it2);
				}
				return phi3(*this, it);
			}
			EMR_Map3 emr = EMR_Map3(*this);
			Dart d3 = phi3(*this, d);
			uint32 l_d3 = dart_level(d3);
			uint32 l_d = dart_level(d);
			emr.current_level_ = std::max(l_d, l_d3);
			Dart test = phi1(emr, d);
			if (dart_is_visible(test))
				return test;
			if (l_d3 > l_d)
			{
				emr.current_level_ = l_d3 - 1;
				Dart tmp = phi3(emr, d);
				emr.current_level_ = l_d3;
				return phi3(emr, tmp);
			}
			emr.current_level_ = l_d3;
			Dart d_tmp = phi3(emr, d3);
			emr.current_level_--;
			d_tmp = phi3(emr, d_tmp);
			emr.current_level_++;
			return phi3(emr, d_tmp);

			if (l_d3 == l_d)
			{
				emr.current_level_ = l_d;
				return phi1(emr, d);
			}
			if (l_d3 > l_d)
			{
				emr.current_level_ = l_d3;
				Dart tmp = phi1(emr, d);
				if (get_dart_visibility(tmp) <= current_level_)
					return tmp;
				emr.current_level_--;
				tmp = phi1(emr, d);
				emr.current_level_++;
				return phi_1(emr, tmp);
			}

			// emr.current_level_ = l_d - 1;
			// Dart d_1 = phi3(emr, d3);

			////////test////////
			emr.current_level_ = l_d;
			Dart tmp = phi1(emr, d);
			if (get_dart_visibility(tmp) <= current_level_)
				return tmp;
			Dart d_1 = phi3(emr, d);
			emr.current_level_ = l_d3;
			d_1 = phi3(emr, d_1);
			////////////////////
			/*tmp = phi1(emr, d_1);
			if (get_dart_visibility(tmp) <= current_level_)
				return tmp;*/
			emr.current_level_ = l_d3 - 1;
			tmp = phi1(emr, d_1);
			emr.current_level_++;
			return phi_1(emr, tmp);
		};
		std::get<3>(buffer) = fn();
		std::get<0>(buffer) = m_.clock_;
		std::get<1>(buffer) = clock_views_;
		std::get<2>(buffer) = current_level_;
	}
	return std::get<3>(buffer);
}

Dart EMR_Map3_Adaptative::get_phi2_buffer(Dart d) const
{
	if (current_level_ == maximum_level_)
		return (*((*m_.MR_phi2_)[current_level_]))[d.index];
	auto& buffer = (*phi2_buffer_)[d.index];

	if (get_parent() && clock_parent_ != get_parent()->clock_views_)
	{
		clock_views_++;
		clock_parent_ = get_parent()->clock_views_;
	}

	if (std::get<0>(buffer) != m_.clock_ || std::get<1>(buffer) != clock_views_ ||
		std::get<2>(buffer) != current_level_)
	{
		if (current_level_ == maximum_level_)
		{
			std::get<3>(buffer) = (*((*m_.MR_phi2_)[current_level_]))[d.index];
		}
		else
		{
			EMR_Map3 emr = EMR_Map3(*this);
			Dart d3 = phi3(*this, d);
			emr.current_level_ = std::max(dart_level(d), dart_level(d3));
			Dart result = phi2(emr, d);
			while (get_dart_visibility(result) > current_level_)
			{
				result = phi2(emr, phi3(emr, result));
			}
			std::get<3>(buffer) = result;
		}
		std::get<0>(buffer) = m_.clock_;
		std::get<1>(buffer) = clock_views_;
		std::get<2>(buffer) = current_level_;
	}
	return std::get<3>(buffer);
}

Dart EMR_Map3_Adaptative::get_phi3_buffer(Dart d) const
{
	if (current_level_ == maximum_level_)
		return (*((*m_.MR_phi3_)[current_level_]))[d.index];
	auto& buffer = (*phi3_buffer_)[d.index];

	if (get_parent() && clock_parent_ != get_parent()->clock_views_)
	{
		clock_views_++;
		clock_parent_ = get_parent()->clock_views_;
	}

	if (std::get<0>(buffer) != m_.clock_ || std::get<1>(buffer) != clock_views_ ||
		std::get<2>(buffer) != current_level_)
	{
		Dart result;
		if (current_level_ == maximum_level_)
		{
			result = (*((*m_.MR_phi3_)[current_level_]))[d.index];
		}
		else
		{
			uint32 d_level = dart_level(d);
			if (d_level == maximum_level_)
			{
				result = (*((*m_.MR_phi3_)[d_level]))[d.index];
			}
			else
			{

				for (int i = maximum_level_; i >= int(d_level); --i)
				{
					result = (*((*m_.MR_phi3_)[i]))[d.index];
					if (get_dart_visibility(result) <= current_level_)
					{
						break;
					}
				}
			}
		}
		std::get<3>(buffer) = result;
		std::get<0>(buffer) = m_.clock_;
		std::get<1>(buffer) = clock_views_;
		std::get<2>(buffer) = current_level_;
	}
	return std::get<3>(buffer);
}

EMR_Map3_Adaptative* EMR_Map3_Adaptative::get_child()
{
	EMR_Map3_Adaptative* result = get_copy();
	result->parent = this;
	return result;
}
EMR_Map3_Adaptative* EMR_Map3_Adaptative::get_copy()
{
	auto result = new EMR_Map3_Adaptative(m_);
	for (Dart it = m_.begin(); it != m_.end(); it = m_.next(it))
	{
		(*result->dart_visibility_)[it.index] = (*dart_visibility_)[it.index];
	}
	result->parent = parent;
	return result;
}

void EMR_Map3_Adaptative::copy_visibility(const EMR_Map3_Adaptative& other)
{
	for (Dart it = m_.begin(); it != m_.end(); it = m_.next(it))
	{
		(*dart_visibility_)[it.index] = (*other.dart_visibility_)[it.index];
	}
	clock_views_++;
};

EMR_Map3_Adaptative* EMR_Map3_Adaptative::get_parent() const
{
	if (parent != nullptr)
		return parent;
	if (topology_ != this)
		return topology_;
	return nullptr;
}

Dart EMR_Map3_Adaptative::get_representative(Dart d) const
{
	if (get_parent())
		return get_parent()->get_representative(d);
	auto& p = (*dart_representative_)[d.index];
	if (!p.first)
	{
		uint32 d_level = this->dart_level(d);
		Dart d3 = (*((*m_.MR_phi3_)[d_level]))[d.index];
		bool same_side = false;
		while (dart_level(d3) != d_level)
		{
			d_level = dart_level(d3);
			d3 = (*((*m_.MR_phi3_)[d_level]))[d3.index];
			same_side = !same_side;
		}
		if (!same_side)
			d3 = (*((*m_.MR_phi3_)[d_level]))[d3.index];
		p.first = true;
		p.second = d3;
	}
	return p.second;
}

uint32 EMR_Map3_Adaptative::get_dart_visibility(Dart d) const
{
	uint32 d_level = this->dart_level(d);
	if (d_level <= this->current_level_)
		return d_level;

	auto& buffer = (*dart_visibility_buffer_)[d.index];

	if (get_parent() && clock_parent_ != get_parent()->clock_views_)
	{
		clock_views_++;
		clock_parent_ = get_parent()->clock_views_;
	}

	if (std::get<0>(buffer) != m_.clock_ || std::get<1>(buffer) != clock_views_ ||
		std::get<2>(buffer) != current_level_)
	{
		auto p = (*dart_visibility_)[d.index];

		uint32 result = d_level;
		if (p.first)
			result = std::min(result, p.second);
		if (!is_boundary(*this, d))
		{
			/* auto p = (*dart_visibility_)[d.index];
			if (p.first)
				result = std::min(result, p.second);*/
			if (get_parent() != nullptr)
			{
				result = std::min(get_parent()->get_dart_visibility(d), result);
			}
			Dart r = get_representative(d);
			if (r != d)
				result = std::max(get_dart_visibility(r), result);
		}
		else
		{

			if (get_parent() != nullptr)
			{
				result = std::min(get_parent()->get_dart_visibility(d), result);
			}
			/* Dart tmp = (*((*m_.MR_phi3_)[d_level]))[d.index];
			if (tmp == d)
			{
				int d_level2 = d_level + 1;
				while (tmp == d)
				{
					tmp = (*((*m_.MR_phi3_)[d_level2]))[d.index];
				}
			}
			else
			{
				if (dart_level(tmp) != d_level)
				{
					tmp = (*((*m_.MR_phi3_)[d_level - 1]))[tmp.index];
					tmp = (*((*m_.MR_phi3_)[d_level]))[tmp.index];
				}
			}
			result = std::min(get_dart_visibility(tmp),result);*/

			/* Dart tmp = (*((*m_.MR_phi3_)[d_level]))[d.index];
			if (tmp == d)
			{
				auto p = (*dart_visibility_)[d.index];
				if (p.first)
					result = std::min(result, p.second);
			}
			else
			{
				if (dart_level(tmp) != d_level)
				{
					tmp = (*((*m_.MR_phi3_)[d_level - 1]))[tmp.index];
					tmp = (*((*m_.MR_phi3_)[d_level]))[tmp.index];
				}
				if (tmp == d)
					std::cout << "problem" << std::endl;
				result = get_dart_visibility(tmp);
			}*/
		}
		std::get<3>(buffer) = result;
		std::get<0>(buffer) = m_.clock_;
		std::get<1>(buffer) = clock_views_;
		std::get<2>(buffer) = current_level_;
	}
	return std::get<3>(buffer);
}

uint32 EMR_Map3_Adaptative::get_dart_visibility_fast(Dart d) const
{
	uint32 d_level = this->dart_level(d);
	if (d_level == 0)
		return 0;

	auto p = (*dart_visibility_)[d.index];
	uint32 result = d_level;
	if (p.first)
		result = std::min(result, p.second);

	if (get_parent() != nullptr)
	{
		result = std::min(get_parent()->get_dart_visibility_fast(d), result);
	}
	return result;
	/*uint32 d_level = this->dart_level(d);
	if (d_level == 0)
		return 0;
	uint32 result = UINT_MAX;
	if (parent != nullptr)
		result = parent->get_dart_visibility_fast(d);
	auto p = (*dart_visibility_)[d.index];
	if (p.first)
		result = std::min(result, p.second);

	if (d_level < result)
		return d_level;

	return result;*/
}

void EMR_Map3_Adaptative::set_dart_visibility(Dart d, uint32 v)
{
	(*dart_visibility_)[d.index].first = true;
	(*dart_visibility_)[d.index].second = v;
}

bool EMR_Map3_Adaptative::dart_is_visible(Dart d) const
{
	return get_dart_visibility(d) <= current_level_;
}

Dart EMR_Map3_Adaptative::begin() const
{
	Dart d(m_.darts_->first_index());
	uint32 lastidx = m_.darts_->last_index();
	while (d.index < lastidx && get_dart_visibility(d) > current_level_)
		d = Dart(m_.darts_->next_index(d.index));
	return d;
}

Dart EMR_Map3_Adaptative::end() const
{
	return Dart(m_.darts_->last_index());
}

Dart EMR_Map3_Adaptative::next(Dart d) const
{
	uint32 lastidx = m_.darts_->last_index();
	do
	{
		d = Dart(m_.darts_->next_index(d.index));
	} while (d.index < lastidx && get_dart_visibility(d) > current_level_);
	return d;
}

bool EMR_Map3_Adaptative::vertex_is_visible(Dart d) const
{
	return dart_is_visible(d);
}

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

Dart EMR_Map3_Adaptative::edge_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (is_boundary(*this, d))
	{
		return edge_youngest_dart(phi3(*this, d));
	}
	Dart it = phi2(*this, d);
	if (m_.dart_level(d) >= m_.dart_level(it))
		return d;
	return it;
}

Dart EMR_Map3_Adaptative::edge_oldest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = phi2(*this, d);
	if (m_.dart_level(d) < m_.dart_level(it))
		return d;
	return it;
}

uint32 EMR_Map3_Adaptative::edge_level(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	return m_.dart_level(edge_youngest_dart(d));
}

bool EMR_Map3_Adaptative::edge_is_subdivided(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	uint32 e_level = edge_level(d);
	if (e_level == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = e_level + 1;
	if (phi1(*this, d) == phi1(m2, d))
	{
		return false;
	}
	return true;
}

bool EMR_Map3_Adaptative::is_topologycal_edge(Dart d) const
{
	uint32 e_level = edge_level(d);
	if (topology_ == nullptr)
	{
		return e_level == 0;
	}
	if (topology_->dart_is_visible(d))
		return e_level == topology_->edge_level(d);
	return false;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

Dart EMR_Map3_Adaptative::face_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");

	if (edge_level(d) == 0)
	{
		return d;
	}

	if (false && is_indexed<Face>(*this))
	{
		std::unordered_set<uint32> cell_id;
		Dart result = d;
		foreach_dart_of_orbit(*this, Face2(d), [&](Dart it) -> bool {
			auto p = cell_id.insert(index_of(static_cast<const EMR_Map3::MAP&>(*this), Face(it)));
			if (!p.second)
			{
				result = it;
				return false;
			}

			return true;
		});
		return result;
	}

	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	Dart it, it2;
	it = d;
	it2 = phi_1(*this, d);
	do
	{
		if (dart_level(it) == dart_level(it2))
			return it;
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		it2 = it;
		it = phi1(*this, it);
	} while (it != d);
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = true;
		Dart result_young = phi1(m2, old);
		result = marker.is_marked(result_young);
		if (result)
		{
			return result_young;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo face level");
	} while (!result);
	return d;
}

Dart EMR_Map3_Adaptative::face_oldest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = phi1(*this, d);
	Dart result = d;
	while (it != d)
	{
		if (m_.dart_level(it) < m_.dart_level(result))
			result = it;
		it = phi1(*this, it);
	}

	return result;
}

bool EMR_Map3_Adaptative::face_is_subdivided(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = face_level(d);
	return m2.face_is_subdivided(face_oldest_dart(d));
}

uint32 EMR_Map3_Adaptative::face_level(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (edge_level(d) == 0)
	{
		return 0;
	}
	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	foreach_dart_of_orbit(*this, Face2(d), [&](Dart it) -> bool {
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		return true;
	});
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = marker.is_marked(phi1(m2, old));
		if (result)
		{
			return m2.current_level_;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo face level");
	} while (!result);
	return m2.current_level_;
}

bool EMR_Map3_Adaptative::is_topologycal_face(Dart d) const
{
	uint32 f_level = face_level(d);
	if (topology_ == nullptr)
	{
		return f_level == 0;
	}
	if (topology_->dart_is_visible(d))
		return f_level == topology_->face_level(d);
	return false;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

Dart EMR_Map3_Adaptative::volume_youngest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");

	if (get_parent() && clock_parent_ != get_parent()->clock_views_)
	{
		clock_views_++;
		clock_parent_ = get_parent()->clock_views_;
	}

	auto& buffer = (*volume_dart_buffer_)[d.index];

	if (std::get<0>(buffer) != m_.clock_ || std::get<1>(buffer) != clock_views_ ||
		std::get<2>(buffer) != current_level_)
	{

		if (edge_level(d) == 0)
		{
			return d;
		}
		Dart old = d;
		DartMarkerStore<EMR_Map3> marker(*this);
		foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
			marker.mark(it);
			if (dart_level(it) < dart_level(old))
				old = it;
			return true;
		});
		EMR_Map3 m2(m_);
		m2.current_level_ = dart_level(old);

		auto fn = [&]() -> Dart {
			bool result = false;
			do
			{
				result = true;
				Dart result_young = phi1(m2, old);
				result = marker.is_marked(result_young);
				if (result)
				{
					return result_young;
				}
				m2.current_level_++;
				cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo volume youngest");
			} while (!result);
			return d;
		};

		Dart result = fn();

		for (Dart dd : marker.marked_darts())
		{
			auto& buffer2 = (*volume_dart_buffer_)[dd.index];
			std::get<3>(buffer2) = result;
			std::get<0>(buffer2) = m_.clock_;
			std::get<1>(buffer2) = clock_views_;
			std::get<2>(buffer2) = current_level_;
		}
	}
	return std::get<3>(buffer);
}

Dart EMR_Map3_Adaptative::volume_oldest_dart(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	Dart result = d;
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		if (m_.dart_level(it) < m_.dart_level(result))
		{
			result = it;
		}
		return true;
	});
	return result;
}
bool EMR_Map3_Adaptative::volume_is_subdivided(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = volume_level(d);
	if (dart_level(d) > m2.current_level_)
		return false;
	return m2.volume_is_subdivided(d);
}

uint32 EMR_Map3_Adaptative::volume_level(Dart d) const
{
	cgogn_message_assert(get_dart_visibility(d) <= current_level_, "Access to a dart introduced after current level");
	if (edge_level(d) == 0)
	{
		return 0;
	}

	Dart old = d;
	DartMarkerStore<EMR_Map3> marker(*this);
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		marker.mark(it);
		if (dart_level(it) < dart_level(old))
			old = it;
		return true;
	});
	EMR_Map3 m2(m_);
	m2.current_level_ = dart_level(old);

	bool result = false;
	do
	{
		result = marker.is_marked(phi1(m2, old));
		if (result)
		{
			return m2.current_level_;
		}
		m2.current_level_++;
		cgogn_message_assert(m2.current_level_ <= maximum_level_, "Pb algo volume level");
	} while (!result);
	return m2.current_level_;
}

bool EMR_Map3_Adaptative::is_topologycal_volume(Dart d) const
{
	uint32 v_level = volume_level(d);
	if (topology_ == nullptr)
	{
		return v_level == 0;
	}
	if (topology_->dart_is_visible(volume_oldest_dart(d)))
		return v_level == topology_->volume_level(volume_oldest_dart(d));
	return false;
}

/***************************************************
 *            ADAPTATIVE SUBDIVISION               *
 ***************************************************/

void EMR_Map3_Adaptative::activate_edge_subdivision(Edge e)
{
	if (!edge_is_subdivided(e.dart))
		return;

	Edge e2 = Edge(phi2(*this, e.dart));
	EMR_Map3 m2(m_);
	uint32 e_level = edge_level(e.dart);
	m2.current_level_ = e_level + 1;

	Dart d2 = phi2(m2, e.dart);
	Dart it = d2;
	do
	{
		set_dart_visibility(it, current_level_);
		it = phi<2, 3>(m2, it);
	} while (it != d2);

	d2 = phi2(m2, e2.dart);
	it = d2;
	do
	{
		set_dart_visibility(it, current_level_);
		it = phi<2, 3>(m2, it);
	} while (it != d2);
	clock_views_++;
}

void EMR_Map3_Adaptative::activate_face_subdivision(Face f)
{
	if (!face_is_subdivided(f.dart))
	{
		return;
	}

	Dart d = face_oldest_dart(f.dart);
	EMR_Map3 m2(m_);
	m2.current_level_ = face_level(f.dart);
	std::vector<Dart> vect_edges;
	Dart it = d;
	do
	{
		vect_edges.push_back(it);
		it = phi1(m2, it);
	} while (it != d);
	m2.current_level_++;

	for (auto e : vect_edges)
	{
		if (edge_level(e) < m2.current_level_)
			activate_edge_subdivision(Edge(e));
		Dart it2 = phi1(m2, e);
		while (get_dart_visibility(it2) > current_level_)
		{
			set_dart_visibility(it2, current_level_);
			set_dart_visibility(phi3(m2, it2), current_level_);
			it2 = phi1(m2, it2);
		}
	}
	clock_views_++;
}
bool EMR_Map3_Adaptative::activate_volume_subdivision(Volume v)
{

	EMR_Map3 m2(m_);
	Dart d = volume_youngest_dart(v.dart);

	uint32 v_level = dart_level(d);

	/*m2.current_level_ = v_level;
	if (dart_level(v.dart) > m2.current_level_)
		return false;
	if (!m2.volume_is_subdivided(v.dart))
	{
		return false;
	}*/

	// Dart d = volume_oldest_dart(v.dart);

	m2.current_level_ = v_level;

	if (v_level == maximum_level_)
		return false;

	std::vector<Vertex> vect_vertices;
	std::vector<Dart> vect_edge;

	foreach_incident_vertex(m2, Volume(d), [&vect_vertices](Vertex w) -> bool {
		vect_vertices.push_back(w);
		return true;
	});
	/*foreach_incident_edge(m2, Volume(d), [this, &v_level](Edge e) -> bool {
		if (edge_level(e.dart) == v_level)
			activate_edge_subdivision(e);
		return true;
	});*/

	foreach_dart_of_orbit(m2, Volume(d), [&vect_edge](Dart d2) -> bool {
		vect_edge.push_back(d2);
		return true;
	});

	m2.current_level_++;
	for (Vertex w : vect_vertices)
	{
		foreach_dart_of_orbit(m2, Volume(w.dart), [this, &m2](Dart d) -> bool {
			set_dart_visibility(d, current_level_);
			set_dart_visibility(phi3(m2, d), current_level_);
			return true;
		});
	}
	for (Dart d : vect_edge)
	{
		Dart d2 = phi2(m2, d);
		Dart it = d2;
		do
		{
			set_dart_visibility(it, current_level_);
			it = phi<2, 3>(m2, it);
		} while (it != d2);
	}
	clock_views_++;
	return true;
}

// In this function we assume that v.dart is the youngest dart allow more pre computation
bool EMR_Map3_Adaptative::activate_volume_subdivision_fast(Volume v)
{

	EMR_Map3 m2(m_);
	Dart d = v.dart;

	uint32 v_level = dart_level(d);

	m2.current_level_ = v_level;

	std::vector<Vertex> vect_vertices;
	std::vector<Dart> vect_edge;

	foreach_incident_vertex(m2, Volume(d), [&vect_vertices](Vertex w) -> bool {
		vect_vertices.push_back(w);
		return true;
	});

	foreach_dart_of_orbit(m2, Volume(d), [&vect_edge](Dart d2) -> bool {
		vect_edge.push_back(d2);
		return true;
	});

	m2.current_level_++;
	for (Vertex w : vect_vertices)
	{
		foreach_dart_of_orbit(m2, Volume(w.dart), [this, &m2](Dart d) -> bool {
			set_dart_visibility(d, current_level_);
			set_dart_visibility(phi3(m2, d), current_level_);
			return true;
		});
	}
	for (Dart d : vect_edge)
	{
		Dart d2 = phi2(m2, d);
		Dart it = d2;
		do
		{
			set_dart_visibility(it, current_level_);
			it = phi<2, 3>(m2, it);
		} while (it != d2);
	}
	clock_views_++;
	return true;
}

bool EMR_Map3_Adaptative::disable_edge_subdivision(Edge e)
{
	uint32 e_level = edge_level(e.dart);
	if (e_level <= current_level_)
		return false;
	if (is_topologycal_edge(e.dart))
	{
		return false;
	}
	Dart old = edge_oldest_dart(e.dart);
	if (dart_level(old) == e_level)
		return false;
	EMR_Map3 m2(m_);
	m2.current_level_ = e_level - 1;

	// equilibrage des deux parties de l'arete
	Dart d2 = phi3(m2, old);
	while (edge_level(d2) != e_level)
		if (!disable_edge_subdivision(Edge(d2)))
			return false;

	// test des faces adjacentes
	Dart test = old;
	do
	{
		if (phi3(*this, phi1(*this, test)) != phi3(m2, test))
			return false;
		test = phi2(*this, phi3(*this, test));
	} while (test != old);

	// deactivation de l'arete
	m2.current_level_ = e_level;
	Dart it, it2;
	it = old;
	it2 = d2;
	do
	{
		Dart tmp = phi2(m2, it);
		Dart tmp2 = phi2(m2, it2);
		set_dart_visibility(tmp, UINT_MAX);
		set_dart_visibility(tmp2, UINT_MAX);
		it = phi3(m2, tmp);
		it2 = phi3(m2, tmp2);
	} while (it != old);
	clock_views_++;
	return true;
}
bool EMR_Map3_Adaptative::disable_face_subdivision(Face f, bool disable_edge, bool disable_subface)
{
	uint32 f_level = face_level(f.dart);
	if (f_level == 0)
		return false;
	if (is_topologycal_face(f.dart))
	{
		return false;
	}
	EMR_Map3 m2(m_);
	m2.current_level_ = f_level;
	Dart old = face_oldest_dart(f.dart);
	Dart test = phi1(m2, old);
	// Check that the two adjacents volumes are not subdivide
	if (phi<2, 3, 2, 3>(*this, test) != test)
		return false;

	std::vector<Dart> vec_vertices;
	m2.current_level_ = f_level - 1;
	Dart it = old;
	// Simplification of the subfaces
	do
	{
		vec_vertices.push_back(it);
		if (face_level(it) != f_level)
		{
			if (!disable_subface)
			{
				return false;
			}
			if (!disable_face_subdivision(Face(it), disable_edge, disable_subface))
				return false;
		}
		it = phi1(m2, it);
	} while (it != old);

	m2.current_level_ = f_level;

	std::vector<Dart> list_dart_disable;
	for (Dart d : vec_vertices)
	{
		Dart it = phi1(m2, d);
		Dart d11 = phi<1, 1>(m2, d);
		while (it != d11)
		{
			Dart d3 = phi3(*this, it);
			Dart it2 = d3;
			m2.current_level_ = std::max(dart_level(it), dart_level(d3));
			do
			{
				list_dart_disable.push_back(it2);
				it2 = phi2(m2, it2);
				list_dart_disable.push_back(it2);
				it2 = phi3(m2, it2);
			} while (it2 != d3);
			m2.current_level_ = f_level;
			it = phi1(*this, it);
		}
	}
	for (Dart d : list_dart_disable)
	{
		set_dart_visibility(d, UINT_MAX);
	}
	clock_views_++;
	if (disable_edge)
	{
		for (Dart d : vec_vertices)
		{
			disable_edge_subdivision(Edge(d));
		}
	}
	return true;
}
bool EMR_Map3_Adaptative::disable_volume_subdivision(Volume v, bool disable_face)
{
	uint32 v_level = volume_level(v.dart);
	if (v_level <= current_level_)
		return false;
	if (is_topologycal_volume(v.dart))
	{
		return false;
	}

	EMR_Map3 m2(m_);
	m2.current_level_ = v_level;
	Dart old;
	foreach_dart_of_orbit(m2, v, [&](Dart d) -> bool {
		if (dart_level(d) <= v_level - 1)
		{
			old = d;
			return false;
		}
		return true;
	});
	DartMarkerStore<EMR_Map3_Adaptative> dm(*this);
	CellMarkerStore<EMR_Map3_Adaptative, Vertex> vm(*this);
	std::vector<Dart> vect_vertices;
	std::vector<Dart> vect_volume;

	m2.current_level_ = v_level - 1;
	foreach_dart_of_orbit(m2, Volume(old), [&](Dart d) -> bool {
		dm.mark(d);
		vect_volume.push_back(d);
		if (!vm.is_marked(Vertex(d)))
		{
			vm.mark(Vertex(d));
			vect_vertices.push_back(d);
		}
		return true;
	});
	for (Dart d : vect_vertices)
	{
		while (volume_level(d) != v_level)
			disable_volume_subdivision(Volume(d), disable_face);
	}

	m2.current_level_ = v_level;
	std::vector<Dart> vect_dart;
	for (Dart d : vect_volume)
	{
		Dart tmp = phi<1, 2>(m2, d);
		if (dm.is_marked(tmp))
			continue;
		while (face_level(tmp) != v_level)
			disable_face_subdivision(Face(tmp), true, true);
		Dart it = tmp;
		do
		{
			Dart it2 = phi3(m2, it);
			vect_dart.push_back(it);
			dm.mark(it);
			vect_dart.push_back(it2);
			dm.mark(it2);
			it = phi1(m2, it);
		} while (it != tmp);
	}
	for (Dart d : vect_dart)
	{
		set_dart_visibility(d, UINT_MAX);
	}
	clock_views_++;
	if (disable_face)
	{
		for (Dart d : vect_volume)
		{
			while (face_level(d) != v_level - 1)
			{
				if (!disable_face_subdivision(Face(d), true, true))
					break;
			}
		}
	}

	return true;
}

bool EMR_Map3_Adaptative::disable_volume_subdivision_fast(Volume v, bool disable_face)
{
	uint32 v_level = dart_level(v.dart);
	if (v_level <= current_level_)
		return false;

	EMR_Map3 m2(m_);
	m2.current_level_ = v_level;
	Dart old;
	foreach_dart_of_orbit(m2, v, [&](Dart d) -> bool {
		if (dart_level(d) <= v_level - 1)
		{
			old = d;
			return false;
		}
		return true;
	});

	static DartMarker<EMR_Map3> dm(*this);
	static CellMarker<EMR_Map3, Vertex> vm(*this);

	std::vector<Dart> vect_vertices;
	std::vector<Dart> vect_volume;

	m2.current_level_ = v_level - 1;
	foreach_dart_of_orbit(m2, Volume(old), [&](Dart d) -> bool {
		dm.mark(d);
		vect_volume.push_back(d);
		if (!vm.is_marked(Vertex(d)))
		{
			vm.mark(Vertex(d));
			vect_vertices.push_back(d);
		}
		return true;
	});

	for (Dart d : vect_vertices)
	{
		while (volume_level(d) != v_level)
			disable_volume_subdivision(Volume(d), disable_face);
	}

	m2.current_level_ = v_level;
	std::vector<Dart> vect_dart;
	for (Dart d : vect_volume)
	{
		Dart tmp = phi<1, 2>(m2, d);
		if (dm.is_marked(tmp))
			continue;
		while (face_level(tmp) != v_level)
			disable_face_subdivision(Face(tmp), true, true);
		Dart it = tmp;
		do
		{
			Dart it2 = phi3(m2, it);
			vect_dart.push_back(it);
			dm.mark(it);
			vect_dart.push_back(it2);
			dm.mark(it2);
			it = phi1(m2, it);
		} while (it != tmp);
	}
	for (Dart d : vect_dart)
	{
		set_dart_visibility(d, UINT_MAX);
	}
	clock_views_++;
	if (disable_face)
	{
		for (Dart d : vect_volume)
		{
			while (face_level(d) != v_level - 1)
			{
				if (!disable_face_subdivision(Face(d), true, true))
					break;
			}
		}
	}

	for (Dart d : vect_volume)
	{
		dm.unmark(d);
	}
	for (Dart d : vect_vertices)
	{
		vm.unmark(Vertex(d));
	}

	return true;
}

} // namespace cgogn
