#include "cph3_adaptative.h"
#include <cgogn/core/types/cmap/phi.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/edge.h>

namespace cgogn {

int CPH3_adaptative::id = 0;

std::shared_ptr<CPH3_adaptative> CPH3_adaptative::get_child()
{
	std::shared_ptr<CPH3_adaptative> result = std::make_shared<CPH3_adaptative>(CPH3_adaptative(CPH3(*this)));
	result->father_ = std::shared_ptr<CPH3_adaptative>(this);
	return result;
}

uint32 CPH3_adaptative::get_dart_visibility_level(Dart d)const
{
	uint32 result = dart_level(d);
	if(result == 0)
		return result;
	if(representative_is_visible(d)){
		auto it = (*dart_visibility_level_)[d.index].begin();
		if(it != (*dart_visibility_level_)[d.index].end() && *it < result)
			result = *it;
	}
	if(father_ != nullptr){
		return std::min(father_->get_dart_visibility_level(d),result);
	}
	return result;
}

void CPH3_adaptative::set_representative(Dart d,Dart r){
	(*representative_)[d.index] = r;
}

Dart CPH3_adaptative::get_representative(Dart d)const{
	if(dart_level(d) == 0){
		return d;
	}
	return (*representative_)[d.index];
}

uint32 CPH3_adaptative::get_representative_visibility_level(Dart d)const
{
	if(dart_level(d) == 0)
		return 0;
	Dart r = get_representative(d);
	auto it = (*representative_visibility_level_)[r.index].begin();
	uint32 result = dart_level(r);
	if(it != (*representative_visibility_level_)[r.index].end() && *it < result)
		result = *it;
	if(father_ != nullptr){
		return std::min(father_->get_representative_visibility_level(d),result);
	}
	return result;
}

void CPH3_adaptative::set_representative_visibility_level(Dart d,uint32 l){
	(*representative_visibility_level_)[d.index].insert(l);
}
void CPH3_adaptative::unset_representative_visibility_level(Dart d,uint32 l){
	(*representative_visibility_level_)[d.index].erase(l);
}

bool CPH3_adaptative::representative_is_visible(Dart d)const
{
	return get_representative_visibility_level(d) <= current_level_;
}

bool CPH3_adaptative::dart_is_visible(Dart d)const
{
	return get_dart_visibility_level(d) <= current_level_;
}

void CPH3_adaptative::set_visibility_level(Dart d,uint32 l)
{
	(*dart_visibility_level_)[d.index].insert(l);
}

void CPH3_adaptative::unset_visibility_level(Dart d,uint32 l)
{
	(*dart_visibility_level_)[d.index].erase(l);
}

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

uint32 CPH3_adaptative::edge_level(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	return std::max(dart_level(d), dart_level(phi2(*this, d)));
}

Dart CPH3_adaptative::edge_oldest_dart(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	Dart d2 = phi2(*this, d);
	if (dart_level(d) < dart_level(d2))
		return d;
	return d2;
}

Dart CPH3_adaptative::edge_youngest_dart(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	Dart d2 = phi2(*this, d);
	if (dart_level(d) > dart_level(d2))
		return d;
	return d2;
}

bool CPH3_adaptative::edge_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	if (current_level_ == maximum_level_)
		return false;

	Dart d1 = phi1(*this, d);
	CPH3_adaptative m2(*this);
	m2.current_level_ = edge_level(d) + 1;
	Dart d1_l = phi1(m2, d);
	if (d1 != d1_l)
		return true;
	else
		return false;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

uint32 CPH3_adaptative::face_level(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	Dart it = d;
	Dart old = it;
	uint32 l_old = dart_level(old);
	uint32 fLevel = edge_level(it);
	do
	{
		it = phi1(*this, it);
		uint32 dl = dart_level(it);
		uint32 l = edge_level(it);
		
		if(l<fLevel){
			fLevel = l;
			old = it;
			l_old = dl;
		}else{
			if(l == fLevel && dl < l_old){
				old = it;
				l_old = dl;
			}
		}

		// compute the oldest dart of the face in the same time
		if (dl < l_old)
		{
			old = it;
			l_old = dl;
		}
		
		fLevel = l < fLevel ? l : fLevel;
	} while (it != d);

	uint32 nbSubd = 0;
	it = old;
	uint32 eId = edge_id(old);
	uint32 init_dart_level = dart_level(it);
	do
	{
		++nbSubd;
		it = phi1(*this, it);
	} while (edge_id(it) == eId && dart_level(it) != init_dart_level);

	while (nbSubd > 1)
	{
		nbSubd /= 2;
		--fLevel;
	}

	return fLevel;
	
	/*Dart it = d;
	Dart old = it;
	uint32 l_old = dart_level(old);
	uint32 fLevel = edge_level(it);
	do
	{
		it = phi1(*this, it);
		uint32 dl = dart_level(it);

		// compute the oldest dart of the face in the same time
		if (dl < l_old)
		{
			old = it;
			l_old = dl;
		}
		uint32 l = edge_level(it);
		fLevel = l < fLevel ? l : fLevel;
	} while (it != d);

	CPH3 m2(CPH3(*this));
	m2.current_level_ = fLevel;

	uint32 nbSubd = 0;
	it = old;
	uint32 eId = m2.edge_id(old);
	uint32 init_dart_level = m2.dart_level(it);
	do
	{
		++nbSubd;
		it = phi1(m2, it);
	} while (m2.edge_id(it) == eId && m2.dart_level(it) != init_dart_level);

	while (nbSubd > 1)
	{
		nbSubd /= 2;
		--fLevel;
	}

	return fLevel;*/
}

Dart CPH3_adaptative::face_oldest_dart(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	Dart it = d;
	Dart oldest = it;
	uint32 l_old = dart_level(oldest);
	do
	{
		uint32 l = dart_level(it);
		if (l == 0)
			return it;

		if (l < l_old)
		//		if(l < l_old || (l == l_old && it < oldest))
		{
			oldest = it;
			l_old = l;
		}
		it = phi1(*this, it);
	} while (it != d);

	return oldest;
}

Dart CPH3_adaptative::face_youngest_dart(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	Dart it = d;
	Dart youngest = it;
	uint32 l_young = UINT32_MAX;
	Dart it2 = phi_1(*this,it);
	do{
		if(edge_id(it) != edge_id(it2)){
			uint32 l = dart_level(it);
			if (l > l_young)
			{
				youngest = it;
				l_young = l;
			}else{
				if(l == l_young && it.index <  youngest.index)
					youngest = it;
			}
		}
		it2=it;
		it = phi1(*this,it);
	}while(it != d);
	
	it = phi<31>(*this,youngest);
	if(it.index < youngest.index)
		youngest = it;
	
	/*do
	{
		uint32 l = dart_level(it);
		if (l > l_young)
		{
			youngest = it;
			l_young = l;
		}
		it = phi1(*this, it);
	} while (it != d);*/

	return youngest;
}

Dart CPH3_adaptative::face_origin(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	Dart p = d;
	uint32 pLevel = dart_level(p);
	CPH3 m2(CPH3(*this));
	do
	{
		m2.current_level_ = pLevel;
		p = m2.face_oldest_dart(p);
		pLevel = dart_level(p);
	} while (pLevel > 0);
	return p;
}

bool CPH3_adaptative::face_is_subdivided(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;
	uint32 fLevel = face_level(d) ;
	if(fLevel < current_level_)
		return false ;

	do{
		d = phi1(*this,d);
	}while(dart_level(d) > fLevel);

	bool subd = false ;

	CPH3 m2(CPH3(*this));
	m2.current_level_ = fLevel+1;
	if(dart_level(phi1(m2,d)) == m2.current_level_
		&& edge_id(phi1(m2,d)) != edge_id(d))
		subd = true;

	return subd;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

uint32 CPH3_adaptative::volume_level(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;

	// The level of a volume is the
	// minimum of the levels of its faces

	Dart oldest = d;
	uint32 lold = dart_level(oldest);
	uint32 vLevel = std::numeric_limits<uint32>::max();

	foreach_incident_face(*this, CPH3_adaptative::CMAP::Volume(d), [&](CPH3_adaptative::CMAP::Face f) -> bool {
		uint32 fLevel = face_level(f.dart);
		vLevel = fLevel < vLevel ? fLevel : vLevel;
		Dart old = face_oldest_dart(f.dart);
		if (dart_level(old) < lold)
		{
			oldest = old;
			lold = dart_level(old);
		}
		return true;
	});

	CPH3 m2(CPH3(*this));
	m2.current_level_ = vLevel;

	uint32 nbSubd = 0;
	Dart it = oldest;
	uint32 eId = edge_id(oldest);
	do
	{
		++nbSubd;
		it = phi<121>(m2, it);
	} while (edge_id(it) == eId && lold != dart_level(it));

	while (nbSubd > 1)
	{
		nbSubd /= 2;
		--vLevel;
	}


	return vLevel;
}

Dart CPH3_adaptative::volume_oldest_dart(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;

	Dart oldest = d;
	uint32 l_old = dart_level(oldest);
	foreach_incident_face(*this, CPH3_adaptative::CMAP::Volume(oldest), [&](CPH3_adaptative::CMAP::Face f) -> bool {
		Dart old = face_oldest_dart(f.dart);
		uint32 l = dart_level(old);
		if (l < l_old)
		{
			oldest = old;
			l_old = l;
		}
		return true;
	});

	return oldest;
}

Dart CPH3_adaptative::volume_youngest_dart(Dart d)const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;

	Dart youngest = d;
	uint32 l_young = dart_level(youngest);
	foreach_incident_face(*this, CPH3_adaptative::CMAP::Volume(youngest), [&](CPH3_adaptative::CMAP::Face f) -> bool {
		Dart young = face_youngest_dart(f.dart);
		uint32 l = dart_level(young);
		if (l > l_young)
		{
			youngest = young;
			l_young = l;
		}
		return true;
	});

	return youngest;
}

bool CPH3_adaptative::volume_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_is_visible(d),
						 "Access to a dart not visible at this level") ;

	uint vLevel = volume_level(d);
	if (vLevel < current_level_)
		return false;

	CPH3 m2(CPH3(*this));
	d = volume_oldest_dart(d);
	m2.current_level_ = vLevel;
	return m2.volume_is_subdivided(d);
}

/***************************************************
 *            ADAPTATIVE RESOLUTION                *
 ***************************************************/

void CPH3_adaptative::activated_edge_subdivision(CMAP::Edge e){
	if(!edge_is_subdivided(e.dart))
		return;
	CMAP::Edge e2 = CMAP::Edge(phi2(*this,e.dart));
	CPH3 m2(CPH3(*this));
	m2.current_level_ = edge_level(e.dart)+1;
	foreach_dart_of_orbit(m2,e,[&](Dart d)->bool{
		set_visibility_level(d,current_level_);
		return true;
	});
	foreach_dart_of_orbit(m2,e2,[&](Dart d)->bool{
		set_visibility_level(d,current_level_);
		return true;
	});
}

void CPH3_adaptative::activated_face_subdivision(CMAP::Face f){
	if(!face_is_subdivided(f.dart))
		return;
	Dart d = face_oldest_dart(f.dart);
	CPH3 m2(CPH3(*this));
	m2.current_level_ = face_level(f.dart);
	std::vector<Dart> vect_edges;
	Dart it = d;
	do{
		vect_edges.push_back(it);
		it = phi1(m2,it);
	}while(it != d);
	m2.current_level_++;
	
	for(auto e : vect_edges){
		if(edge_level(e) < m2.current_level_)
			activated_edge_subdivision(CMAP::Edge(e));
		foreach_dart_of_orbit(m2,CMAP::Edge(phi1(m2,e)),[&](Dart dd)->bool{
			set_visibility_level(dd,current_level_);
			return true;
		});
		Dart it = phi1(m2,e);
		set_representative_visibility_level(it,current_level_);
		set_representative_visibility_level(phi3(m2,it),current_level_);
		it = phi1(m2,it);
		set_representative_visibility_level(it,current_level_);
		set_representative_visibility_level(phi3(m2,it),current_level_);
	}
}

}

