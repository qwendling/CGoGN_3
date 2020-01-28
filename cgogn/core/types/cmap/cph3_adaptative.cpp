#include "cph3_adaptative.h"

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
	auto it = (*representative_visibility_level_)[d.index].lower_bound(l);
	if(it != (*representative_visibility_level_)[d.index].end() && *it == l)
		(*representative_visibility_level_)[d.index].erase(it);
}

bool CPH3_adaptative::representative_is_visible(Dart d)const
{
	return get_representative_visibility_level(d) <= std::max(visibility_level_,current_level_);
}

bool CPH3_adaptative::dart_is_visible(Dart d)const
{
	return get_dart_visibility_level(d) <= std::max(visibility_level_,current_level_);
}

void CPH3_adaptative::set_visibility_level(Dart d,uint32 l)
{
	(*dart_visibility_level_)[d.index].insert(l);
}

void CPH3_adaptative::unset_visibility_level(Dart d,uint32 l)
{
	auto it = (*dart_visibility_level_)[d.index].lower_bound(l);
	if(it != (*dart_visibility_level_)[d.index].end() && *it == l)
		(*dart_visibility_level_)[d.index].erase(it);
}

}

