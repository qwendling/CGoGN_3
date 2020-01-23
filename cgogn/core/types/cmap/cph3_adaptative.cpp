#include "cph3_adaptative.h"

namespace cgogn {

std::shared_ptr<CPH3_adaptative> CPH3_adaptative::get_child()
{
	std::shared_ptr<CPH3_adaptative> result = std::make_shared<CPH3_adaptative>(CPH3_adaptative(this->cph_));
	result->father_ = std::shared_ptr<CPH3_adaptative>(this);
	return result;
}

uint32 CPH3_adaptative::get_dart_visibility_level(Dart d)const
{
	auto it = (*dart_visibility_level_)[d.index].begin();
	uint32 result;
	if(it != (*dart_visibility_level_)[d.index].end())
		result = std::min(cph_.dart_level(d),*it);
	else
		result = cph_.dart_level(d);
	if(father_ != nullptr){
		return std::min(father_->get_dart_visibility_level(d),result);
	}
	return result;
}

uint32 CPH3_adaptative::get_representative_visibility_level(Dart d)const
{
	Dart r = (*representative_)[d.index];
	auto it = (*representative_visibility_level_)[r.index].begin();
	uint32 result;
	if(it != (*representative_visibility_level_)[r.index].end())
		result = std::min(cph_.dart_level(r),*it);
	else
		result = cph_.dart_level(r);
	if(father_ != nullptr){
		return std::min(father_->get_representative_visibility_level(d),result);
	}
	return result;
}

bool CPH3_adaptative::representative_is_visible(Dart d)const
{
	return get_representative_visibility_level(d) >= visibility_level_;
}

bool CPH3_adaptative::dart_is_visible(Dart d)const
{
	return get_dart_visibility_level(d) >= visibility_level_ && representative_is_visible(d);
}

void CPH3_adaptative::set_representative_visibility_level(Dart d,uint32 l)
{
	(*representative_visibility_level_)[(*representative_)[d.index].index].insert(l);
}
void CPH3_adaptative::unset_representative_visibility_level(Dart d,uint32 l)
{
	Dart r = (*representative_)[d.index];
	auto it = (*representative_visibility_level_)[r.index].lower_bound(l);
	if(it != (*representative_visibility_level_)[r.index].end() && *it == l)
		(*representative_visibility_level_)[r.index].erase(it);
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

