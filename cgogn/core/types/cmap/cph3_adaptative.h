#ifndef CGOGN_CORE_TYPES_CMAP_CPH3_ADAPTATIVE_H
#define CGOGN_CORE_TYPES_CMAP_CPH3_ADAPTATIVE_H

#include <cgogn/core/types/cmap/cph3.h>
#include <set>

namespace cgogn {

struct CPH3_adaptative : public CPH3
{
	using Self = CPH3_adaptative;
	
	std::shared_ptr<Self> father_;
	std::shared_ptr<Attribute<std::multiset<uint32>>> dart_visibility_level_;
	std::shared_ptr<Attribute<std::multiset<uint32>>> representative_visibility_level_;
	std::shared_ptr<Attribute<Dart>> representative_;
	static int id;
	uint32 visibility_level_;
	
	CPH3_adaptative() = delete;
	CPH3_adaptative(CMAP& m):CPH3(m),father_(nullptr),visibility_level_(0){
		dart_visibility_level_ = m_.darts_.add_attribute<std::multiset<uint32>>("dart_visibility_level"+id);
		representative_visibility_level_ = m_.darts_.add_attribute<std::multiset<uint32>>("representative_visibility_level"+id);
		representative_ = m_.darts_.get_attribute<Dart>("representative");
		if (!representative_)
			representative_ = m_.darts_.add_attribute<Dart>("representative");
		id++;
	}
	CPH3_adaptative(const CPH3& cph):CPH3(cph),father_(nullptr),visibility_level_(0){
		dart_visibility_level_ = m_.darts_.add_attribute<std::multiset<uint32>>("dart_visibility_level"+id);
		representative_visibility_level_ = m_.darts_.add_attribute<std::multiset<uint32>>("representative_visibility_level"+id);
		representative_ = m_.darts_.get_attribute<Dart>("representative");
		if (!representative_)
			representative_ = m_.darts_.add_attribute<Dart>("representative");
		id++;
	}
	CPH3_adaptative(const CPH3_adaptative& other):CPH3(other),father_(other.father_),
		dart_visibility_level_(other.dart_visibility_level_),
		representative_visibility_level_(other.representative_visibility_level_),
		representative_(other.representative_), visibility_level_(other.visibility_level_)
	{}
	
	std::shared_ptr<CPH3_adaptative> get_child();
	
	void set_representative(Dart d,Dart r);
	Dart get_representative(Dart d)const;
	uint32 get_dart_visibility_level(Dart d)const;
	uint32 get_representative_visibility_level(Dart d)const;
	bool representative_is_visible(Dart d)const;
	bool dart_is_visible(Dart d)const;
	
	void set_representative_visibility_level(Dart d,uint32 l);
	void unset_representative_visibility_level(Dart d,uint32 l);
	void set_visibility_level(Dart d,uint32 l);
	void unset_visibility_level(Dart d,uint32 l);
};

}



#endif // CGOGN_CORE_TYPES_CMAP_CPH3_ADAPTATIVE_H
