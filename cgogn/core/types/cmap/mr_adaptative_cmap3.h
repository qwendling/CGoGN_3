#ifndef CGOGN_CORE_TYPES_CMAP_MR_ADAPTATIVE_CMAP3_H
#define CGOGN_CORE_TYPES_CMAP_MR_ADAPTATIVE_CMAP3_H

#include <cgogn/core/types/cmap/mr_cmap3.h>

namespace cgogn{
/*
class MRadaptativeCMap3 : public MRCmap3{
public:
	using Self = MRadaptativeCMap3;
	using Inherit = MRCmap3;

	using Vertex = MRCmap3::Vertex;
	using Vertex2 = MRCmap3::Vertex2;
	using Edge = MRCmap3::Edge;
	using Edge2 = MRCmap3::Edge2;
	using Face = MRCmap3::Face;
	using Face2 = MRCmap3::Face2;
	using Volume = MRCmap3::Volume;
	using CC = MRCmap3::CC;

	using Cells = std::tuple<Vertex, Vertex2, Edge, Edge2, Face, Face2, Volume>;

	using AttributeContainer = AttributeContainerT<ChunkArray>;
	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	uint32 current_level_;
	
	static const bool is_mesh_view = true;
	
	std::shared_ptr<CPH3<CMap3>> base_cph_;

	MRadaptativeCMap3():current_level_(0u){
		base_cph_ = std::make_shared<CPH3<CMap3>>();
	}
	MRadaptativeCMap3(std::shared_ptr<CMap3> m):current_level_(0u){
		base_cph_ = std::make_shared<CPH3<CMap3>>(m);
	}
	MRadaptativeCMap3(const std::shared_ptr<CPH3<CMap3>>& m):current_level_(0u){
		base_cph_ = m;
	}
	MRadaptativeCMap3(const MRCmap3& mr2):MRCmap3(mr2.base_cph_){
		this->current_level(mr2.current_level());
	}
	MRadaptativeCMap3(const MRCmap3& mr2,uint32 l):MRCmap3(mr2.base_cph_){
		this->current_level(l);
	}

	inline const CPH3<CMap3>& cph() const {return *base_cph_;}
	inline CPH3<CMap3>& cph(){return *base_cph_;}
	
	inline const CMap3& mesh() const {return cph().mesh();}
	inline CMap3& mesh(){return cph().mesh();}
};
*/
}


#endif // CGOGN_CORE_TYPES_CMAP_MR_CMAP3_ADAPTATIVE_H
