#ifndef CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
#define CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class Linked_volumes
{
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MAP>::Vertex;
    using Face = typename mesh_traits<MAP>::Face;
    using Volume = typename mesh_traits<MAP>::Volume;

    MAP& m_;
    std::shared_ptr<Attribute<Vec3>> centroid_;

public:

    Linked_volumes(MAP& m):m_(m){
        centroid_ = get_attribute<Vec3, Volume>(m, "centroid");
        if (centroid_ == nullptr)
            centroid_ = add_attribute<Vec3, Volume>(m, "centroid");
    }

    void compute_cut(Attribute<Vec3> pos,bool compute_centroid = true){
        //Ajout plan de coupe en parametre
        if(compute_centroid){
           geometry::compute_centroid<Vec3,Volume>(m_,pos,centroid_.get());
        }
        foreach_cell(m_,[&](Face f){
            //Si pas à la frontière
            //Check si centroid des volumes de chaque coté du plan
            //Si oui sépartion de la face
        })
    }

};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
