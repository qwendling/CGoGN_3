#include "EMR_Map3.h"
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

template <>
Dart EMR_Map3::edge_youngest_dart(Dart d) const
{
	Dart it = phi2(*this, d);
	if (this->dart_level(d) < this->dart_level(it))
		return d;
	return it;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

template <>
Dart EMR_Map3::face_youngest_dart(Dart d) const
{
	Dart it = phi1(*this, d);
	Dart result = d;
	while (it != d)
	{
		if (this->dart_level(it) < this->dart_level(result))
			result = it;
		it = phi1(*this, d);
	}

	return result;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

template <>
Dart EMR_Map3::volume_youngest_dart(Dart d) const
{
	Dart result = d;
	foreach_dart_of_orbit(*this, Volume(d), [&](Dart it) -> bool {
		if (this->dart_level(it) < this->dart_level(result))
		{
			result = it;
		}
		return true;
	});
	return result;
}

} // namespace cgogn
