#include "EMR_Map3.h"
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

Dart EMR_Map3::edge_youngest_dart(Dart d) const
{
	Dart it = phi2(*this, d);
	if (m_.dart_level(d) < m_.dart_level(it))
		return d;
	return it;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

Dart EMR_Map3::face_youngest_dart(Dart d) const
{
	Dart it = phi1(*this, d);
	Dart result = d;
	while (it != d)
	{
		if (m_.dart_level(it) < m_.dart_level(result))
			result = it;
		it = phi1(*this, d);
	}

	return result;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

Dart EMR_Map3::volume_youngest_dart(Dart d) const
{
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

} // namespace cgogn
