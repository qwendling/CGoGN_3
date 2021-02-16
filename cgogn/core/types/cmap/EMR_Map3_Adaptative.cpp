#include "EMR_Map3_Adaptative.h"

namespace cgogn
{
uint32 EMR_Map3_Adaptative::get_dart_visibility(Dart d) const
{
	return std::min(this->dart_level(d), (*dart_visibility_)[d.index]);
}

} // namespace cgogn
