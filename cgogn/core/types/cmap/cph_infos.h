#ifndef CPH_INFOS_H
#define CPH_INFOS_H

#include <cgogn/core/types/cmap/cph3.h>
#include <cgogn/core/types/cmap/phi.h>

#include <cgogn/core/functions/traversals/face.h>


namespace cgogn{
/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

template<typename MESH>
auto edge_level(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,uint32>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	return std::max(m.dart_level(d), m.dart_level(phi1(m, d)));
}

template<typename MESH>
auto edge_youngest_dart(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,Dart>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	Dart d2 = phi2(m, d);
	if (m.dart_level(d) > m.dart_level(d2))
		return d;
	return d2;
}

template<typename MESH>
auto edge_is_subdivided(const MESH& m,Dart d) 
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,bool>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	if (m.current_level_ == m.maximum_level_)
		return false;

	Dart d1 = phi1(m, d);
	MESH m2(m);
	m2.current_level_++;
	Dart d1_l = phi1(m2, d);
	if (d1 != d1_l)
		return true;
	else
		return false;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

template<typename MESH>
auto face_level(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,uint32>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	if (m.current_level_ == 0)
		return 0;

	Dart it = d;
	Dart old = it;
	uint32 l_old = m.dart_level(old);
	uint32 fLevel = edge_level(m,it);
	do
	{
		it = phi1(m, it);
		uint32 dl = m.dart_level(it);

		// compute the oldest dart of the face in the same time
		if (dl < l_old)
		{
			old = it;
			l_old = dl;
		}
		uint32 l = edge_level(m,it);
		fLevel = l < fLevel ? l : fLevel;
	} while (it != d);

	MESH m2(m);
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

	return fLevel;
}

template<typename MESH>
auto face_oldest_dart(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,Dart>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	Dart it = d;
	Dart oldest = it;
	uint32 l_old = m.dart_level(oldest);
	do
	{
		uint32 l = m.dart_level(it);
		if (l == 0)
			return it;

		if (l < l_old)
		//		if(l < l_old || (l == l_old && it < oldest))
		{
			oldest = it;
			l_old = l;
		}
		it = phi1(m, it);
	} while (it != d);

	return oldest;
}

template<typename MESH>
auto face_youngest_dart(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,Dart>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	Dart it = d;
	Dart youngest = it;
	uint32 l_young = m.dart_level(youngest);
	do
	{
		uint32 l = m.dart_level(it);
		if (l == m.current_level_)
			return it;

		if (l > l_young)
		//		if(l < l_young || (l == l_young && it < youngest))
		{
			youngest = it;
			l_young = l;
		}
		it = phi1(m, it);
	} while (it != d);

	return youngest;
}

template<typename MESH>
auto face_origin(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,Dart>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	Dart p = d;
	uint32 pLevel = m.dart_level(p);
	CPH3 m2(m);
	do
	{
		m2.current_level_ = pLevel;
		p = face_oldest_dart(m2,p);
		pLevel = m2.dart_level(p);
	} while (pLevel > 0);
	return p;
}

template<typename MESH>
auto face_is_subdivided(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,bool>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");
	uint32 fLevel = face_level(m,d);
	if (fLevel < m.current_level_)
		return false;

	bool subd = false;
	MESH m2(m);
	m2.current_level_++;
	if (m2.dart_level(phi1(m2, d)) == m.current_level_ && m2.edge_id(phi1(m2, d)) != m2.edge_id(d))
		subd = true;
	return subd;
}

template<typename MESH>
auto face_is_subdivided_once(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,bool>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	uint32 fLevel = face_level(m,d);
	if (fLevel < m.current_level_)
		return false;

	uint32 degree = 0;
	bool subd = false;
	bool subdOnce = true;
	Dart fit = d;
	CPH3 m2(m),m3(m);
	m2.current_level_ = m.current_level_ + 1;
	m3.current_level_ = m.current_level_ + 2;
	do
	{
		if (m.dart_level(phi1(m2, fit)) == m.current_level_ && m.edge_id(phi1(m2, fit)) != m.edge_id(fit))
		{
			subd = true;
			if (m.dart_level(phi1(m3, fit)) == m.current_level_ && m.edge_id(phi1(m3, fit)) != m.edge_id(fit))
				subdOnce = false;
		}
		++degree;
		fit = phi1(m, fit);
	} while (subd && subdOnce && fit != d);

	if (degree == 3 && subd)
	{
		Dart cf = phi2(m2, phi1(m2, d));
		if (m.dart_level(phi1(m3, cf)) == m.current_level_ && m.edge_id(phi1(m3, cf)) != m.edge_id(cf))
			subdOnce = false;
	}

	return subd && subdOnce;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

template<typename MESH>
auto volume_level(const MESH& m,Dart d) 
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,uint32>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	// The level of a volume is the
	// minimum of the levels of its faces

	Dart oldest = d;
	uint32 lold = m.dart_level(oldest);
	uint32 vLevel = std::numeric_limits<uint32>::max();

	foreach_incident_face(m, CPH3::CMAP::Volume(d), [&](CPH3::CMAP::Face f) -> bool {
		uint32 fLevel = face_level(m,f.dart);
		vLevel = fLevel < vLevel ? fLevel : vLevel;
		Dart old = face_oldest_dart(m,f.dart);
		if (m.dart_level(old) < lold)
		{
			oldest = old;
			lold = m.dart_level(old);
		}
		return true;
	});

	CPH3 m2(m);
	m2.current_level_ = vLevel;

	uint32 nbSubd = 0;
	Dart it = oldest;
	uint32 eId = m.edge_id(oldest);
	do
	{
		++nbSubd;
		it = phi<121>(m2, it);
	} while (m.edge_id(it) == eId && lold != m.dart_level(it));

	while (nbSubd > 1)
	{
		nbSubd /= 2;
		--vLevel;
	}


	return vLevel;
}

template<typename MESH>
auto volume_oldest_dart(const MESH& m,Dart d) 
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,Dart>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	Dart oldest = d;
	uint32 l_old = m.dart_level(oldest);
	foreach_incident_face(m, CPH3::CMAP::Volume(oldest), [&](CPH3::CMAP::Face f) -> bool {
		Dart old = face_oldest_dart(m,f.dart);
		uint32 l = m.dart_level(old);
		if (l < l_old)
		{
			oldest = old;
			l_old = l;
		}
		return true;
	});

	return oldest;
}

template<typename MESH>
auto volume_youngest_dart(const MESH& m,Dart d)
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,Dart>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	Dart youngest = d;
	uint32 l_young = m.dart_level(youngest);
	foreach_incident_face(m, CPH3::CMAP::Volume(youngest), [&](CPH3::CMAP::Face f) -> bool {
		Dart young = face_youngest_dart(m,f.dart);
		uint32 l = m.dart_level(young);
		if (l > l_young)
		{
			youngest = young;
			l_young = l;
		}
		return true;
	});

	return youngest;
}

template<typename MESH>
auto volume_is_subdivided(const MESH& m,Dart d) 
-> std::enable_if_t<std::is_convertible_v<MESH&, CPH3&>,bool>
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	uint vLevel = volume_level(m,d);
	if (vLevel < m.current_level_)
		return false;

	bool faceAreSubdivided = face_is_subdivided(m,d);

	foreach_incident_face(m, CPH3::CMAP::Volume(d), [&](CPH3::CMAP::Face f) -> bool {
		faceAreSubdivided &= face_is_subdivided(m,f.dart);
		return true;
	});

	bool subd = false;
	CPH3 m2(m);
	m2.current_level_++;
	if (faceAreSubdivided && m.dart_level(phi<112>(m2, d)) == m.current_level_ &&
		m.face_id(phi<112>(m2, d)) != m.face_id(d))
		subd = true;

	return subd;
}

}


#endif // CPH_INFOS_H
