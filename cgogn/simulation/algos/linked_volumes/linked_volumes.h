#ifndef CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
#define CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_solver_multiresolution.h>
#include <cgogn/simulation/algos/linked_volumes/cutting_tools.h>

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

	MAP* m_;
	std::shared_ptr<Attribute<Vec3>> centroid_;
	std::shared_ptr<Attribute<double>> distance_plan_;

public:
	Linked_volumes() : m_(nullptr), centroid_(nullptr), distance_plan_(nullptr)
	{
	}

	void init_mesh(MAP* m)
	{
		m_ = m;
		centroid_ = get_attribute<Vec3, Volume>(*m, "centroid");
		if (centroid_ == nullptr)
			centroid_ = add_attribute<Vec3, Volume>(*m, "centroid");
		distance_plan_ = get_attribute<double, Volume>(*m, "distance_plan");
		if (distance_plan_ == nullptr)
			distance_plan_ = add_attribute<double, Volume>(*m, "distance_plan");
	}

	template <typename FUNC>
	void compute_cut_plan(Vec3 dir_plan, double w, Attribute<Vec3>* pos, const FUNC& callback_vertices,
						  bool compute_centroid = true)
	{
		// Ajout plan de coupe en parametre
		if (compute_centroid)
		{
			geometry::compute_centroid<Vec3, Volume>(*m_, pos, centroid_.get());
		}
		parallel_foreach_cell(*m_, [&](Volume v) -> bool {
			value<double>(*m_, this->distance_plan_.get(), v) =
				dir_plan.dot(value<Vec3>(*m_, this->centroid_.get(), v));
			return true;
		});
		std::vector<Face> face_vect;
		foreach_cell(*m_, [&](Face f) -> bool {
			if (is_incident_to_boundary(*m_, f))
			{
				return true;
			}
			double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(f.dart)) - w;
			double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, f.dart))) - w;
			if (v1 * v2 < 0)
			{
				face_vect.push_back(f);
			}
			return true;
		});

		// unsew faces
		for (auto f : face_vect)
		{
			unsew_volume(*m_, f, callback_vertices, true);
		}
	}

	template <typename FUNC>
	void compute_cut_plan_in_framework(Vec3 dir_plan, double w, Attribute<Vec3>* pos, const FUNC& callback_vertices,
									   Simulation_solver_multiresolution<MAP>* ssm, bool compute_centroid = true)
	{
		// Ajout plan de coupe en parametre
		if (compute_centroid)
		{
			geometry::compute_centroid<Vec3, Volume>(*m_, pos, centroid_.get());
		}
		parallel_foreach_cell(*m_, [&](Volume v) -> bool {
			value<double>(*m_, this->distance_plan_.get(), v) =
				dir_plan.dot(value<Vec3>(*m_, this->centroid_.get(), v));
			return true;
		});
		std::vector<Face> face_vect;
		foreach_cell(*m_, [&](Face f) -> bool {
			if (is_incident_to_boundary(*m_, f))
			{
				return true;
			}
			double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(f.dart)) - w;
			double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, f.dart))) - w;
			if (v1 * v2 < 0)
			{
				face_vect.push_back(f);
			}
			return true;
		});

		// unsew faces
		for (auto f : face_vect)
		{
			unsew_volume(*m_, f, callback_vertices, true);
		}
	}
};

template <>
template <typename FUNC>
inline void Linked_volumes<EMR_Map3_Adaptative>::compute_cut_plan(Vec3 dir_plan, double w, Attribute<Vec3>* pos,
																  const FUNC& callback_vertices, bool compute_centroid)
{
	// Ajout plan de coupe en parametre
	if (compute_centroid)
	{
		geometry::compute_centroid<Vec3, Volume>(*m_, pos, centroid_.get());
	}
	parallel_foreach_cell(*m_, [&](Volume v) -> bool {
		value<double>(*m_, this->distance_plan_.get(), v) = dir_plan.dot(value<Vec3>(*m_, this->centroid_.get(), v));
		return true;
	});
	CellMarker<EMR_Map3_Adaptative, Face> face_marker(*m_);
	std::vector<Face> face_vect;
	foreach_cell(*m_, [&](Face f) -> bool {
		if (is_incident_to_boundary(*m_, f))
		{
			return true;
		}
		double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(f.dart)) - w;
		double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, f.dart))) - w;
		if (v1 * v2 < 0)
		{
			face_vect.push_back(f);
		}
		face_marker.mark(f);
		return true;
	});
	CellMarker<EMR_Map3_Adaptative, Volume> vol_marker(*m_);
	std::vector<Volume> volume_vect;
	std::vector<Volume> vect_new_volume;

	while (!face_vect.empty())
	{
		for (auto f : face_vect)
		{
			if (!vol_marker.is_marked(Volume(f.dart)))
			{
				vol_marker.mark(Volume(f.dart));
				volume_vect.push_back(Volume(f.dart));
			}
			if (!vol_marker.is_marked(Volume(phi3(*m_, f.dart))))
			{
				vol_marker.mark(Volume(phi3(*m_, f.dart)));
				volume_vect.push_back(Volume(phi3(*m_, f.dart)));
			}
		}
		for (auto v : volume_vect)
		{
			foreach_incident_vertex(*m_, v, [&](Vertex w) -> bool {
				vect_new_volume.push_back(Volume(w.dart));
				return true;
			});
			m_->activate_volume_subdivision(v);
		}
		geometry::compute_centroid<Vec3, Volume>(*m_, pos, centroid_.get());
		parallel_foreach_cell(*m_, [&](Volume v) -> bool {
			value<double>(*m_, this->distance_plan_.get(), v) =
				dir_plan.dot(value<Vec3>(*m_, this->centroid_.get(), v));
			return true;
		});
		face_vect.clear();

		foreach_cell(*m_, [&](Face f) -> bool {
			if (is_incident_to_boundary(*m_, f) || face_marker.is_marked(f))
			{
				return true;
			}
			double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(f.dart)) - w;
			double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, f.dart))) - w;
			if (v1 * v2 < 0)
			{
				face_vect.push_back(f);
			}
			face_marker.mark(f);
			return true;
		});
	}

	face_vect.clear();
	foreach_cell(*m_, [&](Face f) -> bool {
		if (is_incident_to_boundary(*m_, f))
		{
			return true;
		}

		Dart y = m_->face_youngest_dart(f.dart);
		double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(y)) - w;
		double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, y))) - w;
		if (v1 * v2 < 0)
		{

			face_vect.push_back(f);
		}
		return true;
	});

	// unsew faces
	for (auto f : face_vect)
	{
		unsew_volume(*m_, f, callback_vertices, true);
	}
}

template <>
template <typename FUNC>
inline void Linked_volumes<EMR_Map3_Adaptative>::compute_cut_plan_in_framework(
	Vec3 dir_plan, double w, Attribute<Vec3>* pos, const FUNC& callback_vertices,
	Simulation_solver_multiresolution<EMR_Map3_Adaptative>* ssm, bool compute_centroid)
{
	init_mesh(ssm->topology_->get_copy());

	// Ajout plan de coupe en parametre
	if (compute_centroid)
	{
		geometry::compute_centroid<Vec3, Volume>(*m_, pos, centroid_.get());
	}
	parallel_foreach_cell(*m_, [&](Volume v) -> bool {
		value<double>(*m_, this->distance_plan_.get(), v) = dir_plan.dot(value<Vec3>(*m_, this->centroid_.get(), v));
		return true;
	});
	CellMarker<EMR_Map3_Adaptative, Face> face_marker(*m_);
	std::vector<Face> face_vect;
	foreach_cell(*m_, [&](Face f) -> bool {
		if (is_incident_to_boundary(*m_, f))
		{
			return true;
		}
		double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(f.dart)) - w;
		double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, f.dart))) - w;
		if (v1 * v2 < 0)
		{
			face_vect.push_back(f);
		}
		face_marker.mark(f);
		return true;
	});
	CellMarker<EMR_Map3_Adaptative, Volume> vol_marker(*m_);
	std::vector<Volume> volume_vect;
	std::vector<Volume> vect_new_volume;

	while (!face_vect.empty())
	{
		for (auto f : face_vect)
		{
			if (!vol_marker.is_marked(Volume(f.dart)))
			{
				vol_marker.mark(Volume(f.dart));
				volume_vect.push_back(Volume(f.dart));
			}
			if (!vol_marker.is_marked(Volume(phi3(*m_, f.dart))))
			{
				vol_marker.mark(Volume(phi3(*m_, f.dart)));
				volume_vect.push_back(Volume(phi3(*m_, f.dart)));
			}
		}
		for (auto v : volume_vect)
		{
			/*foreach_incident_vertex(*m_, v, [&](Vertex w) -> bool {
				vect_new_volume.push_back(Volume(w.dart));
				return true;
			});*/
			if (m_->volume_level(v.dart) < m_->maximum_level_ - 1)
			{
				m_->activate_volume_subdivision(v);
			}
		}
		geometry::compute_centroid<Vec3, Volume>(*m_, pos, centroid_.get());
		parallel_foreach_cell(*m_, [&](Volume v) -> bool {
			value<double>(*m_, this->distance_plan_.get(), v) =
				dir_plan.dot(value<Vec3>(*m_, this->centroid_.get(), v));
			return true;
		});
		face_vect.clear();

		foreach_cell(*m_, [&](Face f) -> bool {
			if (is_incident_to_boundary(*m_, f) || face_marker.is_marked(f))
			{
				return true;
			}
			double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(f.dart)) - w;
			double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, f.dart))) - w;
			if (v1 * v2 < 0)
			{
				face_vect.push_back(f);
			}
			face_marker.mark(f);
			return true;
		});
	}

	face_vect.clear();
	foreach_cell(*m_, [&](Face f) -> bool {
		if (is_incident_to_boundary(*m_, f))
		{
			return true;
		}

		Dart y = m_->face_youngest_dart(f.dart);
		double v1 = value<double>(*m_, this->distance_plan_.get(), Volume(y)) - w;
		double v2 = value<double>(*m_, this->distance_plan_.get(), Volume(phi3(*m_, y))) - w;
		if (v1 * v2 < 0)
		{

			face_vect.push_back(f);
		}
		return true;
	});

	ssm->update_tree_volume(*m_, pos);
	std::cout << "debut decoupe topo" << std::endl;

	std::clock_t start;
	double duration;
	start = std::clock();
	// unsew faces
	for (auto f : face_vect)
	{
		unsew_volume(*m_, f, callback_vertices, true);
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\033[1;32m tmps decoupe topo : \033[0m" << duration << std::endl;

	ssm->sc_coarse_->update_topo(*ssm->coarse_meca_mesh_, {});
	ssm->sc_->update_topo(*ssm->mecanical_mesh_, {});
	ssm->sc_fine_->update_topo(*ssm->fine_meca_mesh_, {});
}

} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_LINKED_VOLUMES_LINKED_VOLUMES_H
