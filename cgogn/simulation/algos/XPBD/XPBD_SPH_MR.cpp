#include "XPBD_SPH_MR.h"
#include <Eigen/SVD>
#include <cgogn/geometry/algos/centroid.h>
#include <forward_list>

namespace cgogn
{
namespace simulation
{

void XPBD_SPH_Multiresolution::init_particule_MR(MAP& m, Attribute<Vec3>* pos)
{
	MAP m2(m.m_);
	m2.current_level_ = m2.maximum_level_;
	MAP m3(m.m_);
	m3.current_level_ = std::max(0u, m3.maximum_level_ - 1);

	particule_volume_ = add_attribute<Particule_SPH_MR*, Volume>(
		m2, "SPH_particule_constraint_solver_particule_volume_" + std::to_string(id));

	initial_volume_ =
		add_attribute<double, Volume>(m2, "SPH_particule_constraint_solver_initial_volume_" + std::to_string(id));
	initial_centroid_volume_ = add_attribute<Vec3, Volume>(
		m2, "SPH_particule_constraint_solver_initial_centroid_volume_" + std::to_string(id));

	geometry::compute_centroid<Vec3, Volume>(m2, pos, initial_centroid_volume_.get());
	geometry::compute_volume(m2, pos, initial_volume_.get());

	foreach_cell(m2, [&](Volume v) -> bool {
		Particule_SPH_MR* p = new Particule_SPH_MR(value<Vec3>(m2, initial_centroid_volume_.get(), v),
												   DENSITY_SPH * value<double>(m2, initial_volume_.get(), v));
		p->initial_volume_ = value<double>(m2, initial_volume_.get(), v);
		value<Particule_SPH_MR*>(m2, particule_volume_.get(), v) = p;
		return true;
	});
	CellMarker<MAP, Volume> marker(m2);
	while (m2.current_level_ > 0)
	{
		foreach_cell(m3, [&](Volume v) -> bool {
			if (m3.volume_level(v.dart) == m2.volume_level(v.dart))
				return true;
			Particule_SPH_MR* p = new Particule_SPH_MR(Vec3(0, 0, 0), 0);
			value<Particule_SPH_MR*>(m3, particule_volume_.get(), v) = p;
			foreach_dart_of_orbit(m3, v, [&](Dart d) -> bool {
				if (marker.is_marked(Volume(d)))
					return true;
				marker.mark(Volume(d));
				p->child_.push_front(value<Particule_SPH_MR*>(m2, particule_volume_.get(), Volume(d)));
				return true;
			});
			double masse = 0;
			double volume = 0;
			Vec3 pos = Vec3::Zero();
			for (auto c : p->child_)
			{
				masse += c->masse_;
				volume += c->initial_volume_;
				pos += c->masse_ * c->initial_position_;
			}
			p->initial_position_ = pos / masse;
			p->current_position_ = pos / masse;
			p->masse_ = masse;
			p->initial_volume_ = volume;
			return true;
		});
		m2.current_level_--;
		m3.current_level_--;
	}
	foreach_cell(m, [&](Volume v) -> bool {
		Particule_SPH_MR* p = value<Particule_SPH_MR*>(m, particule_volume_.get(), v);
		particules_.push_back(p);
		return true;
	});
}
void XPBD_SPH_Multiresolution::compute_neighborhood_Volume(MAP& m, Attribute<Vec3>* pos)
{
	init_particule_MR(m, pos);
	std::function<void(Particule_SPH_MR*)> fn;
	fn = [&](Particule_SPH_MR* p) {
		if (p->child_.empty())
			return;
		for (auto c : p->child_)
		{
			c->h_ = 0.85 * p->h_;
			fn(c);
		}
	};
	parallel_foreach_cell(m, [&](Volume v) -> bool {
		Particule_SPH_MR& p = *value<Particule_SPH_MR*>(m, particule_volume_.get(), v);
		p.initial_volume_ = value<double>(m, initial_volume_.get(), v);
		CellMarker<MAP, Volume> marker(m);
		std::forward_list<Particule_SPH_MR*>& n = p.neighborhood_;
		n.clear();
		double& h = p.h_;
		h = 0;
		Vec3 pos_v1 = p.initial_position_;
		marker.mark(v);
		n.push_front(&p);
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			foreach_incident_volume(m, w, [&](Volume v2) -> bool {
				foreach_incident_vertex(m, v2, [&](Vertex w2) -> bool {
					foreach_incident_volume(m, w2, [&](Volume v3) -> bool {
						if (!marker.is_marked(v3))
						{
							Particule_SPH_MR* p3 = value<Particule_SPH_MR*>(m, particule_volume_.get(), v3);
							Vec3 pos_v3 = p3->initial_position_;
							if (h < 1e-9)
								h = (pos_v1 - pos_v3).norm();
							h = std::min(h, (pos_v1 - pos_v3).norm());
							n.push_front(p3);
							marker.mark(v3);
						}
						return true;
					});
					return true;
				});
				return true;
			});
			return true;
		});
		h *= 4;
		fn(&p);
		return true;
	});
}

void XPBD_SPH_Multiresolution::init_solver(MAP& m, std::shared_ptr<Attribute<Vec3>> pos)
{
	compute_neighborhood_Volume(m, pos.get());
	compute_corrected_matrix();
}

double XPBD_SPH_Multiresolution::Kernel_W(double dist, double h) const
{
	double q = dist / h;
	if (q > 1)
		return 0;
	double m1 = (1.0f - q);
	double m2 = (4.0f * q + 1.0f);
	double h3 = h * h * h;
	double alpha_d = NORMALIZE_TERM * (1 / h3);

	return alpha_d * m1 * m1 * m1 * m1 * m2;
}

XPBD_SPH_Multiresolution::Vec3 XPBD_SPH_Multiresolution::gradient(XPBD_SPH_Multiresolution::Vec3 xij, double dist,
																  double h) const
{
	if (dist < 1e-12)
		return Vec3(0, 0, 0);
	double tmp = 0;
	double q = dist / h;

	if (q <= 1)
	{
		double m_1 = -210 / (M_PI * h * h * h);
		double m2 = 1 - q;
		tmp = m_1 * (1 / (dist * h)) * m2 * m2 * m2;
	}

	return xij * tmp;
}

XPBD_SPH_Multiresolution::Mat3d XPBD_SPH_Multiresolution::Corrected_matrix(Particule_SPH_MR& p, double h) const
{
	Mat3d Li = Mat3d::Zero();
	std::forward_list<Particule_SPH_MR*>& n = p.neighborhood_;
	Vec3 xi0 = p.initial_position_;
	for (Particule_SPH_MR* w : n)
	{
		double init_vol = w->initial_volume_;
		Vec3 xj0 = w->initial_position_;
		Vec3 xji0 = xj0 - xi0;
		Vec3 grad = gradient(xji0, xji0.norm(), h);
		Li -= init_vol * grad * xji0.transpose();
	}
	bool inversible = false;
	Mat3d L;
	Li.computeInverseWithCheck(L, inversible, 1e-9);
	if (!inversible)
		std::cout << "L pas inversible " << std::endl;
	return L;
}

void XPBD_SPH_Multiresolution::compute_corrected_matrix()
{
	for (Particule_SPH_MR* p : particules_)
	{
		p->corrected_matrix_ = Corrected_matrix(*p, p->h_);
	}
}

XPBD_SPH_Multiresolution::Vec3 XPBD_SPH_Multiresolution::Corrected_gradient(Particule_SPH_MR& vi, Particule_SPH_MR& vj,
																			double h) const
{
	Vec3 xi = vi.initial_position_;
	Vec3 xj = vj.initial_position_;
	Vec3 xij = xi - xj;

	return vi.corrected_matrix_ * gradient(xij, xij.norm(), h);
}

void XPBD_SPH_Multiresolution::activate_remove_volume(MAP& m, std::vector<Volume>& list_activate,
													  std::vector<Volume>& list_remove)
{

	std::vector<Volume> new_volumes;
	for (Volume v : list_activate)
	{
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		t->F_ = value<Mat3d>(m, F_, v);
		t->init_cm_ = value<Vec3>(m, init_cm_, v);
		double masse_vol = 0;
		t->cm_ = Vec3(0, 0, 0);
		t->v_cm_ = Vec3(0, 0, 0);
		std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		for (Vertex w : inc_vertices)
		{
			masse_vol += value<double>(m, masse_, w);
			t->cm_ += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
			t->v_cm_ += value<double>(m, masse_, w) * value<Vec3>(m, speed_.get(), w);
		}
		t->cm_ /= masse_vol;
		t->v_cm_ /= masse_vol;
		t->for_each_child([&](tree_volume* c) -> bool {
			new_volumes.push_back(Volume(c->volume_dart));
			return true;
		});
	}
	for (Volume v : list_remove)
	{
		m.disable_volume_subdivision(v, true);
	}
	CellMarkerStore<MAP, Volume> marked_Volumes(m);
	CellMarkerStore<MAP, Vertex> marked_Vertices(m);
	for (Volume v : list_remove)
	{
		marked_Volumes.mark(v);
	}
	for (Volume v : list_activate)
	{
		marked_Volumes.mark(v);
	}
	std::vector<Volume> impacted_volumes;
	for (Volume v : list_activate)
	{
		foreach_adjacent_volume_through_edge(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				impacted_volumes.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}
	for (Volume v : list_remove)
	{
		foreach_adjacent_volume_through_edge(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				impacted_volumes.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}

	std::vector<Volume> volume_need_update;
	for (Volume v : impacted_volumes)
	{
		foreach_incident_vertex(m, v, [&](Vertex inc_vert) -> bool {
			marked_Vertices.mark(inc_vert);
			value<double>(m, masse_, inc_vert) = 0;
			bool res_nested_lambda = true;
			foreach_incident_volume(m, inc_vert, [&](Volume inc_vol) -> bool {
				if (!is_boundary(m, inc_vol.dart))
				{
					if (!marked_Volumes.is_marked(inc_vol))
					{
						volume_need_update.push_back(inc_vol);
						marked_Volumes.mark(inc_vol);
					}
					res_nested_lambda = true;
				}
				return res_nested_lambda;
			});
			return res_nested_lambda;
		});
	}
	for (Volume v : impacted_volumes)
	{
		foreach_adjacent_volume_through_vertex(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				volume_need_update.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}

	for (Volume v : list_activate)
	{
		m.activate_volume_subdivision(v);
	}

	for (Volume v : list_remove)
	{
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<double>(m, masse_, w) = 0;
			return true;
		});
	}
	for (Volume v : new_volumes)
	{
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<double>(m, masse_, w) = 0;
			return true;
		});
	}

	auto compute_masse = [&](std::vector<Volume> list_volumes) {
		for (Volume v : list_volumes)
		{
			std::vector<Vertex> vertices;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				vertices.push_back(w);
				return true;
			});
			double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
			for (auto w : vertices)
			{
				value<double>(m, masse_, w) += masse;
			}
		}
	};

	compute_masse(new_volumes);
	compute_masse(list_remove);
	compute_masse(impacted_volumes);
	for (Volume v : volume_need_update)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			if (marked_Vertices.is_marked(w))
				value<double>(m, masse_, w) += masse;
		}
	}
	auto fn = [&](const std::vector<Volume>& l_vol) {
		for (Volume v : l_vol)
		{
			double masse = 0;
			Vec3 cm = Vec3(0, 0, 0);
			std::vector<Vertex> inc_vertices;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				masse += value<double>(m, masse_, w);
				cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
				inc_vertices.push_back(w);
				return true;
			});
			value<Vec3>(m, init_cm_, v) = cm / masse;
			Mat3d Q = Mat3d::Zero();
			for (Vertex w : inc_vertices)
			{
				Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
				double masse = value<double>(m, masse_, w);
				Q += masse * init_r_i * init_r_i.transpose();
			}
			double s = 1.0 / Q.sum();
			value<double>(m, s_, v) = s;
			value<Mat3d>(m, inv_Q_, v) = (s * Q).eval().inverse();
		};
	};

	fn(list_remove);
	fn(new_volumes);
	fn(impacted_volumes);
	fn(volume_need_update);
	for (Volume v : new_volumes)
	{
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		t = t->pere;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<Vec3>(m, speed_, w) = t->v_cm_;
			Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - t->init_cm_;
			value<Vec3>(m, pos_, w) = t->cm_ + t->F_ * init_r_i;
			return true;
		});
	}
}

void XPBD_SPH_Multiresolution::activate_volume(MAP& m, std::vector<Volume>& list_Volumes)
{
	CellMarkerStore<MAP, Volume> marked_Volumes(m);
	for (Volume v : list_Volumes)
	{
		marked_Volumes.mark(v);
	}
	std::vector<Volume> impacted_volumes;
	for (Volume v : list_Volumes)
	{
		foreach_adjacent_volume_through_edge(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				impacted_volumes.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}

	std::vector<Volume> new_volumes;
	for (Volume v : list_Volumes)
	{
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		t->F_ = value<Mat3d>(m, F_, v);
		t->init_cm_ = value<Vec3>(m, init_cm_, v);
		double masse_vol = 0;
		t->cm_ = Vec3(0, 0, 0);
		t->v_cm_ = Vec3(0, 0, 0);
		std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		for (Vertex w : inc_vertices)
		{
			masse_vol += value<double>(m, masse_, w);
			t->cm_ += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
			t->v_cm_ += value<double>(m, masse_, w) * value<Vec3>(m, speed_.get(), w);
		}
		t->cm_ /= masse_vol;
		t->v_cm_ /= masse_vol;
		t->for_each_child([&](tree_volume* c) -> bool {
			new_volumes.push_back(Volume(c->volume_dart));
			return true;
		});
	}

	for (Volume v : impacted_volumes)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			value<double>(m, masse_, w) -= masse;
		}
	}
	for (Volume v : list_Volumes)
	{
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<double>(m, masse_, w) = 0;
			return true;
		});
	}
	std::vector<Volume> volume_need_update;
	for (Volume v : impacted_volumes)
	{
		foreach_adjacent_volume_through_vertex(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				volume_need_update.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}
	for (Volume v : list_Volumes)
	{
		m.activate_volume_subdivision(v);
	}
	for (Volume v : new_volumes)
	{
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<double>(m, masse_, w) = 0;
			return true;
		});
	}

	for (Volume v : new_volumes)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			value<double>(m, masse_, w) += masse;
		}
	}
	for (Volume v : impacted_volumes)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			value<double>(m, masse_, w) += masse;
		}
	}

	auto fn = [&](const std::vector<Volume>& l_vol) {
		for (Volume v : l_vol)
		{
			double masse = 0;
			Vec3 cm = Vec3(0, 0, 0);
			std::vector<Vertex> inc_vertices;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				masse += value<double>(m, masse_, w);
				cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
				inc_vertices.push_back(w);
				return true;
			});
			value<Vec3>(m, init_cm_, v) = cm / masse;
			Mat3d Q = Mat3d::Zero();
			for (Vertex w : inc_vertices)
			{
				Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
				double masse = value<double>(m, masse_, w);
				Q += masse * init_r_i * init_r_i.transpose();
			}
			double s = 1.0 / Q.sum();
			value<double>(m, s_, v) = s;
			value<Mat3d>(m, inv_Q_, v) = (s * Q).eval().inverse();
		};
	};

	fn(new_volumes);
	fn(impacted_volumes);
	fn(volume_need_update);

	for (Volume v : new_volumes)
	{
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		t = t->pere;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<Vec3>(m, speed_, w) = t->v_cm_;
			Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - t->init_cm_;
			value<Vec3>(m, pos_, w) = t->cm_ + t->F_ * init_r_i;
			return true;
		});
	}
}

void XPBD_SPH_Multiresolution::remove_volume(MAP& m, std::vector<Volume>& list_Volumes)
{
	for (Volume v : list_Volumes)
	{
		m.disable_volume_subdivision(v, true);
	}
	CellMarkerStore<MAP, Volume> marked_Volumes(m);
	CellMarkerStore<MAP, Vertex> marked_Vertices(m);
	for (Volume v : list_Volumes)
	{
		marked_Volumes.mark(v);
	}
	std::vector<Volume> impacted_volumes;
	for (Volume v : list_Volumes)
	{
		foreach_adjacent_volume_through_edge(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				impacted_volumes.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}

	for (Volume v : impacted_volumes)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			marked_Vertices.mark(w);
			value<double>(m, masse_, w) = 0;
			return true;
		});
	}
	for (Volume v : list_Volumes)
	{
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<double>(m, masse_, w) = 0;
			return true;
		});
	}
	std::vector<Volume> volume_need_update;
	for (Volume v : impacted_volumes)
	{
		foreach_adjacent_volume_through_vertex(m, v, [&](Volume w) -> bool {
			if (!marked_Volumes.is_marked(w))
			{
				volume_need_update.push_back(w);
				marked_Volumes.mark(w);
			}
			return true;
		});
	}

	for (Volume v : list_Volumes)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			value<double>(m, masse_, w) += masse;
		}
	}
	for (Volume v : impacted_volumes)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			value<double>(m, masse_, w) += masse;
		}
	}

	for (Volume v : volume_need_update)
	{
		std::vector<Vertex> vertices;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vertices.push_back(w);
			return true;
		});
		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			if (marked_Vertices.is_marked(w))
				value<double>(m, masse_, w) += masse;
		}
	}

	auto fn = [&](const std::vector<Volume>& l_vol) {
		for (Volume v : l_vol)
		{
			double masse = 0;
			Vec3 cm = Vec3(0, 0, 0);
			std::vector<Vertex> inc_vertices;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				masse += value<double>(m, masse_, w);
				cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
				inc_vertices.push_back(w);
				return true;
			});
			value<Vec3>(m, init_cm_, v) = cm / masse;
			Mat3d Q = Mat3d::Zero();
			for (Vertex w : inc_vertices)
			{
				Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
				double masse = value<double>(m, masse_, w);
				Q += masse * init_r_i * init_r_i.transpose();
			}
			double s = 1.0 / Q.sum();
			value<double>(m, s_, v) = s;
			value<Mat3d>(m, inv_Q_, v) = (s * Q).eval().inverse();
		};
	};

	fn(list_Volumes);
	fn(impacted_volumes);
	fn(volume_need_update);
}

void XPBD_SPH_Multiresolution::update_topo(MAP& m)
{
	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		value<double>(m, masse_, v) = 0;
		return true;
	});
	foreach_cell(m, [&](Volume v) -> bool {
		std::vector<Vertex>& vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		vertices.clear();
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<std::vector<Vertex>>(m, inc_vertices_.get(), v).push_back(w);
			return true;
		});

		double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
		for (auto w : vertices)
		{
			value<double>(m, masse_, w) += masse;
		}
		return true;
	});

	foreach_cell(m, [&](Volume v) -> bool {
		double masse = 0;
		Vec3 cm = Vec3(0, 0, 0);
		std::vector<Vertex>& vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		for (Vertex w : vertices)
		{
			masse += value<double>(m, masse_, w);
			cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
		}
		value<Vec3>(m, init_cm_, v) = cm / masse;
		return true;
	});

	foreach_cell(m, [&](Volume v) -> bool {
		Mat3d Q = Mat3d::Zero();
		std::vector<Vertex>& vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		for (Vertex w : vertices)
		{
			Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			double masse = value<double>(m, masse_, w);
			Q += masse * init_r_i * init_r_i.transpose();
		}
		double s = 1.0 / Q.sum();
		value<double>(m, s_, v) = s;
		value<Mat3d>(m, inv_Q_, v) = (s * Q).eval().inverse();
		return true;
	});
}

void XPBD_SPH_Multiresolution::compute_deformation_gradient_particule(Particule_SPH_MR* p)
{
	double h = p->h_;
	Mat3d Ftemp = Mat3d::Identity();
	std::forward_list<Particule_SPH_MR*>& n = p->neighborhood_;
	Vec3 xi = p->current_position_;
	Vec3 xi0 = p->initial_position_;
	for (Particule_SPH_MR* w : n)
	{
		double Vj0 = w->initial_volume_;
		Vec3 xj = w->current_position_;
		Vec3 xji = xj - xi;
		Vec3 W = Corrected_gradient(*p, *w, h);
		Ftemp += Vj0 * (xji)*W.transpose();
	}
	p->deformation_gradient_ = Ftemp;
}

void XPBD_SPH_Multiresolution::compute_deformation_gradient_particules()
{
	for (Particule_SPH_MR* p : particules_)
	{
		compute_deformation_gradient_particule(p);
	};
}

void XPBD_SPH_Multiresolution::constraint_Neo_Hookean_H(MAP& m, Particule_SPH_MR& p, double h)
{
	double h_smooth = p.h_;
	std::forward_list<Particule_SPH_MR*>& n = p.neighborhood_;
	for (Particule_SPH_MR* w : n)
	{
		compute_deformation_gradient_particule(w);
	}

	Mat3d F = p.deformation_gradient_;

	for (Particule_SPH_MR* w : n)
	{
		double det_F = w->deformation_gradient_.determinant();
		w->CH = det_F - (1. + LAME_MU / LAME_LAMBDA);
	}
	// Compute Volume
	// double Ve = geometry::volume(m, v, pos_.get());
	double Ve = fabs(p.deformation_gradient_.determinant() * p.initial_volume_);

	// Compute alpha_H = 1/(LAMBDA*V)
	double alpha_h = 1.0 / (LAME_LAMBDA * Ve);

	// Compute denum
	double denum = 0;
	const Mat3d T_inv_Q = inv_Q.transpose();
	Mat3d tmp;
	tmp.col(0) = F.col(1).cross(F.col(2));
	tmp.col(1) = F.col(2).cross(F.col(0));
	tmp.col(2) = F.col(0).cross(F.col(1));
	tmp = tmp * T_inv_Q;
	Vec3 GC;
	for (Particule_SPH_MR* w : n)
	{
		double V_j = w->masse_;
		Vec3 xi0 = p.initial_position_;
		Vec3 xj0 = w->initial_position_;
		GC += V_j * w->CH * gradient((xi0 - xj0), (xi0 - xj0).norm(), h_smooth);
	}
	for (Vertex w : inc_vertices)
	{
		// Compute Grad_C_i
		double m_i = value<double>(m, masse_, w);
		Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
		Vec3 GC = m_i * tmp * init_r_i;
		GC = value<double>(m, s_, v) * GC;
		value<Vec3>(m, Grad_C_i_, w) = GC;
		denum += 1.0 / m_i * GC.squaredNorm();
	}
	denum += alpha_h / (h * h);
	// Compute lambda
	double lambda = -C / denum;
	for (Vertex w : inc_vertices)
	{
		if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), w))
		{
			continue;
		}
		double m_i = value<double>(m, masse_, w);
		Vec3 delta_x = lambda * (1 / m_i) * value<Vec3>(m, Grad_C_i_, w);
		value<Vec3>(m, pos_.get(), w) += delta_x;
	}
}
void XPBD_SPH_Multiresolution::constraint_Neo_Hookean_D(MAP& m, Particule_SPH_MR& p, double h)
{
	// Compute center of mass
	double masse_vol = 0;
	Vec3 cm = Vec3(0, 0, 0);
	std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
	for (Vertex w : inc_vertices)
	{
		masse_vol += value<double>(m, masse_, w);
		cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
	}
	cm /= masse_vol;

	// Compute P
	Mat3d P = Mat3d::Zero();
	for (Vertex w : inc_vertices)
	{
		double m_i = value<double>(m, masse_, w);
		Vec3 r_i = value<Vec3>(m, pos_.get(), w) - cm;
		Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
		P.noalias() += m_i * r_i * init_r_i.transpose();
	}

	Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
	P = (value<double>(m, s_, v) * P).eval();
	// compute F = P*Q^-1
	Mat3d F = P * inv_Q;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (i == j)
			{
				if (fabs(F(i, j) - 1) < EPS)
					F(i, j) = 1;
			}
			else
			{
				if (fabs(F(i, j)) < EPS)
					F(i, j) = 0;
			}
		}
	}
	F = value<Mat3d>(m, F_, v);

	// Compute C  = sqrt(tr(F^T*F))
	double C = sqrt((F.transpose() * F).trace());

	// Compute Volume
	// double Ve = geometry::volume(m, v, pos_.get());
	double Ve = fabs(F.determinant() * value<double>(m, init_volume_, v));

	// Compute alpha_D = 1/(MU*V)
	double alpha_d = 1.0 / (LAME_MU * Ve);

	// Compute denum
	double denum = 0;
	const Mat3d T_inv_Q = inv_Q.transpose();
	double r = sqrt(F.col(0).squaredNorm() + F.col(1).squaredNorm() + F.col(2).squaredNorm());
	const Mat3d tmp = F * T_inv_Q;
	for (Vertex w : inc_vertices)
	{
		// Compute Grad_C_i
		double m_i = value<double>(m, masse_, w);
		Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
		Vec3 GC = m_i / r * tmp * init_r_i;
		GC = value<double>(m, s_, v) * GC;
		value<Vec3>(m, Grad_C_i_, w) = GC;
		denum += 1.0 / m_i * GC.squaredNorm();
	}
	denum += alpha_d / (h * h);
	// Compute lambda
	double lambda = -C / denum;

	for (Vertex w : inc_vertices)
	{
		if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), w))
		{
			continue;
		}
		double m_i = value<double>(m, masse_, w);
		Vec3 delta_x = lambda * (1. / m_i) * value<Vec3>(m, Grad_C_i_, w);
		value<Vec3>(m, pos_.get(), w) += delta_x;
	}
}
void XPBD_SPH_Multiresolution::constraint_Zero_Energy(MAP& m, Particule_SPH_MR& p, double)
{
	// Compute center of mass
	double masse_vol = 0;
	Vec3 cm = Vec3(0, 0, 0);
	std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
	for (Vertex w : inc_vertices)
	{
		masse_vol += value<double>(m, masse_, w);
		cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
	}
	cm /= masse_vol;

	// Compute P
	Mat3d P = Mat3d::Zero();
	for (Vertex w : inc_vertices)
	{
		double m_i = value<double>(m, masse_, w);
		Vec3 r_i = value<Vec3>(m, pos_.get(), w) - cm;
		Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
		P.noalias() += m_i * r_i * init_r_i.transpose();
	}

	Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
	P = (value<double>(m, s_, v) * P).eval();
	// compute F = P*Q^-1
	Mat3d F = P * inv_Q;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (i == j)
			{
				if (fabs(F(i, j) - 1) < EPS)
					F(i, j) = 1;
			}
			else
			{
				if (fabs(F(i, j)) < EPS)
					F(i, j) = 0;
			}
		}
	}

	// F = value<Mat3d>(m, F_, v);
	value<Mat3d>(m, F_, v) = F;
	value<Vec3>(m, centroid_, v) = cm;
	// Pour l'erreur et la Visu
	// double det_F = F.determinant();
	/*double energy = (det_F - (1 + LAME_MU / LAME_LAMBDA));
	energy *= energy;
	energy = LAME_LAMBDA / 2.0 * energy + LAME_MU / 2.0 * ((F.transpose() * F).trace() - 3);
	double vol = value<double>(m, init_volume_, v) * det_F;
	value<double>(m, Det_F_Volume_, v) = log(1 + fabs(vol * energy));*/

	for (Vertex w : inc_vertices)
	{
		if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), w))
		{
			continue;
		}
		Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
		value<Vec3>(m, pos_.get(), w) = cm + F * init_r_i;
	}
}
void XPBD_SPH_Multiresolution::applyDamping(MAP& m, Volume v, double damping_coeff, double time_step)
{
	Vec3 x_cm = Vec3::Zero();
	Vec3 v_cm = Vec3::Zero();
	Vec3 L = Vec3::Zero();
	Mat3d I = Mat3d::Zero();
	double sm_i = 0;
	foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
		double m_i = value<double>(m, masse_, w);
		x_cm += m_i * value<Vec3>(m, pos_.get(), w);
		v_cm += m_i * value<Vec3>(m, speed_, w);
		sm_i += m_i;
		return true;
	});
	x_cm /= sm_i;
	v_cm /= sm_i;
	foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
		Vec3 r_i = value<Vec3>(m, pos_.get(), w) - x_cm;
		double m_i = value<double>(m, masse_, w);
		L += r_i.cross(m_i * value<Vec3>(m, speed_, w));
		Mat3d Ri;
		Ri << 0, -r_i(2), r_i(1), r_i(2), 0, -r_i(0), -r_i(1), r_i(0), 0;
		I += Ri * Ri.transpose() * m_i;
		return true;
	});
	Vec3 omega = I.inverse() * L;
	foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
		Vec3 r_i = value<Vec3>(m, pos_.get(), w) - x_cm;
		Vec3 new_v_i = v_cm + omega.cross(r_i);
		value<Vec3>(m, speed_, w) += std::min(damping_coeff * time_step, 1.0) * (new_v_i - value<Vec3>(m, speed_, w));
		return true;
	});
}

void XPBD_SPH_Multiresolution::solve_surface(MAP& m, MAP& geom, Volume v)
{
	bool result = false;
	foreach_dart_of_orbit(m, v, [&m, &result](Dart d) -> bool {
		if (is_boundary(m, phi3(m, d)))
		{
			result = true;
			return false;
		}
		return true;
	});
	if (!result)
		return;
	// Compute center of mass
	/*double masse_vol = 0;
	Vec3 cm = Vec3(0, 0, 0);
	std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
	for (Vertex w : inc_vertices)
	{
		masse_vol += value<double>(m, masse_, w);
		cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
	}
	cm /= masse_vol;*/
	Vec3 cm = value<Vec3>(m, centroid_, v);
	Mat3d F = value<Mat3d>(m, F_, v);

	foreach_incident_vertex(geom, v, [&](Vertex w) -> bool {
		if (!m.vertex_is_visible(w.dart))
		{
			Vec3 init_r_i = value<Vec3>(geom, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
			value<Vec3>(geom, pos_.get(), w) = cm + F * init_r_i;
		}
		return true;
	});
}

#define TEST_ERROR 1
void XPBD_SPH_Multiresolution::compute_error(MAP& m, std::vector<Volume>& volume_activate,
											 std::vector<Volume>& volume_disable)
{
#if TEST_ERROR
	std::forward_list<tree_volume*> list_volume_coarse;
	std::forward_list<tree_volume*> list_volume_fine;
	uint32 nb_fine = 0, nb_coarse = 0;
	static uint32 clock_error = 0;
	clock_error++;

	foreach_cell(m, [&](Volume v) -> bool {
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		if (t->type == CURRENT && t->fils != nullptr)
		{
			list_volume_fine.push_front(t);
			double detF = value<double>(m, this->Det_F_Volume_, Volume(t->volume_dart));
			double error = fabs((detF - 1) * value<double>(m, this->init_volume_, Volume(t->volume_dart)));
			error = fabs(log2(detF));
			value<double>(m, error_volume_, Volume(t->volume_dart)) = error;
			t->error = error;
			nb_fine++;
		}
		if (!t->is_topo && t->pere != nullptr && t->pere->type == COARSE)
		{
			if (t->pere->clock == clock_error)
				return true;
			t->pere->clock = clock_error;
			list_volume_coarse.push_front(t);
			t->pere->error = 0;
			t->pere->for_each_child([&](tree_volume* c) -> bool {
				double detF = value<double>(m, this->Det_F_Volume_, Volume(c->volume_dart));
				double error = fabs((detF - 1) * value<double>(m, this->init_volume_, Volume(c->volume_dart)));
				error = fabs(log2(detF));
				t->pere->error += error;
				return true;
			});
			t->pere->error /= 8;
			nb_coarse++;
		}
		return true;
	});
	list_volume_fine.sort([&](tree_volume* t1, tree_volume* t2) { return t1->error > t2->error; });
	list_volume_coarse.sort([&](tree_volume* t1, tree_volume* t2) { return t1->pere->error < t2->pere->error; });

	for (uint32 i = 0; i < nb_fine / 10; i++)
	{
		tree_volume* t = list_volume_fine.front();
		list_volume_fine.pop_front();
		double error = t->error;
		if (error < 0.25)
			break;

		volume_activate.push_back(Volume(t->volume_dart));
		if (t->pere && t->pere->type != ROOT)
		{
			t->pere->type = NONE;
		}
		t->type = COARSE;
		t->for_each_child([&](tree_volume* c) -> bool {
			c->type = CURRENT;
			return true;
		});
	}

	for (uint32 i = 0; i < nb_coarse / 10; i++)
	{
		tree_volume* t = list_volume_coarse.front();
		list_volume_coarse.pop_front();
		if (t->pere->type != COARSE)
			continue;
		double error = t->error;
		if (error > 1.5)
			break;
		volume_disable.push_back(Volume(t->pere->volume_dart));
		t->pere->type = CURRENT;
		t->pere->for_each_child([&](tree_volume* c) -> bool {
			c->type = NONE;
			return true;
		});
		if (t->pere->pere != nullptr)
		{
			bool result = true;
			t->pere->pere->for_each_child([&](tree_volume* c) -> bool {
				if (c->type != CURRENT)
				{
					result = false;
				}
				return result;
			});
			if (result)
			{
				t->pere->pere->type = COARSE;
			}
		}
	}
#else
	foreach_cell(m, [&](Volume v) -> bool {
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		if (t->type == CURRENT && t->fils != nullptr)
		{
			if (std::rand() / double(RAND_MAX + 1u) < 0.1)
			{
				volume_activate.push_back(Volume(t->volume_dart));
				if (t->pere && t->pere->type != ROOT)
				{
					t->pere->type = NONE;
				}
				t->type = COARSE;
				t->for_each_child([&](tree_volume* c) -> bool {
					c->type = CURRENT;
					return true;
				});
				return true;
			}
		}
		if (!t->is_topo && t->pere != nullptr && t->pere->type == COARSE)
		{
			if (std::rand() / double(RAND_MAX + 1u) < 0.1)
			{
				volume_disable.push_back(Volume(t->pere->volume_dart));
				t->pere->type = CURRENT;
				t->pere->for_each_child([&](tree_volume* c) -> bool {
					c->type = NONE;
					return true;
				});
			}
		}

		return true;
	});

	for (Volume v : volume_disable)
	{
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		if (t->pere->pere != nullptr)
		{
			bool result = true;
			t->pere->pere->for_each_child([&](tree_volume* c) -> bool {
				if (c->type != CURRENT)
				{
					result = false;
				}
				return result;
			});
			if (result)
			{
				t->pere->pere->type = COARSE;
			}
		}
	}
#endif
}

void XPBD_SPH_Multiresolution::compute_error_point(MAP& m, const Vec3& p, double quotat)
{
	std::vector<Volume> volume_activate;
	std::vector<Volume> volume_disable;
	std::forward_list<tree_volume*> list_volume_coarse;
	std::forward_list<tree_volume*> list_volume_fine;
	uint32 nb_fine = 0, nb_coarse = 0;
	static uint32 clock_error = 0;
	clock_error++;

	foreach_cell(m, [&](Volume v) -> bool {
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		if (t->type == CURRENT && t->fils != nullptr)
		{
			list_volume_fine.push_front(t);
			double error = (geometry::centroid<Vec3>(m, Volume(t->volume_dart), pos_.get()) - p).squaredNorm();
			value<double>(m, error_volume_, Volume(t->volume_dart)) = error;
			t->error = error;
			nb_fine++;
		}
		if (!t->is_topo && t->pere != nullptr && t->pere->type == COARSE)
		{
			if (t->pere->clock == clock_error)
				return true;
			t->pere->clock = clock_error;
			list_volume_coarse.push_front(t);
			t->pere->error = 0;
			t->pere->for_each_child([&](tree_volume* c) -> bool {
				double error = (geometry::centroid<Vec3>(m, Volume(c->volume_dart), pos_.get()) - p).squaredNorm();
				t->pere->error += error;
				return true;
			});
			t->pere->error /= 8;
			nb_coarse++;
		}
		return true;
	});
	list_volume_fine.sort([&](tree_volume* t1, tree_volume* t2) { return t1->error < t2->error; });
	list_volume_coarse.sort([&](tree_volume* t1, tree_volume* t2) { return t1->pere->error > t2->pere->error; });

	auto it = list_volume_coarse.begin();
	double e_max = 0;
	if (!list_volume_coarse.empty())
	{
		e_max = list_volume_coarse.front()->pere->error;
	}

	while (1)
	{
		if (list_volume_fine.empty())
			break;

		tree_volume* t = list_volume_fine.front();
		list_volume_fine.pop_front();
		double e = t->error;
		if (e > e_max && nb_volume_current > quotat * nb_volume_init)
		{
			break;
		}
		if (e < e_max)
		{
			it++;
			e_max = 0;
			if (it != list_volume_coarse.end())
			{
				e_max = (*it)->pere->error;
			}
		}

		nb_volume_current += 7;
		volume_activate.push_back(Volume(t->volume_dart));
		if (t->pere && t->pere->type != ROOT)
		{
			t->pere->type = NONE;
		}
		t->type = COARSE;
		t->for_each_child([&](tree_volume* c) -> bool {
			c->type = CURRENT;
			return true;
		});
	}

	while (nb_volume_current > quotat * nb_volume_init)
	{
		if (list_volume_coarse.empty())
			break;
		tree_volume* t = list_volume_coarse.front();
		list_volume_coarse.pop_front();
		if (t->pere->type != COARSE)
			continue;
		nb_volume_current -= 7;
		volume_disable.push_back(Volume(t->pere->volume_dart));
		t->pere->type = CURRENT;
		t->pere->for_each_child([&](tree_volume* c) -> bool {
			c->type = NONE;
			return true;
		});
		if (t->pere->pere != nullptr)
		{
			bool result = true;
			t->pere->pere->for_each_child([&](tree_volume* c) -> bool {
				if (c->type != CURRENT)
				{
					result = false;
				}
				return result;
			});
			if (result)
			{
				t->pere->pere->type = COARSE;
			}
		}
	}
	activate_remove_volume(m, volume_activate, volume_disable);
}

#define SHOW_PERFORMANCE_LOG 1
void XPBD_SPH_Multiresolution::solver(MAP& m, MAP* geom, double timestep, bool allow_modif_topo)
{
	std::clock_t start;
	double duration;
	start = std::clock();
	double h = timestep / NUM_SUBSTEP;
	std::vector<Volume> vec_volume;
	std::vector<Vertex> vec_vertices;
	foreach_cell(m, [&](Vertex v) -> bool {
		vec_vertices.push_back(v);
		return true;
	});
	foreach_cell(m, [&](Volume v) -> bool {
		vec_volume.push_back(v);
		std::vector<Vertex>& vector_inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		vector_inc_vertices.clear();
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			vector_inc_vertices.push_back(w);
			return true;
		});
		return true;
	});
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
	std::cout << "\033[1;32mtime init XPBD : \033[0m" << duration << std::endl;
	std::cout << "\033[1;32mnb DOFs : \033[0m" << vec_vertices.size() << std::endl;
#endif

	start = std::clock();
	for (int i = 0; i < NUM_SUBSTEP; i++)
	{
		// Initialisation sub step
		for (Vertex v : vec_vertices)
		{
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
			{
				continue;
			}
			value<Vec3>(m, pos_prev_, v) = value<Vec3>(m, pos_, v);
			value<Vec3>(m, speed_, v) += h * value<Vec3>(m, f_ext_, v) / value<double>(m, masse_, v);
			value<Vec3>(m, pos_, v) += h * value<Vec3>(m, speed_, v);
		}
		// Constraint solver
		for (Volume v : vec_volume)
		{
			constraint_Neo_Hookean_H(m, v, h);
			constraint_Neo_Hookean_D(m, v, h);
			constraint_Zero_Energy(m, v, h);
		}
		// Speed Solver
		for (Vertex v : vec_vertices)
		{
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
			{
				continue;
			}
			Vec3 new_v = (value<Vec3>(m, pos_, v) - value<Vec3>(m, pos_prev_, v)) / h;
			for (int i = 0; i < 3; i++)
				if (fabs(new_v[i]) < EPS)
					new_v[i] = 0;
			// value<Vec3>(m, speed_, v) = (1 - (0.005 * h)) * new_v;
			value<Vec3>(m, speed_, v) = new_v;
		}
		// Damping
		foreach_cell(m, [&](Volume v) -> bool {
			applyDamping(m, v, 0.1, timestep);
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
	std::cout << "\033[1;36mtime resolve XPBD : \033[0m" << duration << std::endl;
#endif
	start = std::clock();
	if (geom != nullptr)
	{
		for (Volume v : vec_volume)
		{
			solve_surface(m, *geom, v);
		}
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
	std::cout << "\033[1;36mtime update surface XPBD : \033[0m" << duration << std::endl;
#endif
	std::vector<Volume> vol_activate;
	std::vector<Volume> vol_disable;

	if (allow_modif_topo)
	{
		start = std::clock();
		compute_error(m, vol_activate, vol_disable);
		activate_remove_volume(m, vol_activate, vol_disable);

		// activate_volume(m, vol_activate);
		// remove_volume(m, vol_disable);

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
		std::cout << "\033[1;36mtime update topo XPBD : \033[0m" << duration << std::endl;
#endif
	}
}

} // namespace simulation
} // namespace cgogn
