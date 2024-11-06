#include "XPBD_multiresolution.h"
#include <Eigen/SVD>
#include <cgogn/geometry/algos/centroid.h>
#include <forward_list>

namespace cgogn
{
namespace simulation
{
void XPBD_Multiresolution::init_solver(MAP& m, std::shared_ptr<Attribute<Vec3>> pos)
{
	nb_volume_init = 0;
	pos_ = pos;
	init_pos_ = add_attribute<Vec3, Vertex>(m, "XPBD_Init_pos");
	init_cm_ = add_attribute<Vec3, Volume>(m, "XPBD_init_cm");
	centroid_ = add_attribute<Vec3, Volume>(m, "XPBD_centroid");
	distance_plan_ = add_attribute<double, Volume>(m, "XPBD_distance_plan");

	masse_ = add_attribute<double, Vertex>(m, "XPBD_masse");
	s_ = add_attribute<double, Volume>(m, "XPBD_s_normalize");

	init_volume_ = add_attribute<double, Volume>(m, "XPBD_init_volume");
	Det_F_Volume_ = add_attribute<double, Volume>(m, "XPBD_Det_F_volume");
	F_ = add_attribute<Mat3d, Volume>(m, "XPBD_F_volume");
	inv_Q_ = add_attribute<Mat3d, Volume>(m, "XPBD_inv_Q");
	inc_vertices_ = add_attribute<std::vector<Vertex>, Volume>(m, "XPBD_inc_vertices_vector_");

	pos_prev_ = add_attribute<Vec3, Vertex>(m, "XPBD_pos_prev_");
	speed_ = add_attribute<Vec3, Vertex>(m, "XPBD_speed_");

	f_ext_ = add_attribute<Vec3, Vertex>(m, "XPBD_f_ext");

	Grad_C_i_ = add_attribute<Vec3, Vertex>(m, "XPBD_Grad_C_i");
	Grad_C2_i_ = add_attribute<Vec3, Vertex>(m, "XPBD_Grad_C2_i");
	error_volume_ = add_attribute<double, Volume>(m, "XPBD_error_volume");

	using MR_Base = typename MAP::Inherit;
	MR_Base tmp(m);
	tmp.current_level_ = 0;
	std::function<void(tree_volume*, MR_Base&)> progress_tree;

	int nb_node = 0;

	progress_tree = [&](tree_volume* p, MR_Base& cph) -> void {
		if (cph.current_level_ == cph.maximum_level_)
			return;
		foreach_incident_vertex(cph, Volume(p->volume_dart), [&](Vertex v) -> bool {
			cph.current_level_++;
			if (cph.current_level_ != cph.volume_level(v.dart))
				return false;

			tree_volume* t = new tree_volume();
			t->id = nb_node++;
			t->volume_dart = cph.volume_oldest_dart(v.dart);
			t->frere = p->fils;
			p->fils = t;
			t->pere = p;
			t->type = NONE;
			value<tree_volume*>(cph, hierarchy_node_, Volume(v.dart)) = t;
			if (cph.volume_is_subdivided(v.dart))
				progress_tree(t, cph);

			cph.current_level_--;

			return true;
		});
	};

	hierarchy_node_ = get_attribute<tree_volume*, Volume>(m, "Solver_multiresolution_hierarchy_node");
	if (hierarchy_node_ == nullptr)
		hierarchy_node_ = add_attribute<tree_volume*, Volume>(m, "Solver_multiresolution_hierarchy_node");

	hierarchy_ = new tree_volume();
	hierarchy_->type = ROOT;
	hierarchy_->is_topo = true;
	hierarchy_->id = nb_node++;

	foreach_cell(tmp, [&](typename MR_Base::Volume v) -> bool {
		MR_Base tmp2(tmp);
		tree_volume* t = new tree_volume();
		t->id = nb_node++;
		t->frere = hierarchy_->fils;
		hierarchy_->fils = t;
		t->pere = hierarchy_;
		t->volume_dart = tmp.volume_oldest_dart(v.dart);
		t->type = TOPOLOGY;
		t->is_topo = true;
		value<tree_volume*>(tmp, hierarchy_node_, v) = t;
		progress_tree(t, tmp2);
		return true;
	});

	uint32 cur = m.current_level_;
	m.current_level_ = m.maximum_level_;
	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		value<Vec3>(m, init_pos_, v) = value<Vec3>(m, pos_, v);
		value<double>(m, masse_, v) = 0;
		value<Vec3>(m, speed_, v) = Vec3(0, 0, 0);
		value<Vec3>(m, f_ext_, v) = Vec3(0, 0, 0);
		return true;
	});
	geometry::compute_volume(m, pos_.get(), init_volume_.get());
	m.current_level_ = cur;

	tmp.current_level_ = m.maximum_level_ - 1;
	for (int i = m.maximum_level_ - 1; i >= 0; i--)
	{
		tmp.current_level_ = i;
		foreach_cell(tmp, [&](Volume v) -> bool {
			tree_volume* t = value<tree_volume*>(tmp, hierarchy_node_, v);
			if (t->fils != nullptr)
			{
				double vol = 0;
				tmp.current_level_++;
				t->for_each_child([&](tree_volume* c) -> bool {
					vol += value<double>(tmp, init_volume_, Volume(c->volume_dart));
					return true;
				});
				tmp.current_level_--;
				value<double>(tmp, init_volume_, v) = vol;
			}
			return true;
		});
	}

	foreach_cell(m, [&](Volume v) -> bool {
		nb_volume_init++;
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
		return true;
	});
	nb_volume_current = nb_volume_init;
	nb_dof_init = get_nb_dof(m);
	nb_dof_current = nb_dof_init;

	parallel_foreach_cell(m, [&](Volume v) -> bool {
		double masse = 0;
		Vec3 cm = Vec3(0, 0, 0);
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			masse += value<double>(m, masse_, w);
			cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
			return true;
		});
		value<Vec3>(m, init_cm_, v) = cm / masse;
		return true;
	});

	parallel_foreach_cell(m, [&](Volume v) -> bool {
		Mat3d Q = Mat3d::Zero();
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			double masse = value<double>(m, masse_, w);
			Q.noalias() += masse * init_r_i * init_r_i.transpose();
			return true;
		});
		double s = 1.0 / Q.sum();
		value<double>(m, s_, v) = s;
		// s = 1.;
		Q = s * Q;
		value<Mat3d>(m, inv_Q_, v) = Q.eval().inverse();
		return true;
	});
	foreach_cell(m, [&](Volume v) -> bool {
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		t->type = CURRENT;
		return true;
	});
	foreach_cell(m, [&](Volume v) -> bool {
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
		if (t->pere != nullptr && t->pere->type != COARSE && t->pere->type != ROOT)
		{
			t->pere->type = COARSE;
		}
		return true;
	});
}

void XPBD_Multiresolution::activate_remove_volume(MAP& m, std::vector<Volume>& list_activate,
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

void XPBD_Multiresolution::activate_volume_tree(MAP& m, std::vector<Volume>& list_Volumes)
{
	for (Volume v : list_Volumes)
	{
		tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);

		nb_volume_current += 7;
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

	activate_volume(m, list_Volumes);
}

void XPBD_Multiresolution::activate_volume(MAP& m, std::vector<Volume>& list_Volumes)
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
			c->F_ = t->F_;

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
		std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
		inc_vertices.clear();
		t = t->pere;
		value<Mat3d>(m, F_, v) = t->F_;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			value<Vec3>(m, speed_, w) = t->v_cm_;
			Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - t->init_cm_;
			value<Vec3>(m, pos_, w) = t->cm_ + t->F_ * init_r_i;
			inc_vertices.push_back(w);
			return true;
		});
	}
}

void XPBD_Multiresolution::remove_volume(MAP& m, std::vector<Volume>& list_Volumes)
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

void XPBD_Multiresolution::update_topo(MAP& m)
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

void XPBD_Multiresolution::constraint_Neo_Hookean_H(MAP& m, Volume v, double h)
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
	value<Vec3>(m, centroid_, v) = cm;

	// Compute P
	Mat3d P = Mat3d::Zero();
	for (Vertex w : inc_vertices)
	{
		double m_i = value<double>(m, masse_, w);
		Vec3 r_i = value<Vec3>(m, pos_, w) - cm;
		Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
		P.noalias() += m_i * r_i * init_r_i.transpose();
	}

	Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
	// compute F = P*Q^-1
	double s = value<double>(m, s_, v);
	// P = (s * P).eval();
	Mat3d F = (s * P * inv_Q).eval();

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

	value<Mat3d>(m, F_, v) = F;

	// Compute C  = det(F) - (1+MU/LAMBDA)
	double det_F = F.determinant();
	// HERE - OR + ?
	double C = det_F - (1. + LAME_MU / LAME_LAMBDA);

	// Compute Volume
	// double Ve = geometry::volume(m, v, pos_.get());
	double Ve = fabs(det_F * value<double>(m, init_volume_, v));

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
void XPBD_Multiresolution::constraint_Neo_Hookean_D(MAP& m, Volume v, double h)
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
void XPBD_Multiresolution::constraint_Zero_Energy(MAP& m, Volume v, double)
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
void XPBD_Multiresolution::applyDamping(MAP& m, Volume v, double damping_coeff, double time_step)
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

void XPBD_Multiresolution::solve_surface(MAP& m, MAP& geom, Volume v)
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
void XPBD_Multiresolution::compute_error(MAP& m, std::vector<Volume>& volume_activate,
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

void XPBD_Multiresolution::compute_error_point(MAP& m, const Vec3& p, double quotat)
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
void XPBD_Multiresolution::solver(MAP& m, MAP* geom, double timestep, bool allow_modif_topo)
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
	std::cout << "\033[1;32mnb DOFs topo : \033[0m" << get_nb_dof(*m.topology_) << std::endl;
#endif
	srand(1547989);
	std::random_shuffle(vec_volume.begin(), vec_volume.end());
	start = std::clock();
	for (int i = 0; i < NUM_SUBSTEP; i++)
	{
		// Initialisation sub step
		for (Vertex v : vec_vertices)
		{
			value<Vec3>(m, pos_prev_, v) = value<Vec3>(m, pos_, v);
			if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
			{
				continue;
			}

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
		for (Volume v : vec_volume)
		{
			applyDamping(m, v, 0.1, timestep);
		}
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

template <typename FUNC>
void XPBD_Multiresolution::apply_cut(MAP& m, Vec3 dir_plan, double w, const FUNC& callback_vertices)
{

	geometry::compute_centroid<Vec3, Volume>(m, pos_.get(), centroid_.get());

	parallel_foreach_cell(m, [&](Volume v) -> bool {
		value<double>(m, this->distance_plan_.get(), v) = dir_plan.dot(value<Vec3>(m, this->centroid_.get(), v));
		return true;
	});
	CellMarker<EMR_Map3_Adaptative, Face> face_marker(m);
	std::vector<Face> face_vect;
	foreach_cell(m, [&](Face f) -> bool {
		if (is_incident_to_boundary(m, f))
		{
			return true;
		}
		double v1 = value<double>(m, this->distance_plan_.get(), Volume(f.dart)) - w;
		double v2 = value<double>(m, this->distance_plan_.get(), Volume(phi3(m, f.dart))) - w;
		if (v1 * v2 < 0)
		{
			face_vect.push_back(f);
		}
		face_marker.mark(f);
		return true;
	});
	CellMarker<EMR_Map3_Adaptative, Volume> vol_marker(m);
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
			if (!vol_marker.is_marked(Volume(phi3(m, f.dart))))
			{
				vol_marker.mark(Volume(phi3(m, f.dart)));
				volume_vect.push_back(Volume(phi3(m, f.dart)));
			}
		}
		for (auto v : volume_vect)
		{
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				vect_new_volume.push_back(Volume(w.dart));
				return true;
			});
			m.activate_volume_subdivision(v);
		}
		geometry::compute_centroid<Vec3, Volume>(m, pos_.get(), centroid_.get());
		parallel_foreach_cell(m, [&](Volume v) -> bool {
			value<double>(m, this->distance_plan_.get(), v) = dir_plan.dot(value<Vec3>(m, this->centroid_.get(), v));
			return true;
		});
		face_vect.clear();

		foreach_cell(m, [&](Face f) -> bool {
			if (is_incident_to_boundary(m, f) || face_marker.is_marked(f))
			{
				return true;
			}
			double v1 = value<double>(m, this->distance_plan_.get(), Volume(f.dart)) - w;
			double v2 = value<double>(m, this->distance_plan_.get(), Volume(phi3(m, f.dart))) - w;
			if (v1 * v2 < 0)
			{
				face_vect.push_back(f);
			}
			face_marker.mark(f);
			return true;
		});
	}

	face_vect.clear();
	foreach_cell(m, [&](Face f) -> bool {
		if (is_incident_to_boundary(m, f))
		{
			return true;
		}

		Dart y = m.face_youngest_dart(f.dart);
		double v1 = value<double>(m, this->distance_plan_.get(), Volume(y)) - w;
		double v2 = value<double>(m, this->distance_plan_.get(), Volume(phi3(m, y))) - w;
		if (v1 * v2 < 0)
		{

			face_vect.push_back(f);
		}
		return true;
	});

	// unsew faces
	for (auto f : face_vect)
	{
		unsew_volume(m, f, callback_vertices, true);
	}
}
} // namespace simulation
} // namespace cgogn
