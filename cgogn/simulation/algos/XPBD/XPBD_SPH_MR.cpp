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
}

void XPBD_SPH_Multiresolution::activate_volume(MAP& m, std::vector<Volume>& list_Volumes)
{
}

void XPBD_SPH_Multiresolution::remove_volume(MAP& m, std::vector<Volume>& list_Volumes)
{
}

void XPBD_SPH_Multiresolution::update_topo(MAP& m)
{
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
	compute_deformation_gradient_particule(&p);

	Mat3d F = p.deformation_gradient_;

	p.CH = F.determinant() - (1. + LAME_MU / LAME_LAMBDA);

	// Compute Volume
	// double Ve = geometry::volume(m, v, pos_.get());
	double Ve = fabs(p.deformation_gradient_.determinant() * p.initial_volume_);

	// Compute alpha_H = 1/(LAMBDA*V)
	double alpha_h = 1.0 / (LAME_LAMBDA * Ve);

	// Compute denum
	double denum = 0;
	Mat3d dCHdF;
	dCHdF.col(0) = F.col(1).cross(F.col(2));
	dCHdF.col(1) = F.col(2).cross(F.col(0));
	dCHdF.col(2) = F.col(0).cross(F.col(1));

	denum = 1. / p.masse_ * (-dCHdF * p.dFdX).squaredNorm();

	Vec3 GC;

	for (Particule_SPH_MR* w : n)
	{
		double m_j = w->masse_;
		double V_j = w->initial_volume_;
		Vec3 xi0 = p.initial_position_;
		Vec3 xj0 = w->initial_position_;

		denum += 1. / m_j * (dCHdF * V_j * Corrected_gradient(p, *w, p.h_)).squaredNorm();
	}

	denum += alpha_h / h / h;

	// Compute lambda
	double lambda = -p.CH / denum;

	if (!p.is_fixed)
	{
		Vec3 deltaX = lambda * (1. / p.masse_) * (-dCHdF * p.dFdX);
		p.current_position_ += deltaX;
	}

	for (Particule_SPH_MR* w : n)
	{
		if (w->is_fixed)
			continue;
		Vec3 deltaX = lambda * (1. / w->masse_) * (dCHdF * w->initial_volume_ * Corrected_gradient(p, *w, h_smooth));
		w->current_position_ += deltaX;
	}
}
void XPBD_SPH_Multiresolution::constraint_Neo_Hookean_D(MAP& m, Particule_SPH_MR& p, double h)
{
	double h_smooth = p.h_;
	std::forward_list<Particule_SPH_MR*>& n = p.neighborhood_;
	compute_deformation_gradient_particule(&p);

	Mat3d F = p.deformation_gradient_;

	p.CD = sqrt((F.transpose() * F).trace());

	// Compute Volume
	// double Ve = geometry::volume(m, v, pos_.get());
	double Ve = fabs(p.deformation_gradient_.determinant() * p.initial_volume_);

	// Compute alpha_H = 1/(LAMBDA*V)
	double alpha_d = 1.0 / (LAME_MU * Ve);

	// Compute denum
	double denum = 0;
	Mat3d dCDdF = F / p.CD;

	denum = 1. / p.masse_ * (-dCDdF * p.dFdX).squaredNorm();

	for (Particule_SPH_MR* w : n)
	{
		double m_j = w->masse_;
		double V_j = w->initial_volume_;
		Vec3 xi0 = p.initial_position_;
		Vec3 xj0 = w->initial_position_;

		denum += 1. / m_j * (dCDdF * V_j * Corrected_gradient(p, *w, p.h_)).squaredNorm();
	}

	denum += alpha_d / h / h;

	// Compute lambda
	double lambda = -p.CH / denum;

	if (!p.is_fixed)
	{
		Vec3 deltaX = lambda * (1. / p.masse_) * (-dCDdF * p.dFdX);
		p.current_position_ += deltaX;
	}

	for (Particule_SPH_MR* w : n)
	{
		if (w->is_fixed)
			continue;
		Vec3 deltaX = lambda * (1. / w->masse_) * (dCDdF * w->initial_volume_ * Corrected_gradient(p, *w, h_smooth));
		w->current_position_ += deltaX;
	}
}
void XPBD_SPH_Multiresolution::constraint_Zero_Energy(MAP& m, Particule_SPH_MR& p, double)
{
}
void XPBD_SPH_Multiresolution::applyDamping(MAP& m, Volume v, double damping_coeff, double time_step)
{
}

void XPBD_SPH_Multiresolution::solve_surface(MAP& m, MAP& geom, Volume v)
{
}

void XPBD_SPH_Multiresolution::compute_error(MAP& m, std::vector<Volume>& volume_activate,
											 std::vector<Volume>& volume_disable)
{
}

void XPBD_SPH_Multiresolution::compute_error_point(MAP& m, const Vec3& p, double quotat)
{
}

#define SHOW_PERFORMANCE_LOG 1
void XPBD_SPH_Multiresolution::solver(MAP& m, MAP* geom, double timestep, bool allow_modif_topo)
{
	std::clock_t start;
	double duration;
	start = std::clock();
	double h = timestep / NUM_SUBSTEP;
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
	std::cout << "\033[1;32mtime init XPBD SPH : \033[0m" << duration << std::endl;
	std::cout << "\033[1;32mnb DOFs : \033[0m" << particules_.size() << std::endl;
#endif

	start = std::clock();
	for (int i = 0; i < NUM_SUBSTEP; i++)
	{
		// Initialisation sub step
		for (Particule_SPH_MR* p : particules_)
		{
			if (p->is_fixed)
			{
				continue;
			}
			p->previous_position_ = p->current_position_;
			p->speed_ += h * p->force_ / p->masse_;
			p->current_position_ += h * p->speed_;
		}
		// Constraint solver
		for (Particule_SPH_MR* p : particules_)
		{
			constraint_Neo_Hookean_H(m, *p, h);
			constraint_Neo_Hookean_D(m, *p, h);
		}
		// Speed Solver
		for (Particule_SPH_MR* p : particules_)
		{
			if (p->is_fixed)
			{
				continue;
			}
			Vec3 new_v = (p->current_position_ - p->previous_position_) / h;
			for (int i = 0; i < 3; i++)
				if (fabs(new_v[i]) < EPS)
					new_v[i] = 0;
			p->speed_ = new_v;
		}
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
	std::cout << "\033[1;36mtime resolve XPBD SPH : \033[0m" << duration << std::endl;
#endif
	start = std::clock();
	if (geom != nullptr)
	{
		/*for (Volume v : vec_volume)
		{
			solve_surface(m, *geom, v);
		}*/
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
	std::cout << "\033[1;36mtime update surface XPBD SPH : \033[0m" << duration << std::endl;
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
		std::cout << "\033[1;36mtime update topo XPBD SPH : \033[0m" << duration << std::endl;
#endif
	}
}

} // namespace simulation
} // namespace cgogn
