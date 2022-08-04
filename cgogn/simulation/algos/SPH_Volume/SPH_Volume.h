#ifndef CGOGN_SIMULATION_SPH_VOLUME_SPH_VOLUME_H_
#define CGOGN_SIMULATION_SPH_VOLUME_SPH_VOLUME_H_
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>

#define NORMALIZE_TERM (21 / (2 * M_PI))
#define POISSON_RATIO 0.33
#define YOUNG_MODULUS 5e4
#define LAME_MU 2e4
#define LAME_LAMBDA (LAME_MU / 3.0)

#define SHEAR_MODULUS (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define BULK_MODULUS (YOUNG_MODULUS / (3 * (1 - 2 * POISSON_RATIO)))
#define DENSITY_SPH 10

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class SPH_volume_constraint_solver : public Simulation_constraint<MAP>
{
	using Self = SPH_volume_constraint_solver;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Volume = typename mesh_traits<MAP>::Volume;

public:
	static inline int nb_solver = 0;
	int id;
	std::shared_ptr<Attribute<double>> initial_volume_;
	std::shared_ptr<Attribute<Vec3>> initial_centroid_volume_;
	std::shared_ptr<Attribute<Vec3>> centroid_volume_;
	std::shared_ptr<Attribute<std::vector<Volume>>> neighborhood_volume_;
	std::shared_ptr<Attribute<Mat3d>> corrected_matrix_volume_;
	std::shared_ptr<Attribute<Mat3d>> rotation_volume_;
	std::shared_ptr<Attribute<double>> h_volume_;
	std::shared_ptr<Attribute<Mat3d>> stress_tensor_volume_;
	std::shared_ptr<Attribute<Vec3>> force_volume_;

	SPH_volume_constraint_solver()
		: id(nb_solver++), initial_volume_(nullptr), initial_centroid_volume_(nullptr), centroid_volume_(nullptr),
		  neighborhood_volume_(nullptr), corrected_matrix_volume_(nullptr), rotation_volume_(nullptr),
		  h_volume_(nullptr), stress_tensor_volume_(nullptr), force_volume_(nullptr)
	{
	}

	Simulation_constraint<MAP>* get_new_ptr()
	{
		return new SPH_volume_constraint_solver<MAP>();
	}

	double Kernel_W(double dist, double h) const
	{
		double q = dist / h;
		if (q > 1)
			return 0;
		double m1 = (1 - q);
		double m2 = (4 * q + 1);
		double h3 = h * h * h;
		double alpha_d = NORMALIZE_TERM * (1 / h3);

		return alpha_d * m1 * m1 * m1 * m1 * m2;
	}

	double dwdq(double dist, double h) const
	{
		double q = dist / h;
		if (q > 2 || dist <= 1e-12)
		{
			return 0;
		}
		double h3 = h * h * h;
		double alpha_d = NORMALIZE_TERM * (1 / h3);
		double tmp = 1 - 0.5 * q;
		return alpha_d * (-5.0 * q * tmp * tmp * tmp);
	}

	Vec3 gradient(Vec3 xij, double dist, double h) const
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

	Vec3 gradient2(Vec3 xij, double dist, double h) const
	{
		if (dist < 1e-12)
			return Vec3(0, 0, 0);
		double tmp = dwdq(dist, h) / (h * dist);

		return xij * tmp;
	}

	Mat3d Corrected_matrix(const MAP& m, Volume v, double h) const
	{
		Mat3d Li = Mat3d::Zero();
		std::vector<Volume>& n = value<std::vector<Volume>>(m, neighborhood_volume_.get(), v);
		Vec3 xi = value<Vec3>(m, initial_centroid_volume_.get(), v);
		for (Volume w : n)
		{
			double init_vol = value<double>(m, initial_volume_.get(), v);
			Vec3 xj = value<Vec3>(m, initial_centroid_volume_.get(), w);
			Vec3 xij = xi - xj;
			Vec3 xji = -xij;
			Vec3 grad = gradient(xij, xij.norm(), h);
			Li += init_vol * grad * xji.transpose();
		}
		return Li.inverse();
	}

	void compute_corrected_matrix(const MAP& m) const
	{
		parallel_foreach_cell(m, [&](Volume v) -> bool {
			value<Mat3d>(m, corrected_matrix_volume_.get(), v) =
				Corrected_matrix(m, v, value<double>(m, h_volume_.get(), v));
			return true;
		});
	}

	Vec3 Corrected_gradient(const MAP& m, Volume vi, Volume vj, double h) const
	{
		Vec3 xi = value<Vec3>(m, initial_centroid_volume_.get(), vi);
		Vec3 xj = value<Vec3>(m, initial_centroid_volume_.get(), vj);
		Vec3 xij = xi - xj;

		return value<Mat3d>(m, corrected_matrix_volume_.get(), vi) * gradient(xij, xij.norm(), h);
	}

	Mat3d Rotation_extraction(const MAP& m, Volume v, double h) const
	{
		Mat3d R;
		Mat3d F = Mat3d::Zero();

		std::vector<Volume>& n = value<std::vector<Volume>>(m, neighborhood_volume_.get(), v);
		Vec3 xi = value<Vec3>(m, centroid_volume_.get(), v);
		for (Volume w : n)
		{
			double init_vol = value<double>(m, initial_volume_.get(), v);
			Vec3 xj = value<Vec3>(m, centroid_volume_.get(), w);
			Vec3 xij = xi - xj;
			Vec3 W = Corrected_gradient(m, v, w, h);
			// minus because xij = - xji
			F -= init_vol * xij * W.transpose();
		}
		polarDecompositionStable(F, 1.0e-6, R);
		return R;
	}

	void compute_rotated_kernel(const MAP& m)
	{
		cgogn::parallel_foreach_cell(m, [&](Volume v) -> bool {
			double h = value<double>(m, h_volume_.get(), v);
			value<Mat3d>(m, rotation_volume_.get(), v) = Rotation_extraction(m, v, h);
			return true;
		});
	}

	Vec3 Rotated_gradient(const MAP& m, Volume vi, Volume vj, double h) const
	{
		return value<Mat3d>(m, rotation_volume_.get(), vi) * Corrected_gradient(m, vi, vj, h);
	}

	void compute_force_volume(const MAP& m)
	{
		compute_rotated_kernel(m);

		cgogn::parallel_foreach_cell(m, [&](Volume v) -> bool {
			double h = value<double>(m, h_volume_.get(), v);
			Mat3d Ftemp = Mat3d::Identity();
			std::vector<Volume>& n = value<std::vector<Volume>>(m, neighborhood_volume_.get(), v);
			Vec3 xi = value<Vec3>(m, centroid_volume_.get(), v);
			Vec3 xi0 = value<Vec3>(m, initial_centroid_volume_.get(), v);
			Mat3d& Ri = value<Mat3d>(m, rotation_volume_.get(), v);
			for (Volume w : n)
			{
				double Vj0 = value<double>(m, initial_volume_.get(), w);
				Vec3 xj = value<Vec3>(m, centroid_volume_.get(), w);
				Vec3 xji = xj - xi;
				Vec3 xj0 = value<Vec3>(m, initial_centroid_volume_.get(), w);
				Vec3 xji0 = xj0 - xi0;
				Vec3 W = Rotated_gradient(m, v, w, h);
				Ftemp += Vj0 * (xji - Ri * xji0) * W.transpose();
			}
			Mat3d Etemp = 0.5 * (Ftemp + Ftemp.transpose()) - Mat3d::Identity();
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					if (fabs(Etemp(i, j)) < 1e-12)
					{
						Etemp(i, j) = 0;
					}
				}
			}
			Mat3d Pi = 2 * SHEAR_MODULUS * Etemp +
					   (BULK_MODULUS - (2.0 / 3.0) * SHEAR_MODULUS) * Etemp.trace() * Mat3d::Identity();
			value<Mat3d>(m, stress_tensor_volume_.get(), v) = Pi;
			return true;
		});
		cgogn::parallel_foreach_cell(m, [&](Volume v) -> bool {
			double hi = value<double>(m, h_volume_.get(), v);
			Vec3 F = Vec3::Zero();
			std::vector<Volume>& n = value<std::vector<Volume>>(m, neighborhood_volume_.get(), v);
			Mat3d Pi = value<Mat3d>(m, stress_tensor_volume_.get(), v);
			double Vi0 = value<double>(m, initial_volume_.get(), v);
			for (Volume w : n)
			{
				double Vj0 = value<double>(m, initial_volume_.get(), w);
				double hj = value<double>(m, h_volume_.get(), w);
				Mat3d Pj = value<Mat3d>(m, stress_tensor_volume_.get(), w);
				Vec3 Wi = Rotated_gradient(m, v, w, hi);
				Vec3 Wj = Rotated_gradient(m, w, v, hj);
				F += Vj0 * Vi0 * (Pi * Wi - Pj * Wj);
			}
			value<Vec3>(m, force_volume_, v) = F;
			return true;
		});
	}

	void init_solver(MAP& m, const std::shared_ptr<Attribute<Vec3>>& init_pos,
					 const std::shared_ptr<Attribute<double>>& masse)
	{
	}

	void init_solver(MAP& m, Attribute<Vec3>* pos)
	{
		initial_volume_ = add_attribute<double, Volume>(m, "SPH_simulation_constraint_solver_initial_volume_" + id);
		initial_centroid_volume_ =
			add_attribute<Vec3, Volume>(m, "SPH_simulation_constraint_solver_initial_centroid_volume_" + id);
		centroid_volume_ = add_attribute<Vec3, Volume>(m, "SPH_simulation_constraint_solver_centroid_volume_" + id);
		neighborhood_volume_ =
			add_attribute<std::vector<Volume>, Volume>(m, "SPH_simulation_constraint_solver_neighborhood_volume_" + id);
		corrected_matrix_volume_ =
			add_attribute<Mat3d, Volume>(m, "SPH_simulation_constraint_solver_corrected_matrix_volume_" + id);
		rotation_volume_ = add_attribute<Mat3d, Volume>(m, "SPH_simulation_constraint_solver_rotation_volume_" + id);
		h_volume_ = add_attribute<double, Volume>(m, "SPH_simulation_constraint_solver_h_volume_" + id);
		stress_tensor_volume_ =
			add_attribute<Mat3d, Volume>(m, "SPH_simulation_constraint_solver_stress_tensor_volume_" + id);
		force_volume_ = add_attribute<Vec3, Volume>(m, "SPH_simulation_constraint_solver_force_volume_" + id);

		geometry::compute_centroid<Vec3, Volume>(m, pos, initial_centroid_volume_.get());
		geometry::compute_volume(m, pos, initial_volume_.get());

		parallel_foreach_cell(m, [&](Volume v) -> bool {
			CellMarker<MAP, Volume> marker(m);
			std::vector<Volume>& n = value<std::vector<Volume>>(m, neighborhood_volume_.get(), v);
			n.clear();
			double& h = value<double>(m, h_volume_.get(), v);
			h = 0;
			Vec3 centroid_v1 = value<Vec3>(m, initial_centroid_volume_.get(), v);
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				foreach_incident_volume(m, w, [&](Volume v2) -> bool {
					foreach_incident_vertex(m, v2, [&](Vertex w2) -> bool {
						foreach_incident_volume(m, w2, [&](Volume v3) -> bool {
							if (!marker.is_marked(v3))
							{
								Vec3 centroid_v3 = value<Vec3>(m, initial_centroid_volume_.get(), v3);
								h = std::max(h, (centroid_v1 - centroid_v3).norm());
								n.push_back(v3);
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
			h *= 1.2;
			return true;
		});
		compute_corrected_matrix(m);
	}

	void update_topo(const MAP& m, const std::vector<Vertex>&)
	{
	}

	double oneNorm(const Mat3d& A) const
	{
		const double sum1 = fabs(A(0, 0)) + fabs(A(1, 0)) + fabs(A(2, 0));
		const double sum2 = fabs(A(0, 1)) + fabs(A(1, 1)) + fabs(A(2, 1));
		const double sum3 = fabs(A(0, 2)) + fabs(A(1, 2)) + fabs(A(2, 2));
		double maxSum = sum1;
		if (sum2 > maxSum)
			maxSum = sum2;
		if (sum3 > maxSum)
			maxSum = sum3;
		return maxSum;
	}

	double infNorm(const Mat3d& A) const
	{
		const double sum1 = fabs(A(0, 0)) + fabs(A(0, 1)) + fabs(A(0, 2));
		const double sum2 = fabs(A(1, 0)) + fabs(A(1, 1)) + fabs(A(1, 2));
		const double sum3 = fabs(A(2, 0)) + fabs(A(2, 1)) + fabs(A(2, 2));
		double maxSum = sum1;
		if (sum2 > maxSum)
			maxSum = sum2;
		if (sum3 > maxSum)
			maxSum = sum3;
		return maxSum;
	}

	void polarDecompositionStable(const Mat3d& M, double tolerance, Mat3d& R) const
	{
		Mat3d Mt = M.transpose();
		double Mone = oneNorm(M);
		double Minf = infNorm(M);
		double Eone;
		Mat3d MadjTt, Et;
		do
		{
			MadjTt.row(0) = Mt.row(1).cross(Mt.row(2));
			MadjTt.row(1) = Mt.row(2).cross(Mt.row(0));
			MadjTt.row(2) = Mt.row(0).cross(Mt.row(1));

			double det = Mt(0, 0) * MadjTt(0, 0) + Mt(0, 1) * MadjTt(0, 1) + Mt(0, 2) * MadjTt(0, 2);

			if (fabs(det) < 1.0e-12)
			{
				Vec3 len;
				unsigned int index = 0xffffffff;
				for (unsigned int i = 0; i < 3; i++)
				{
					len[i] = MadjTt.row(i).squaredNorm();
					if (len[i] > 1.0e-12)
					{
						// index of valid cross product
						// => is also the index of the vector in Mt that must be exchanged
						index = i;
						break;
					}
				}
				if (index == 0xffffffff)
				{
					R.setIdentity();
					return;
				}
				else
				{
					Mt.row(index) = Mt.row((index + 1) % 3).cross(Mt.row((index + 2) % 3));
					MadjTt.row((index + 1) % 3) = Mt.row((index + 2) % 3).cross(Mt.row(index));
					MadjTt.row((index + 2) % 3) = Mt.row(index).cross(Mt.row((index + 1) % 3));
					Mat3d M2 = Mt.transpose();
					Mone = oneNorm(M2);
					Minf = infNorm(M2);
					det = Mt(0, 0) * MadjTt(0, 0) + Mt(0, 1) * MadjTt(0, 1) + Mt(0, 2) * MadjTt(0, 2);
				}
			}

			const double MadjTone = oneNorm(MadjTt);
			const double MadjTinf = infNorm(MadjTt);

			const double gamma = sqrt(sqrt((MadjTone * MadjTinf) / (Mone * Minf)) / fabs(det));

			const double g1 = gamma * static_cast<double>(0.5);
			const double g2 = static_cast<double>(0.5) / (gamma * det);

			for (unsigned char i = 0; i < 3; i++)
			{
				for (unsigned char j = 0; j < 3; j++)
				{
					Et(i, j) = Mt(i, j);
					Mt(i, j) = g1 * Mt(i, j) + g2 * MadjTt(i, j);
					Et(i, j) -= Mt(i, j);
				}
			}

			Eone = oneNorm(Et);

			Mone = oneNorm(Mt);
			Minf = infNorm(Mt);
		} while (Eone > Mone * tolerance);

		// Q = Mt^T
		R = Mt.transpose();
	}

	void solve_constraint(const MAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces, double) override
	{
		geometry::compute_centroid<Vec3, Volume>(m, pos, centroid_volume_.get());
		compute_force_volume(m);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			int nb_volume = 0;
			Vec3 f = Vec3::Zero();
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				f += value<Vec3>(m, force_volume_, w);
				nb_volume++;
				return true;
			});
			value<Vec3>(m, result_forces, v) += f / double(nb_volume);
			return true;
		});
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_SPH_VOLUME_SPH_VOLUME_H_
