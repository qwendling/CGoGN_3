#ifndef CGOGN_SIMULATION_SPH_PEER_2018_SPH_PEER_2018_H_
#define CGOGN_SIMULATION_SPH_PEER_2018_SPH_PEER_2018_H_
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
class SPH_constraint_solver : public Simulation_constraint<MAP>
{
	using Self = SPH_constraint_solver;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Volume = typename mesh_traits<MAP>::Volume;
	using Quaternion = Eigen::Quaternion<double>;
	using AngleAxisd = Eigen::AngleAxis<double>;

public:
	static inline int nb_solver = 0;
	int id;
	std::shared_ptr<Attribute<double>> initial_vol_;
	std::shared_ptr<Attribute<Vec3>> initial_pos_;
	Attribute<Vec3>* pos_;
	std::shared_ptr<Attribute<std::vector<Vertex>>> neighborhood_;
	std::shared_ptr<Attribute<Mat3d>> corrected_matrix_;
	std::shared_ptr<Attribute<Mat3d>> rotation_;
	std::shared_ptr<Attribute<double>> h_;
	std::shared_ptr<Attribute<Mat3d>> stress_tensor_;
	Attribute<Vec3>* force_;

	SPH_constraint_solver()
		: id(nb_solver++), initial_vol_(nullptr), initial_pos_(nullptr), pos_(nullptr), neighborhood_(nullptr),
		  corrected_matrix_(nullptr), rotation_(nullptr), h_(nullptr), stress_tensor_(nullptr), force_(nullptr)
	{
	}

	Simulation_constraint<MAP>* get_new_ptr()
	{
		return new SPH_constraint_solver<MAP>();
	}

	double Kernel_W(double dist, double h) const
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

	Mat3d Corrected_matrix(const MAP& m, Vertex v, double h) const
	{
		Mat3d Li = Mat3d::Zero();
		std::vector<Vertex>& n = value<std::vector<Vertex>>(m, neighborhood_.get(), v);
		Vec3 xi0 = value<Vec3>(m, initial_pos_.get(), v);
		for (Vertex w : n)
		{
			double init_vol = value<double>(m, initial_vol_.get(), w);
			Vec3 xj0 = value<Vec3>(m, initial_pos_.get(), w);
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

	void compute_corrected_matrix(const MAP& m) const
	{
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Mat3d>(m, corrected_matrix_.get(), v) = Corrected_matrix(m, v, value<double>(m, h_.get(), v));
			return true;
		});
	}

	Vec3 Corrected_gradient(const MAP& m, Vertex vi, Vertex vj, double h) const
	{
		Vec3 xi = value<Vec3>(m, initial_pos_.get(), vi);
		Vec3 xj = value<Vec3>(m, initial_pos_.get(), vj);
		Vec3 xij = xi - xj;

		return value<Mat3d>(m, corrected_matrix_.get(), vi) * gradient(xij, xij.norm(), h);
	}

	Mat3d Rotation_extraction(const MAP& m, Vertex v, double h) const
	{
		Mat3d R;
		Mat3d F = Mat3d::Zero();

		std::vector<Vertex>& n = value<std::vector<Vertex>>(m, neighborhood_.get(), v);
		Vec3 xi = value<Vec3>(m, pos_, v);
		for (Vertex w : n)
		{
			double init_vol = value<double>(m, initial_vol_.get(), w);
			Vec3 xj = value<Vec3>(m, pos_, w);
			Vec3 xji = xj - xi;
			Vec3 W = Corrected_gradient(m, v, w, h);
			F += init_vol * xji * W.transpose();
		}
		polarDecompositionStable(F, 1.0e-6, R);

		/*Quaternion q(value<Mat3d>(m, rotation_.get(), v));
		rotationextraction(F, q, 10);
		R = q.matrix();*/

		return R;
	}

	void compute_rotated_kernel(const MAP& m)
	{
		cgogn::parallel_foreach_cell(m, [&](Vertex v) -> bool {
			double h = value<double>(m, h_.get(), v);
			value<Mat3d>(m, rotation_.get(), v) = Rotation_extraction(m, v, h);
			return true;
		});
	}

	Vec3 Rotated_gradient(const MAP& m, Vertex vi, Vertex vj, double h) const
	{
		return value<Mat3d>(m, rotation_.get(), vi) * Corrected_gradient(m, vi, vj, h);
	}

	void compute_force_vertex(const MAP& m)
	{
		compute_rotated_kernel(m);

		cgogn::parallel_foreach_cell(m, [&](Vertex v) -> bool {
			double h = value<double>(m, h_.get(), v);
			Mat3d Ftemp = Mat3d::Identity();
			std::vector<Vertex>& n = value<std::vector<Vertex>>(m, neighborhood_.get(), v);
			Vec3 xi = value<Vec3>(m, pos_, v);
			Vec3 xi0 = value<Vec3>(m, initial_pos_.get(), v);
			Mat3d& Ri = value<Mat3d>(m, rotation_.get(), v);
			for (Vertex w : n)
			{
				double Vj0 = value<double>(m, initial_vol_.get(), w);
				Vec3 xj = value<Vec3>(m, pos_, w);
				Vec3 xji = xj - xi;
				Vec3 xj0 = value<Vec3>(m, initial_pos_.get(), w);
				Vec3 xji0 = xj0 - xi0;
				Vec3 W = Rotated_gradient(m, v, w, h);
				Ftemp += Vj0 * (xji - Ri * xji0) * W.transpose();
			}
			/*Mat3d Ftemp = Mat3d::Zero();
			std::vector<Vertex>& n = value<std::vector<Vertex>>(m, neighborhood_.get(), v);
			Vec3 xi = value<Vec3>(m, pos_, v);
			for (Vertex w : n)
			{
				double Vj0 = value<double>(m, initial_vol_.get(), w);
				Vec3 xj = value<Vec3>(m, pos_, w);
				Vec3 xji = xj - xi;
				Vec3 W = Rotated_gradient(m, v, w, h);
				Ftemp += Vj0 * xji * W.transpose();
			}*/
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

			value<Mat3d>(m, stress_tensor_.get(), v) = Pi;
			return true;
		});
		cgogn::parallel_foreach_cell(m, [&](Vertex v) -> bool {
			double hi = value<double>(m, h_.get(), v);
			Vec3 F = Vec3::Zero();
			std::vector<Vertex>& n = value<std::vector<Vertex>>(m, neighborhood_.get(), v);
			Mat3d Pi = value<Mat3d>(m, stress_tensor_.get(), v);
			double Vi0 = value<double>(m, initial_vol_.get(), v);
			for (Vertex w : n)
			{
				double Vj0 = value<double>(m, initial_vol_.get(), w);
				double hj = value<double>(m, h_.get(), w);
				Mat3d Pj = value<Mat3d>(m, stress_tensor_.get(), w);
				Vec3 Wi = Rotated_gradient(m, v, w, hi);
				Vec3 Wj = Rotated_gradient(m, w, v, hj);
				F += Vj0 * Vi0 * (Pi * Wi - Pj * Wj);
			}
			value<Vec3>(m, force_, v) += F;
			return true;
		});
	}

	void init_solver(MAP& m, const std::shared_ptr<Attribute<Vec3>>& init_pos,
					 const std::shared_ptr<Attribute<double>>& masse)
	{
	}

	void init_solver(MAP& m, Attribute<Vec3>* pos)
	{
		initial_vol_ = add_attribute<double, Vertex>(m, "SPH_simulation_constraint_solver_initial_vol_" + id);
		initial_pos_ = add_attribute<Vec3, Vertex>(m, "SPH_simulation_constraint_solver_initial_pos_" + id);
		neighborhood_ =
			add_attribute<std::vector<Vertex>, Vertex>(m, "SPH_simulation_constraint_solver_neighborhood_" + id);
		corrected_matrix_ = add_attribute<Mat3d, Vertex>(m, "SPH_simulation_constraint_solver_corrected_matrix_" + id);
		rotation_ = add_attribute<Mat3d, Vertex>(m, "SPH_simulation_constraint_solver_rotation_" + id);
		h_ = add_attribute<double, Vertex>(m, "SPH_simulation_constraint_solver_h_" + id);
		stress_tensor_ = add_attribute<Mat3d, Vertex>(m, "SPH_simulation_constraint_solver_stress_tensor_" + id);
		pos_ = pos;

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<double>(m, initial_vol_.get(), v) = DENSITY_SPH;
			value<Vec3>(m, initial_pos_.get(), v) = value<Vec3>(m, pos_, v);
			return true;
		});

		cgogn::parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Mat3d>(m, rotation_.get(), v) = Mat3d::Identity();
			return true;
		});

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			CellMarker<MAP, Vertex> marker(m);
			std::vector<Vertex>& n = value<std::vector<Vertex>>(m, neighborhood_.get(), v);
			n.clear();
			double& h = value<double>(m, h_.get(), v);
			h = 0;
			Vec3 pos_v1 = value<Vec3>(m, initial_pos_.get(), v);
			marker.mark(v);
			n.push_back(v);
			/*foreach_cell(m, [&](Vertex w) -> bool {
				if (marker.is_marked(w))
					return true;
				Vec3 pos_w = value<Vec3>(m, initial_pos_.get(), w);
				if (h < 1e-9)
					h = (pos_v1 - pos_w).norm();
				h = std::min(h, (pos_v1 - pos_w).norm());
				n.push_back(w);
				return true;
			});*/
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					foreach_incident_volume(m, v2, [&](Volume w2) -> bool {
						foreach_incident_vertex(m, w2, [&](Vertex v3) -> bool {
							if (!marker.is_marked(v3))
							{
								Vec3 pos_v3 = value<Vec3>(m, initial_pos_.get(), v3);
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
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
			// h *= 1.2;
			h *= 4;

			double density = 0;
			for (Vertex w : n)
			{
				Vec3 pos_v2 = value<Vec3>(m, initial_pos_.get(), w);
				density += DENSITY_SPH * Kernel_W((pos_v1 - pos_v2).norm(), h);
			}
			value<double>(m, initial_vol_.get(), v) = DENSITY_SPH / density;

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

	void rotationextraction(const Mat3d& M, Quaternion& q, int maxIter) const
	{
		for (unsigned int iter = 0; iter < maxIter; iter++)
		{
			Mat3d R = q.matrix();
			Vec3 omega =
				(R.col(0).cross(M.col(0)) + R.col(1).cross(M.col(1)) + R.col(2).cross(M.col(2))) *
				(1.0 / fabs(R.col(0).dot(M.col(0)) + R.col(1).dot(M.col(1)) + R.col(2).dot(M.col(2)) + 1.0e-9));
			double w = omega.norm();
			if (w < 1.0e-9)
				break;
			q = Quaternion(AngleAxisd(w, (1.0 / w) * omega)) * q;
			q.normalize();
		}
	}

	void test_kernel()
	{
		using Real = double;
		float eps = 1.0e-4f;
		const double supportRadius = 4.15692;
		const unsigned int numberOfSteps = 50;
		const double stepSize = static_cast<Real>(2.0) * supportRadius / (Real)(numberOfSteps - 1);
		Vec3 xi;
		xi.setZero();
		Real sum = 0.0;
		Vec3 sumV = Vec3::Zero();
		bool positive = true;
		Real V = pow(stepSize, 3);
		for (unsigned int i = 0; i < numberOfSteps; i++)
		{
			for (unsigned int j = 0; j < numberOfSteps; j++)
			{
				for (unsigned int k = 0; k < numberOfSteps; k++)
				{
					const Vec3 xj(-supportRadius + i * stepSize, -supportRadius + j * stepSize,
								  -supportRadius + k * stepSize);
					const Real W = Kernel_W((xi - xj).norm(), supportRadius);
					sum += W * V;
					sumV += gradient(xi - xj, (xi - xj).norm(), supportRadius) * V;
					if (W < -eps)
						positive = false;
				}
			}
		}
		if (fabs(sum - 1.0) < eps)
		{
			std::cout << "Kernel OK" << std::endl;
		}
		if (sumV.norm() < eps)
		{
			std::cout << "Gradient OK" << std::endl;
		}
		if (positive)
		{
			std::cout << "Kernel always positive" << std::endl;
		}
	}

	void solve_constraint(const MAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces, double) override
	{
		pos_ = pos;
		force_ = result_forces;
		compute_force_vertex(m);
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_SPH_PEER_2018_SPH_PEER_2018_H_
