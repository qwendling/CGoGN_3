#ifndef CGOGN_SIMULATION_SHAPE_MATCHING_SHAPE_MATCHING_H_
#define CGOGN_SIMULATION_SHAPE_MATCHING_SHAPE_MATCHING_H_
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class shape_matching_constraint_solver : public Simulation_constraint<MAP>
{
	using Self = shape_matching_constraint_solver;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MAP>::Vertex;

public:
	std::shared_ptr<Attribute<Vec3>> vertex_init_position_;
	std::shared_ptr<Attribute<double>> masse_;
	std::shared_ptr<Attribute<Vec3>> q_;
	std::shared_ptr<Attribute<Vec3>> goals_;
	Vec3 init_cm_;
	double stiffness_;

	shape_matching_constraint_solver(double stiffness)
		: vertex_init_position_(nullptr), q_(nullptr), stiffness_(stiffness)
	{
	}

	void init_solver(MAP& m, const std::shared_ptr<Attribute<Vec3>>& init_pos,
					 const std::shared_ptr<Attribute<double>>& masse)
	{
		vertex_init_position_ = init_pos;
		masse_ = masse;
		q_ = get_attribute<Vec3, Vertex>(m, "shape_matching_constraint_solver_q");
		if (q_ == nullptr)
			q_ = add_attribute<Vec3, Vertex>(m, "shape_matching_constraint_solver_q");
		goals_ = get_attribute<Vec3, Vertex>(m, "shape_matching_constraint_solver_goals");
		if (goals_ == nullptr)
			goals_ = add_attribute<Vec3, Vertex>(m, "shape_matching_constraint_solver_goals");
		init_cm_ = Vec3(0, 0, 0);
		double masse_totale = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			init_cm_ += value<double>(m, masse_, v) * value<Vec3>(m, vertex_init_position_.get(), v);
			masse_totale += value<double>(m, masse_, v);
			return true;
		});
		init_cm_ /= masse_totale;
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, q_.get(), v) = value<Vec3>(m, vertex_init_position_.get(), v) - init_cm_;
			return true;
		});
	}

	void update_topo(const MAP& m)
	{
		init_cm_ = Vec3(0, 0, 0);
		double masse_totale = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			init_cm_ += value<double>(m, masse_, v) * value<Vec3>(m, vertex_init_position_.get(), v);
			masse_totale += value<double>(m, masse_, v);
			return true;
		});
		init_cm_ /= masse_totale;
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, q_.get(), v) = value<Vec3>(m, vertex_init_position_.get(), v) - init_cm_;
			return true;
		});
	}

	double oneNorm(const Mat3d& A)
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

	double infNorm(const Mat3d& A)
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

	void polarDecompositionStable(const Mat3d& M, double tolerance, Mat3d& R)
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

	void solve_constraint(const MAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces, double time_step) override
	{
		Mat3d Apq = Mat3d::Zero();
		Vec3 cm = Vec3(0, 0, 0);
		double masse_total = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			cm += value<double>(m, masse_, v) * value<Vec3>(m, vertex_init_position_.get(), v);
			masse_total += value<double>(m, masse_, v);
			return true;
		});
		cm /= masse_total;
		foreach_cell(m, [&](Vertex v) -> bool {
			Vec3 p = value<Vec3>(m, pos, v) - cm;
			Apq += value<double>(m, masse_, v) * p * value<Vec3>(m, q_.get(), v).transpose();
			return true;
		});
		Mat3d R;
		polarDecompositionStable(Apq, 1.0e-6, R);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, goals_.get(), v) = R * value<Vec3>(m, q_.get(), v) + cm;
			value<Vec3>(m, result_forces, v) = value<double>(m, masse_, v) * stiffness_ *
											   (value<Vec3>(m, goals_.get(), v) - value<Vec3>(m, pos, v)) /
											   (time_step * time_step);
			return true;
		});
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SHAPE_MATCHING_SHAPE_MATCHING_H_
