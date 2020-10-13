#ifndef CGOGN_SIMULATION_LATTICE_SHAPE_MATCHING_LATTICE_SHAPE_MATCHING_H_
#define CGOGN_SIMULATION_LATTICE_SHAPE_MATCHING_LATTICE_SHAPE_MATCHING_H_
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class lattice_shape_matching_constraint_solver : public Simulation_constraint<MAP>
{
	using Self = lattice_shape_matching_constraint_solver;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Volume = typename mesh_traits<MAP>::Volume;

public:
	std::shared_ptr<Attribute<std::vector<Vertex>>> vertex_region_;
	std::shared_ptr<Attribute<double>> masse_region_;
	std::shared_ptr<Attribute<double>> modify_masse_vertex_;
	std::shared_ptr<Attribute<Vec3>> init_cm_region_;
	std::shared_ptr<Attribute<Vec3>> cm_region_;
	std::shared_ptr<Attribute<Vec3>> vertex_init_position_;
	std::shared_ptr<Attribute<Vec3>> goals_;
	std::shared_ptr<Attribute<Mat3d>> A_;
	std::shared_ptr<Attribute<Mat3d>> R_;
	std::shared_ptr<Attribute<Vec3>> translate_;
	std::shared_ptr<Attribute<Vec3>> translate_vertex_;
	double stiffness_;
	int id;
	static inline int nb_solver = 0;

	lattice_shape_matching_constraint_solver(double stiffness)
		: vertex_init_position_(nullptr), vertex_region_(nullptr), masse_region_(nullptr),
		  modify_masse_vertex_(nullptr), init_cm_region_(nullptr), stiffness_(stiffness), id(nb_solver++)
	{
	}

	Simulation_constraint<MAP>* get_new_ptr()
	{
		return new lattice_shape_matching_constraint_solver<MAP>(stiffness_);
	}

	void init_solver(MAP& m, const std::shared_ptr<Attribute<Vec3>>& init_pos,
					 const std::shared_ptr<Attribute<double>>& masse)
	{
		vertex_init_position_ = init_pos;
		this->masse_ = masse;
		init_region(m);
	}

	void init_solver(MAP& m, Attribute<Vec3>* pos)
	{
		vertex_init_position_ =
			get_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_vertex_init_position" + id);
		if (vertex_init_position_ == nullptr)
		{
			vertex_init_position_ =
				add_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_vertex_init_position" + id);
			parallel_foreach_cell(static_cast<CMap3&>(m), [&](Vertex v) -> bool {
				value<Vec3>(m, vertex_init_position_.get(), v) = value<Vec3>(m, pos, v);
				return true;
			});
		}
		this->masse_ = get_attribute<double, Vertex>(m, "lattice_shape_matching_constraint_solver_masse" + id);
		if (this->masse_ == nullptr)
		{
			this->masse_ = add_attribute<double, Vertex>(m, "lattice_shape_matching_constraint_solver_masse" + id);
			foreach_cell(m, [&](Volume v) -> bool {
				double vol = geometry::volume(m, v, vertex_init_position_.get());
				std::vector<Vertex> inc_vertices;
				foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
					inc_vertices.push_back(w);
					return true;
				});
				for (auto w : inc_vertices)
				{
					value<double>(m, this->masse_.get(), w) += vol / inc_vertices.size();
				}
				return true;
			});
		}

		init_region(m);
	}

	void init_region(MAP& m)
	{
		translate_ = get_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_translate" + id);
		if (translate_ == nullptr)
			translate_ = add_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_translate" + id);

		translate_vertex_ =
			get_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_translate_vertex" + id);
		if (translate_vertex_ == nullptr)
			translate_vertex_ =
				add_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_translate_vertex" + id);

		R_ = get_attribute<Mat3d, Vertex>(m, "lattice_shape_matching_constraint_solver_R" + id);
		if (R_ == nullptr)
			R_ = add_attribute<Mat3d, Vertex>(m, "lattice_shape_matching_constraint_solver_R" + id);

		A_ = get_attribute<Mat3d, Vertex>(m, "lattice_shape_matching_constraint_solver_A" + id);
		if (A_ == nullptr)
			A_ = add_attribute<Mat3d, Vertex>(m, "lattice_shape_matching_constraint_solver_A" + id);

		goals_ = get_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_goals" + id);
		if (goals_ == nullptr)
			goals_ = add_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_goals" + id);

		vertex_region_ = get_attribute<std::vector<Vertex>, Vertex>(
			m, "lattice_shape_matching_constraint_solver_vertex_region" + id);
		if (vertex_region_ == nullptr)
			vertex_region_ = add_attribute<std::vector<Vertex>, Vertex>(
				m, "lattice_shape_matching_constraint_solver_vertex_region" + id);

		masse_region_ = get_attribute<double, Vertex>(m, "lattice_shape_matching_constraint_solver_masse_region" + id);
		if (masse_region_ == nullptr)
			masse_region_ =
				add_attribute<double, Vertex>(m, "lattice_shape_matching_constraint_solver_masse_region" + id);

		modify_masse_vertex_ =
			get_attribute<double, Vertex>(m, "lattice_shape_matching_constraint_solver_modify_masse_vertex" + id);
		if (modify_masse_vertex_ == nullptr)
			modify_masse_vertex_ =
				add_attribute<double, Vertex>(m, "lattice_shape_matching_constraint_solver_modify_masse_vertex" + id);

		init_cm_region_ =
			get_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_init_cm_region" + id);
		if (init_cm_region_ == nullptr)
			init_cm_region_ =
				add_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_init_cm_region" + id);

		cm_region_ = get_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_cm_region" + id);
		if (cm_region_ == nullptr)
			cm_region_ = add_attribute<Vec3, Vertex>(m, "lattice_shape_matching_constraint_solver_cm_region" + id);

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<std::vector<Vertex>>(m, vertex_region_.get(), v) = {};
			CellMarkerStore<MAP, Vertex> marker(m);
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					if (!marker.is_marked(v2))
					{
						marker.mark(v2);
						value<std::vector<Vertex>>(m, vertex_region_.get(), v).push_back(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<double>(m, modify_masse_vertex_.get(), v) =
				value<double>(m, this->masse_.get(), v) /
				(double)value<std::vector<Vertex>>(m, vertex_region_.get(), v).size();
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<double>(m, masse_region_.get(), v) = 0;
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<double>(m, masse_region_.get(), v) += value<double>(m, modify_masse_vertex_.get(), v2);
			}
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, init_cm_region_.get(), v) = Vec3::Zero();
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<Vec3>(m, init_cm_region_.get(), v) +=
					value<double>(m, modify_masse_vertex_.get(), v2) * value<Vec3>(m, vertex_init_position_.get(), v);
			}
			value<Vec3>(m, init_cm_region_.get(), v) /= value<double>(m, masse_region_.get(), v);
			return true;
		});
	}

	void update_topo(const MAP& m, const std::vector<Vertex>&)
	{
		std::cout << "________________________________________________________________" << std::endl;
		foreach_cell(m, [&](Volume v) -> bool {
			double vol = geometry::volume(m, v, vertex_init_position_.get());
			std::vector<Vertex> inc_vertices;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				inc_vertices.push_back(w);
				return true;
			});
			for (auto w : inc_vertices)
			{
				value<double>(m, this->masse_.get(), w) += vol / inc_vertices.size();
			}
			return true;
		});

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<std::vector<Vertex>>(m, vertex_region_.get(), v) = {};
			CellMarkerStore<MAP, Vertex> marker(m);
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					if (!marker.is_marked(v2))
					{
						marker.mark(v2);
						value<std::vector<Vertex>>(m, vertex_region_.get(), v).push_back(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<double>(m, modify_masse_vertex_.get(), v) =
				value<double>(m, this->masse_.get(), v) /
				(double)value<std::vector<Vertex>>(m, vertex_region_.get(), v).size();
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<double>(m, masse_region_.get(), v) = 0;
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<double>(m, masse_region_.get(), v) += value<double>(m, modify_masse_vertex_.get(), v2);
			}
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, init_cm_region_.get(), v) = Vec3::Zero();
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<Vec3>(m, init_cm_region_.get(), v) +=
					value<double>(m, modify_masse_vertex_.get(), v2) * value<Vec3>(m, vertex_init_position_.get(), v2);
			}
			value<Vec3>(m, init_cm_region_.get(), v) /= value<double>(m, masse_region_.get(), v);
			std::cout << "value<Vec3>(m, init_cm_region_.get(), v)" << std::endl;
			std::cout << value<Vec3>(m, init_cm_region_.get(), v) << std::endl;
			return true;
		});
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

	void solve_constraint(const MAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* result_forces, double time_step) override
	{
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, cm_region_.get(), v) = Vec3(0, 0, 0);
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<Vec3>(m, cm_region_.get(), v) +=
					value<double>(m, modify_masse_vertex_.get(), v2) * value<Vec3>(m, pos, v2);
			}
			value<Vec3>(m, cm_region_.get(), v) /= value<double>(m, masse_region_.get(), v);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Mat3d>(m, A_.get(), v) = Mat3d::Zero();
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<Mat3d>(m, A_.get(), v) += value<double>(m, modify_masse_vertex_.get(), v2) *
												value<Vec3>(m, pos, v2) *
												value<Vec3>(m, vertex_init_position_.get(), v2).transpose();
			}
			value<Mat3d>(m, A_.get(), v) -= value<double>(m, masse_region_.get(), v) *
											value<Vec3>(m, cm_region_.get(), v) *
											value<Vec3>(m, init_cm_region_.get(), v).transpose();
			polarDecompositionStable(value<Mat3d>(m, A_.get(), v), 1.0e-6, value<Mat3d>(m, A_.get(), v));

			value<Vec3>(m, translate_.get(), v) =
				value<Vec3>(m, cm_region_.get(), v) -
				value<Mat3d>(m, A_.get(), v) * value<Vec3>(m, init_cm_region_.get(), v);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Mat3d>(m, R_.get(), v) = Mat3d::Zero();
			value<Vec3>(m, translate_vertex_.get(), v) = Vec3::Zero();
			auto r = value<std::vector<Vertex>>(m, vertex_region_.get(), v);
			for (Vertex v2 : r)
			{
				value<Mat3d>(m, R_.get(), v) += value<Mat3d>(m, A_.get(), v2);
				value<Vec3>(m, translate_vertex_.get(), v) += value<Vec3>(m, translate_.get(), v2);
			}
			value<Mat3d>(m, R_.get(), v) /= (double)r.size();
			value<Vec3>(m, translate_vertex_.get(), v) /= (double)r.size();
			value<Vec3>(m, goals_.get(), v) =
				value<Mat3d>(m, R_.get(), v) * value<Vec3>(m, vertex_init_position_.get(), v) +
				value<Vec3>(m, translate_vertex_.get(), v);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Vec3 result = value<double>(m, this->masse_, v) * stiffness_ *
						  (value<Vec3>(m, goals_.get(), v) - value<Vec3>(m, pos, v)) / (time_step * time_step);
			if (result.norm() > 1.0e-10)
				value<Vec3>(m, result_forces, v) += result;
			return true;
		});
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_LATTICE_SHAPE_MATCHING_LATTICE_SHAPE_MATCHING_H_
