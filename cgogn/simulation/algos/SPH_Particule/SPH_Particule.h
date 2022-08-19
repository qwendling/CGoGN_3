#ifndef CGOGN_SIMULATION_SPH_PARTICULE_SPH_PARTICULE_H_
#define CGOGN_SIMULATION_SPH_PARTICULE_SPH_PARTICULE_H_
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

struct Particule_SPH
{
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;

	Vec3 initial_position_;
	Vec3 current_position_;
	double initial_volume_;
	std::vector<Particule_SPH*> neighborhood_;
	Mat3d corrected_matrix_;
	Mat3d rotation_;
	double h_;
	Mat3d stress_tensor_;
	Mat3d deformation_gradient_;
	Vec3 force_;
	Vec3 speed_;
	double masse_;
	double shepard_filter_;
	std::array<Vec3, 9> RK_coeff;
	bool is_fixed;
	
	Particule_SPH(Vec3 pos, double masse)
		: initial_position_(pos), current_position_(pos),force_(0,0,0), speed_(0, 0, 0), masse_(masse),is_fixed(false)
	{
		rotation_ = std::move(Mat3d::Identity());
	}
};

template <typename MAP>
class SPH_Particule_constraint_solver : public Simulation_constraint<MAP>
{
	using Self = SPH_Particule_constraint_solver;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Volume = typename mesh_traits<MAP>::Volume;
	using Face = typename mesh_traits<MAP>::Face;
	using Quaternion = Eigen::Quaternion<double>;
	using AngleAxisd = Eigen::AngleAxis<double>;
	std::shared_ptr<Attribute<Particule_SPH*>> particule_vertex_;
	std::shared_ptr<Attribute<Particule_SPH*>> particule_volume_;
	std::shared_ptr<Attribute<Particule_SPH*>> particule_face_;
	
	std::shared_ptr<Attribute<double>> initial_volume_;
	std::shared_ptr<Attribute<Vec3>> initial_centroid_volume_;

	enum type_particule
	{
		VERTEX_PARTICULE,
		VOLUME_PARTICULE,
		VOLUME_FACE_PARTICULE
	};
	type_particule particule_type;

public:
	static inline int nb_solver = 0;
	int id;
	std::vector<Particule_SPH> particules_;

	SPH_Particule_constraint_solver() : id(nb_solver++), particule_type(VOLUME_PARTICULE)
	{
	}

	Simulation_constraint<MAP>* get_new_ptr()
	{
		return new SPH_Particule_constraint_solver<MAP>();
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

	Mat3d Corrected_matrix(Particule_SPH& p, double h) const
	{
		Mat3d Li = Mat3d::Zero();
		std::vector<Particule_SPH*>& n = p.neighborhood_;
		Vec3 xi0 = p.initial_position_;
		for (Particule_SPH* w : n)
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

	void compute_corrected_matrix()
	{
		for (Particule_SPH& p : particules_)
		{
			p.corrected_matrix_ = Corrected_matrix(p, p.h_);
		}
	}

	Vec3 Corrected_gradient(Particule_SPH& vi, Particule_SPH& vj, double h) const
	{
		Vec3 xi = vi.initial_position_;
		Vec3 xj = vj.initial_position_;
		Vec3 xij = xi - xj;

		return vi.corrected_matrix_ * gradient(xij, xij.norm(), h);
	}

	Mat3d Rotation_extraction(Particule_SPH& v, double h) const
	{
		Mat3d R;
		Mat3d F = Mat3d::Zero();

		std::vector<Particule_SPH*>& n = v.neighborhood_;
		Vec3 xi = v.current_position_;
		for (Particule_SPH* w : n)
		{
			double init_vol = w->initial_volume_;
			Vec3 xj = w->current_position_;
			Vec3 xji = xj - xi;
			Vec3 W = Corrected_gradient(v, *w, h);
			F += init_vol * xji * W.transpose();
		}
		polarDecompositionStable(F, 1.0e-6, R);

		/*Quaternion q(value<Mat3d>(m, rotation_.get(), v));
		rotationextraction(F, q, 10);
		R = q.matrix();*/

		return R;
	}

	void compute_rotated_kernel()
	{
		for (Particule_SPH& p : particules_)
		{
			p.rotation_ = Rotation_extraction(p, p.h_);
		}
	}

	Vec3 Rotated_gradient(Particule_SPH& vi, Particule_SPH& vj, double h) const
	{
		return vi.rotation_ * Corrected_gradient(vi, vj, h);
	}

	void compute_force_particules()
	{
		compute_rotated_kernel();

		for (Particule_SPH& p : particules_)
		{
			double h = p.h_;
			Mat3d Ftemp = Mat3d::Identity();
			std::vector<Particule_SPH*>& n = p.neighborhood_;
			Vec3 xi = p.current_position_;
			Vec3 xi0 = p.initial_position_;
			Mat3d& Ri = p.rotation_;
			for (Particule_SPH* w : n)
			{
				double Vj0 = w->initial_volume_;
				Vec3 xj = w->current_position_;
				Vec3 xji = xj - xi;
				Vec3 xj0 = w->initial_position_;
				Vec3 xji0 = xj0 - xi0;
				Vec3 W = Rotated_gradient(p, *w, h);
				Ftemp += Vj0 * (xji - Ri * xji0) * W.transpose();
			}
			p.deformation_gradient_ = Ftemp;
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

			p.stress_tensor_ = Pi;
		};
		for (Particule_SPH& p : particules_)
		{
			double hi = p.h_;
			Vec3 F = Vec3::Zero();
			std::vector<Particule_SPH*>& n = p.neighborhood_;
			Mat3d Pi = p.stress_tensor_;
			double Vi0 = p.initial_volume_;
			for (Particule_SPH* w : n)
			{
				double Vj0 = w->initial_volume_;
				double hj = w->h_;
				Mat3d Pj = w->stress_tensor_;
				Vec3 Wi = Rotated_gradient(p, *w, hi);
				Vec3 Wj = Rotated_gradient(*w, p, hj);
				F += Vj0 * Vi0 * (Pi * Wi - Pj * Wj);
			}
			p.force_ += F;
		};
	}

	void init_solver(MAP& m, const std::shared_ptr<Attribute<Vec3>>& init_pos,
					 const std::shared_ptr<Attribute<double>>& masse)
	{
	}

	void compute_neighborhood_Vertex( MAP& m, Attribute<Vec3>* pos)
	{
		particule_vertex_ =
			add_attribute<Particule_SPH*, Vertex>(m, "SPH_particule_constraint_solver_particule_vertex_" + id);
		particules_.reserve(nb_cells<Vertex>(m));
		foreach_cell(m, [&](Vertex v) -> bool {
			particules_.emplace_back(value<Vec3>(m, pos, v), 1);
			Particule_SPH* tmp = &(particules_.back());
			value<Particule_SPH*>(m, particule_vertex_.get(), v) = tmp;
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Particule_SPH* p = value<Particule_SPH*>(m, particule_vertex_.get(), v);
			CellMarker<MAP, Vertex> marker(m);
			std::vector<Particule_SPH*>& n = p->neighborhood_;
			n.clear();
			double& h = p->h_;
			h = 0;
			Vec3 pos_v1 = p->initial_position_;
			marker.mark(v);
			n.push_back(p);
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					foreach_incident_volume(m, v2, [&](Volume w2) -> bool {
						foreach_incident_vertex(m, w2, [&](Vertex v3) -> bool {
							if (!marker.is_marked(v3))
							{
								Particule_SPH* p3 = value<Particule_SPH*>(m, particule_vertex_.get(), v3);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_back(p3);
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

			double density = 0;
			for (Particule_SPH* w : n)
			{
				Vec3 pos_v2 = w->initial_position_;
				density += DENSITY_SPH * Kernel_W((pos_v1 - pos_v2).norm(), h);
			}
			p->initial_volume_ = DENSITY_SPH / density;
			return true;
		});
	}

	void compute_neighborhood_Volume(MAP& m, Attribute<Vec3>* pos)
	{
		particule_volume_ =
			add_attribute<Particule_SPH*, Volume>(m, "SPH_particule_constraint_solver_particule_volume_" + id);
		
		initial_volume_ = add_attribute<double, Volume>(m, "SPH_particule_constraint_solver_initial_volume_" + id);
		initial_centroid_volume_ =
			add_attribute<Vec3, Volume>(m, "SPH_particule_constraint_solver_initial_centroid_volume_" + id);
		
		geometry::compute_centroid<Vec3, Volume>(m, pos, initial_centroid_volume_.get());
		geometry::compute_volume(m, pos, initial_volume_.get());
		
		particules_.reserve(nb_cells<Volume>(m));
		
		foreach_cell(m, [&](Volume v) -> bool {
			particules_.emplace_back(value<Vec3>(m, initial_centroid_volume_.get(), v), DENSITY_SPH*value<double>(m, initial_volume_.get(), v));
			value<Particule_SPH*>(m, particule_volume_.get(), v) = &(particules_.back());
			return true;
		});
		parallel_foreach_cell(m, [&](Volume v) -> bool {
			Particule_SPH& p = *value<Particule_SPH*>(m, particule_volume_.get(), v);
			p.initial_volume_ = value<double>(m, initial_volume_.get(),v);
			CellMarker<MAP, Volume> marker(m);
			std::vector<Particule_SPH*>& n = p.neighborhood_;
			n.clear();
			double& h = p.h_;
			h = 0;
			Vec3 pos_v1 = p.initial_position_;
			marker.mark(v);
			n.push_back(&p);
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				foreach_incident_volume(m, w, [&](Volume v2) -> bool {
					foreach_incident_vertex(m, v2, [&](Vertex w2) -> bool {
						foreach_incident_volume(m, w2, [&](Volume v3) -> bool {
							if (!marker.is_marked(v3))
							{
								Particule_SPH* p3 = value<Particule_SPH*>(m, particule_volume_.get(), v3);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_back(p3);
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
			return true;
		});
		particule_vertex_ =
			add_attribute<Particule_SPH*, Vertex>(m, "SPH_particule_constraint_solver_particule_vertex_" + id);
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Particule_SPH*>(m, particule_vertex_.get(), v) = new Particule_SPH(value<Vec3>(m, pos, v),1);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Particule_SPH& p = *value<Particule_SPH*>(m, particule_vertex_.get(), v);
			CellMarker<MAP, Volume> marker(m);
			std::vector<Particule_SPH*>& n = p.neighborhood_;
			n.clear();
			double& h = p.h_;
			h = 0;
			Vec3 pos_v1 = p.initial_position_;
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					foreach_incident_volume(m, v2, [&](Volume w2) -> bool {
							if (!marker.is_marked(w2))
							{
								Particule_SPH* p3 = value<Particule_SPH*>(m, particule_volume_.get(), w2);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_back(p3);
								marker.mark(w2);
							}

						return true;
					});
					return true;
				});
				return true;
			});
			h *= 4;

			double sk = 0;
			for (Particule_SPH* w : n)
			{
				Vec3 pos_v2 = w->initial_position_;
				sk += w->initial_volume_* Kernel_W((pos_v1 - pos_v2).norm(), h);
			}
			p.shepard_filter_ = 1.0f / sk;
			return true;
		});
	}

	void compute_neighborhood_Volume_Face(MAP& m, Attribute<Vec3>* pos)
	{
		particule_volume_ =
			add_attribute<Particule_SPH*, Volume>(m, "SPH_particule_constraint_solver_particule_volume_" + id);
		particule_face_ =
			add_attribute<Particule_SPH*, Face>(m, "SPH_particule_constraint_solver_particule_face_" + id);
		
		initial_volume_ = add_attribute<double, Volume>(m, "SPH_particule_constraint_solver_initial_volume_" + id);
		initial_centroid_volume_ =
			add_attribute<Vec3, Volume>(m, "SPH_particule_constraint_solver_initial_centroid_volume_" + id);
		
		geometry::compute_centroid<Vec3, Volume>(m, pos, initial_centroid_volume_.get());
		geometry::compute_volume(m, pos, initial_volume_.get());
		
		particules_.reserve(nb_cells<Volume>(m)+nb_cells<Face>(m));
		
		foreach_cell(m, [&](Volume v) -> bool {
			double masse_volume = DENSITY_SPH*value<double>(m, initial_volume_.get(), v);
			foreach_incident_face(m,v,[&](Face f)->bool{
				if(is_boundary(m,phi3(m,f.dart))){
					particules_.emplace_back(geometry::centroid<Vec3>(m,f,pos), 0.25*masse_volume);
					Particule_SPH* tmp = &(particules_.back());
					value<Particule_SPH*>(m, particule_face_.get(), f) = tmp;
					masse_volume *= 0.75;
					//ajout dans son propre voisinage
					tmp->neighborhood_.push_back(tmp);
					return true;
				}
				return true;
			});
			particules_.emplace_back(value<Vec3>(m, initial_centroid_volume_.get(), v), masse_volume);
			value<Particule_SPH*>(m, particule_volume_.get(), v) = &(particules_.back());
			return true;
		});
		parallel_foreach_cell(m, [&](Volume v) -> bool {
			Particule_SPH& p = *value<Particule_SPH*>(m, particule_volume_.get(), v);
			p.initial_volume_ = value<double>(m, initial_volume_.get(),v);
			CellMarker<MAP, Volume> marker(m);
			std::vector<Particule_SPH*>& n = p.neighborhood_;
			n.clear();
			double& h = p.h_;
			h = 0;
			Vec3 pos_v1 = p.initial_position_;
			marker.mark(v);
			n.push_back(&p);
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				foreach_incident_volume(m, w, [&](Volume v2) -> bool {
					foreach_incident_vertex(m, v2, [&](Vertex w2) -> bool {
						foreach_incident_volume(m, w2, [&](Volume v3) -> bool {
							
							
							
							if (!marker.is_marked(v3))
							{
								foreach_incident_face(m,v3,[&](Face f)->bool{
									if(is_boundary(m,phi3(m,f.dart))){
										Particule_SPH* pf = value<Particule_SPH*>(m, particule_face_.get(),f);
										Vec3 pos_pf = pf->initial_position_;
										if (h < 1e-9)
											h = (pos_v1 - pos_pf).norm();
										h = std::min(h, (pos_v1 - pos_pf).norm());
										n.push_back(pf);
										pf->neighborhood_.push_back(&p);
										return true;
									}
									return true;
								});
								Particule_SPH* p3 = value<Particule_SPH*>(m, particule_volume_.get(), v3);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_back(p3);
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
			return true;
		});
		particule_vertex_ =
			add_attribute<Particule_SPH*, Vertex>(m, "SPH_particule_constraint_solver_particule_vertex_" + id);
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Particule_SPH*>(m, particule_vertex_.get(), v) = new Particule_SPH(value<Vec3>(m, pos, v),1);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Particule_SPH& p = *value<Particule_SPH*>(m, particule_vertex_.get(), v);
			CellMarker<MAP, Volume> marker(m);
			std::vector<Particule_SPH*>& n = p.neighborhood_;
			n.clear();
			double& h = p.h_;
			h = 0;
			Vec3 pos_v1 = p.initial_position_;
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					foreach_incident_volume(m, v2, [&](Volume w2) -> bool {
							if (!marker.is_marked(w2))
							{
								foreach_incident_face(m,w2,[&](Face f)->bool{
									if(is_boundary(m,phi3(m,f.dart))){
										Particule_SPH* pf = value<Particule_SPH*>(m, particule_face_.get(),f);
										Vec3 pos_pf = pf->initial_position_;
										if (h < 1e-9)
											h = (pos_v1 - pos_pf).norm();
										h = std::min(h, (pos_v1 - pos_pf).norm());
										n.push_back(pf);
										return true;
									}
									return true;
								});
								Particule_SPH* p3 = value<Particule_SPH*>(m, particule_volume_.get(), w2);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_back(p3);
								marker.mark(w2);
							}

						return true;
					});
					return true;
				});
				return true;
			});
			h *= 4;

			double sk = 0;
			for (Particule_SPH* w : n)
			{
				Vec3 pos_v2 = w->initial_position_;
				sk += w->initial_volume_* Kernel_W((pos_v1 - pos_v2).norm(), h);
			}
			p.shepard_filter_ = 1.0f / sk;
			return true;
		});
	}

	void init_solver(MAP& m, Attribute<Vec3>* pos)
	{
		switch (particule_type)
		{
		case VERTEX_PARTICULE:
			compute_neighborhood_Vertex(m, pos);
			break;
		case VOLUME_PARTICULE:
			compute_neighborhood_Volume(m, pos);
			break;
		case VOLUME_FACE_PARTICULE:
			compute_neighborhood_Volume_Face(m, pos);
			break;
		}
		compute_corrected_matrix();
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
	
	void SPH_skinning(const MAP& m, Attribute<Vec3>* pos){
		if(particule_type == VERTEX_PARTICULE){
			parallel_foreach_cell(m,[&](Vertex v){
				value<Vec3>(m, pos, v) = value<Particule_SPH*>(m, particule_vertex_.get(), v)->current_position_;
				return true;
			});
		}else{
			parallel_foreach_cell(m,[&](Vertex v){
				Particule_SPH* p = value<Particule_SPH*>(m, particule_vertex_.get(), v);
				Vec3 new_pos = Vec3::Zero();
				for(Particule_SPH* p2 : p->neighborhood_ ){
					new_pos += p2->initial_volume_
							*(p2->deformation_gradient_*(p->initial_position_-p2->initial_position_)+p2->current_position_)
							*Kernel_W((p->initial_position_ - p2->initial_position_).norm(), p->h_);
				}
				value<Vec3>(m, pos, v) = p->shepard_filter_*new_pos;
				return true;
			});
		}
	}
	
	void particule_integration(double time_step){
		for (auto& p : particules_)
		{
			p.RK_coeff[8] = p.force_;
		}
		compute_force_particules();
		for (auto& p : particules_)
		{
			if(p.is_fixed){
				continue;
			}
			// k1
			p.RK_coeff[0] =
				(time_step * p.force_ / p.masse_);
			// j1
			p.RK_coeff[1] =
				time_step * (0.995 * p.speed_ + p.RK_coeff[0]);
	
			p.current_position_ += p.RK_coeff[1] / 2.0f;
	
			p.force_ = p.RK_coeff[8];
		}
		compute_force_particules();
		for (auto& p : particules_)
		{
			if(p.is_fixed){
				continue;
			}
			// k2
			p.RK_coeff[2] =
				(time_step * p.force_ / p.masse_);
			// j2
			p.RK_coeff[3] =
				time_step *
				(0.995 * p.speed_ + p.RK_coeff[2] / 2.0f);
	
			p.current_position_ -= p.RK_coeff[1] / 2.0f;
	
			p.current_position_ += p.RK_coeff[3] / 2.0f;
	
			p.force_ = p.RK_coeff[8];
		}
		compute_force_particules();
		for (auto& p : particules_)
		{
			if(p.is_fixed){
				continue;
			}
			// k3
			p.RK_coeff[4] =
				(time_step * p.force_ / p.masse_);
			// j3
			p.RK_coeff[5] =
				time_step *
				(0.995 * p.speed_ + p.RK_coeff[4] / 2.0f);
	
			p.current_position_ -= p.RK_coeff[3] / 2.0f;
	
			p.current_position_ += p.RK_coeff[5];
	
			p.force_ = p.RK_coeff[8];
		}
		compute_force_particules();
		for (auto& p : particules_)
		{
			if(p.is_fixed){
				continue;
			}
			// k4
			p.RK_coeff[6] =
				(time_step * p.force_ / p.masse_);
			// j4
			p.RK_coeff[7] =
				time_step * (0.995 * p.speed_ + p.RK_coeff[6]);
	
			p.current_position_ -= p.RK_coeff[5];
	
			Vec3 k1 = p.RK_coeff[0];
			Vec3 k2 = p.RK_coeff[2];
			Vec3 k3 = p.RK_coeff[4];
			Vec3 k4 = p.RK_coeff[6];
	
			Vec3 j1 = p.RK_coeff[1];
			Vec3 j2 = p.RK_coeff[3];
			Vec3 j3 = p.RK_coeff[5];
			Vec3 j4 = p.RK_coeff[7];
	
			Vec3 diff_pos = 1. / 6. * (j1 + (2 * j2) + (2 * j3) + j4);
	
			Vec3 diff_speed = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
	
			if (diff_pos.norm() < 1.0e-10)
				diff_pos = Vec3(0, 0, 0);
	
			if (diff_speed.norm() < 1.0e-10)
				diff_speed = Vec3(0, 0, 0);
	
			p.current_position_ += diff_pos;
	
			p.speed_ = 0.995 * p.speed_ + diff_speed;
			p.force_ = Vec3(0, 0, 0);
		}
	}
	
	template <typename FUNC>
	void set_particule_fixed(const FUNC& fn){
		for(Particule_SPH& p:particules_){
			p.is_fixed = fn(p);
		}
	}
	
	template <typename FUNC>
	void set_particule_forces(const FUNC& fn){
		for(Particule_SPH& p:particules_){
			p.force_ = fn(p);
		}
	}
 
	void solve_constraint(const MAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* , double time_step) override
	{
		particule_integration(time_step);
		SPH_skinning(m, pos);
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_SPH_PARTICULE_SPH_PARTICULE_H_
