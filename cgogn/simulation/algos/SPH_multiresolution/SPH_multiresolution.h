#ifndef CGOGN_SIMULATION_SPH_MULTIRESOLUTION_SPH_MULTIRESOLUTION_H_
#define CGOGN_SIMULATION_SPH_MULTIRESOLUTION_SPH_MULTIRESOLUTION_H_
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>
#include <forward_list>

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

struct Particule_SPH_MR
{
	using Vec3 = geometry::Vec3; 
	using Mat3d = geometry::Mat3d;

	Vec3 initial_position_;
	Vec3 current_position_;
	double initial_volume_;
	std::forward_list<Particule_SPH_MR*> neighborhood_;
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
	std::forward_list<Particule_SPH_MR*> child_;
	
	Particule_SPH_MR(Vec3 pos, double masse)
		: initial_position_(pos), current_position_(pos),force_(0,0,0), speed_(0, 0, 0), masse_(masse),is_fixed(false)
	{
		rotation_ = std::move(Mat3d::Identity());
		deformation_gradient_ = std::move(Mat3d::Identity());
	}
};

template <typename MRMAP>
class SPH_Multiresolution_constraint_solver : public Simulation_constraint<MRMAP>
{
	using Self = SPH_Multiresolution_constraint_solver;
	template <typename T>
	using Attribute = typename mesh_traits<MRMAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MRMAP>::Vertex;
	using Volume = typename mesh_traits<MRMAP>::Volume;
	using Face = typename mesh_traits<MRMAP>::Face;
	using Quaternion = Eigen::Quaternion<double>;
	using AngleAxisd = Eigen::AngleAxis<double>;
	std::shared_ptr<Attribute<Particule_SPH_MR*>> particule_vertex_;
	std::shared_ptr<Attribute<Particule_SPH_MR*>> particule_volume_;
	std::shared_ptr<Attribute<Particule_SPH_MR*>> particule_face_;
	
	std::shared_ptr<Attribute<double>> initial_volume_;
	std::shared_ptr<Attribute<Vec3>> initial_centroid_volume_;

	enum type_particule
	{
		VOLUME_PARTICULE,
		VOLUME_FACE_PARTICULE
	};
	type_particule particule_type;

public:
	static inline int nb_solver = 0;
	int id;
	std::vector<Particule_SPH_MR*> particules_;

	SPH_Multiresolution_constraint_solver() : id(nb_solver++), particule_type(VOLUME_PARTICULE)
	{
	}

	Simulation_constraint<MRMAP>* get_new_ptr()
	{
		return new SPH_Multiresolution_constraint_solver<MRMAP>();
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

	Mat3d Corrected_matrix(Particule_SPH_MR& p, double h) const
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

	void compute_corrected_matrix()
	{
		for (Particule_SPH_MR* p : particules_)
		{
			p->corrected_matrix_ = Corrected_matrix(*p, p->h_);
		}
	}

	Vec3 Corrected_gradient(Particule_SPH_MR& vi, Particule_SPH_MR& vj, double h) const
	{
		Vec3 xi = vi.initial_position_;
		Vec3 xj = vj.initial_position_;
		Vec3 xij = xi - xj;

		return vi.corrected_matrix_ * gradient(xij, xij.norm(), h);
	}

	Mat3d Rotation_extraction(Particule_SPH_MR& v, double h) const
	{
		Mat3d R;
		Mat3d F = Mat3d::Zero();

		std::forward_list<Particule_SPH_MR*>& n = v.neighborhood_;
		Vec3 xi = v.current_position_;
		for (Particule_SPH_MR* w : n)
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
		for (Particule_SPH_MR* p : particules_)
		{
			p->rotation_ = Rotation_extraction(*p, p->h_);
		}
	}

	Vec3 Rotated_gradient(Particule_SPH_MR& vi, Particule_SPH_MR& vj, double h) const
	{
		return vi.rotation_ * Corrected_gradient(vi, vj, h);
	}

	void compute_force_particules()
	{
		compute_rotated_kernel();

		for (Particule_SPH_MR* p : particules_)
		{
			double h = p->h_;
			Mat3d Ftemp = Mat3d::Identity();
			std::forward_list<Particule_SPH_MR*>& n = p->neighborhood_;
			Vec3 xi = p->current_position_;
			Vec3 xi0 = p->initial_position_;
			Mat3d& Ri = p->rotation_;
			for (Particule_SPH_MR* w : n)
			{
				double Vj0 = w->initial_volume_;
				Vec3 xj = w->current_position_;
				Vec3 xji = xj - xi;
				Vec3 xj0 = w->initial_position_;
				Vec3 xji0 = xj0 - xi0;
				Vec3 W = Rotated_gradient(*p, *w, h);
				Ftemp += Vj0 * (xji - Ri * xji0) * W.transpose();
			}
			p->deformation_gradient_ = Ftemp;
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

			p->stress_tensor_ = Pi;
		};
		for (Particule_SPH_MR* p : particules_)
		{
			double hi = p->h_;
			Vec3 F = Vec3::Zero();
			std::forward_list<Particule_SPH_MR*>& n = p->neighborhood_;
			Mat3d Pi = p->stress_tensor_;
			double Vi0 = p->initial_volume_;
			for (Particule_SPH_MR* w : n)
			{
				double Vj0 = w->initial_volume_;
				double hj = w->h_;
				Mat3d Pj = w->stress_tensor_;
				Vec3 Wi = Rotated_gradient(*p, *w, hi);
				Vec3 Wj = Rotated_gradient(*w, *p, hj);
				F += Vj0 * Vi0 * (Pi * Wi - Pj * Wj);
			}
			p->force_ += F;
		};
	}

	void init_solver(MRMAP& m, const std::shared_ptr<Attribute<Vec3>>& init_pos,
					 const std::shared_ptr<Attribute<double>>& masse)
	{
	}

	void compute_neighborhood_Volume(MRMAP& m, Attribute<Vec3>* pos)
	{
		init_particule_MR(m,pos);
		std::function<void(Particule_SPH_MR*)> fn;
		fn = [&](Particule_SPH_MR* p){
			if(p->child_.empty())
				return;
			for(auto c:p->child_){
				c->h_ = 0.85*p->h_;
				fn(c);
			}
		};
		parallel_foreach_cell(m, [&](Volume v) -> bool {
			Particule_SPH_MR& p = *value<Particule_SPH_MR*>(m, particule_volume_.get(), v);
			p.initial_volume_ = value<double>(m, initial_volume_.get(),v);
			CellMarker<MRMAP, Volume> marker(m);
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
		particule_vertex_ =
			add_attribute<Particule_SPH_MR*, Vertex>(m, "SPH_particule_constraint_solver_particule_vertex_" + id);
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Particule_SPH_MR*>(m, particule_vertex_.get(), v) = new Particule_SPH_MR(value<Vec3>(m, pos, v),1);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Particule_SPH_MR& p = *value<Particule_SPH_MR*>(m, particule_vertex_.get(), v);
			CellMarker<MRMAP, Volume> marker(m);
			std::forward_list<Particule_SPH_MR*>& n = p.neighborhood_;
			n.clear();
			double& h = p.h_;
			h = 0;
			Vec3 pos_v1 = p.initial_position_;
			foreach_incident_volume(m, v, [&](Volume w) -> bool {
				foreach_incident_vertex(m, w, [&](Vertex v2) -> bool {
					foreach_incident_volume(m, v2, [&](Volume w2) -> bool {
							if (!marker.is_marked(w2))
							{
								Particule_SPH_MR* p3 = value<Particule_SPH_MR*>(m, particule_volume_.get(), w2);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_front(p3);
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
			for (Particule_SPH_MR* w : n)
			{
				Vec3 pos_v2 = w->initial_position_;
				sk += w->initial_volume_* Kernel_W((pos_v1 - pos_v2).norm(), h);
			}
			p.shepard_filter_ = 1.0f / sk;
			return true;
		});
	}

	void compute_neighborhood_Volume_Face(MRMAP& m, Attribute<Vec3>* pos)
	{
		particule_volume_ =
			add_attribute<Particule_SPH_MR*, Volume>(m, "SPH_particule_constraint_solver_particule_volume_" + id);
		particule_face_ =
			add_attribute<Particule_SPH_MR*, Face>(m, "SPH_particule_constraint_solver_particule_face_" + id);
		
		initial_volume_ = add_attribute<double, Volume>(m, "SPH_particule_constraint_solver_initial_volume_" + id);
		initial_centroid_volume_ =
			add_attribute<Vec3, Volume>(m, "SPH_particule_constraint_solver_initial_centroid_volume_" + id);
		
		geometry::compute_centroid<Vec3, Volume>(m, pos, initial_centroid_volume_.get());
		geometry::compute_volume(m, pos, initial_volume_.get());
		
		
		foreach_cell(m, [&](Volume v) -> bool {
			double masse_volume = DENSITY_SPH*value<double>(m, initial_volume_.get(), v);
			foreach_incident_face(m,v,[&](Face f)->bool{
				if(is_boundary(m,phi3(m,f.dart))){
					particules_.push_back(new Particule_SPH_MR(geometry::centroid<Vec3>(m,f,pos), 0.25*masse_volume));
					Particule_SPH_MR* tmp = particules_.back();
					value<Particule_SPH_MR*>(m, particule_face_.get(), f) = tmp;
					masse_volume *= 0.75;
					//ajout dans son propre voisinage
					tmp->neighborhood_.push_front(tmp);
					return true;
				}
				return true;
			});
			particules_.push_back(new Particule_SPH_MR(value<Vec3>(m, initial_centroid_volume_.get(), v), masse_volume));
			value<Particule_SPH_MR*>(m, particule_volume_.get(), v) = particules_.back();
			return true;
		});
		parallel_foreach_cell(m, [&](Volume v) -> bool {
			Particule_SPH_MR& p = *value<Particule_SPH_MR*>(m, particule_volume_.get(), v);
			p.initial_volume_ = value<double>(m, initial_volume_.get(),v);
			CellMarker<MRMAP, Volume> marker(m);
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
								foreach_incident_face(m,v3,[&](Face f)->bool{
									if(is_boundary(m,phi3(m,f.dart))){
										Particule_SPH_MR* pf = value<Particule_SPH_MR*>(m, particule_face_.get(),f);
										Vec3 pos_pf = pf->initial_position_;
										if (h < 1e-9)
											h = (pos_v1 - pos_pf).norm();
										h = std::min(h, (pos_v1 - pos_pf).norm());
										n.push_front(pf);
										pf->neighborhood_.push_front(&p);
										return true;
									}
									return true;
								});
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
			return true;
		});
		particule_vertex_ =
			add_attribute<Particule_SPH_MR*, Vertex>(m, "SPH_particule_constraint_solver_particule_vertex_" + id);
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Particule_SPH_MR*>(m, particule_vertex_.get(), v) = new Particule_SPH_MR(value<Vec3>(m, pos, v),1);
			return true;
		});
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Particule_SPH_MR& p = *value<Particule_SPH_MR*>(m, particule_vertex_.get(), v);
			CellMarker<MRMAP, Volume> marker(m);
			std::forward_list<Particule_SPH_MR*>& n = p.neighborhood_;
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
										Particule_SPH_MR* pf = value<Particule_SPH_MR*>(m, particule_face_.get(),f);
										Vec3 pos_pf = pf->initial_position_;
										if (h < 1e-9)
											h = (pos_v1 - pos_pf).norm();
										h = std::min(h, (pos_v1 - pos_pf).norm());
										n.push_front(pf);
										return true;
									}
									return true;
								});
								Particule_SPH_MR* p3 = value<Particule_SPH_MR*>(m, particule_volume_.get(), w2);
								Vec3 pos_v3 = p3->initial_position_;
								if (h < 1e-9)
									h = (pos_v1 - pos_v3).norm();
								h = std::min(h, (pos_v1 - pos_v3).norm());
								n.push_front(p3);
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
			for (Particule_SPH_MR* w : n)
			{
				Vec3 pos_v2 = w->initial_position_;
				sk += w->initial_volume_* Kernel_W((pos_v1 - pos_v2).norm(), h);
			}
			p.shepard_filter_ = 1.0f / sk;
			return true;
		});
	}

	void init_solver(MRMAP& m, Attribute<Vec3>* pos)
	{
		switch (particule_type)
		{
		case VOLUME_PARTICULE:
			compute_neighborhood_Volume(m, pos);
			break;
		case VOLUME_FACE_PARTICULE:
			compute_neighborhood_Volume_Face(m, pos);
			break;
		}
		compute_corrected_matrix();
	}

	void update_topo(const MRMAP& m, const std::vector<Vertex>&)
	{
	}
	
	void update_topo(const MRMAP& old_view,const MRMAP& new_view, const std::vector<Volume>& new_coarse, const std::vector<Volume>& new_fine)
	{
		std::unordered_set<Particule_SPH_MR*> need_upate_particules;
		for(Volume v : new_coarse){
			Particule_SPH_MR* p = value<Particule_SPH_MR*>(new_view, particule_volume_.get(), v);
			p->neighborhood_.assign(p->child_.front()->neighborhood_.begin(),p->child_.front()->neighborhood_.end());
			p->neighborhood_.remove_if([&p](Particule_SPH_MR* n)->bool{
				return std::count(p->child_.begin(),p->child_.end(),n) > 0;
			});
			for(Particule_SPH_MR* p2 : p->neighborhood_){
				p2->neighborhood_.push_front(p);
				p2->neighborhood_.remove_if([&p](Particule_SPH_MR* n)->bool{
					return std::count(p->child_.begin(),p->child_.end(),n) > 0;
				});
				std::remove_if(particules_.begin(),particules_.end(),[&p](Particule_SPH_MR* n)->bool{
					return std::count(p->child_.begin(),p->child_.end(),n) > 0;
				});
				need_upate_particules.insert(p2);
			}
			need_upate_particules.insert(p);
			p->neighborhood_.push_front(p);
			Vec3 pos = Vec3::Zero();
			Vec3 speed = Vec3::Zero();
			for(Particule_SPH_MR* c : p->child_){
				pos += c->masse_*c->current_position_;
				speed += c->masse_*c->speed_;
			}
			p->current_position_ = pos/p->masse_;
			p->speed_ = speed/p->masse_;
						particules_.push_back(p);
		}
		for(Volume v : new_fine){
			Particule_SPH_MR* p = value<Particule_SPH_MR*>(new_view, particule_volume_.get(), v);
			Particule_SPH_MR* p_old = value<Particule_SPH_MR*>(old_view, particule_volume_.get(), v);
			p->neighborhood_.assign(p_old->neighborhood_.begin(),p_old->neighborhood_.end());
			p->neighborhood_.remove(p_old);
			std::remove(particules_.begin(),particules_.end(),p_old);
			for(Particule_SPH_MR* p2 : p_old->neighborhood_){
				p2->neighborhood_.push_front(p);
				p2->neighborhood_.remove(p_old);
				need_upate_particules.insert(p2);
				for(auto p3:p2->neighborhood_){
					std::cout << p3 << std::endl;
				}
			}
			for(Particule_SPH_MR* p2 : p_old->child_){
				p->neighborhood_.push_front(p2);
				
			}
			need_upate_particules.insert(p);
			p->speed_ = p_old->speed_;
			Vec3 dir_initial_pos = p->initial_position_-p_old->initial_position_;
			Vec3 dir_current_pos = p_old->deformation_gradient_*dir_initial_pos;
			p->current_position_ = p_old->current_position_+dir_current_pos;
						particules_.push_back(p);
		}
		for (Particule_SPH_MR* p : need_upate_particules)
		{
			p->corrected_matrix_ = Corrected_matrix(*p, p->h_);
		}
	}
	
	void propagate_particule(){
		std::vector<Particule_SPH_MR*> fifo_particule;
		for (Particule_SPH_MR* p : particules_)
		{
			fifo_particule.push_back(p);
		};
		for(Particule_SPH_MR* p:fifo_particule){
			if(p->child_.empty())
				continue;
			for(Particule_SPH_MR* c:p->child_){
				Vec3 dir_initial_pos = c->initial_position_-p->initial_position_;
				Vec3 dir_current_pos = p->deformation_gradient_*dir_initial_pos;
				c->current_position_ = p->current_position_+dir_current_pos;
				c->deformation_gradient_ = p->deformation_gradient_;
				if(!c->child_.empty()){
					fifo_particule.push_back(c);
				}
			}
		}
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
	
	
	void init_particule_MR(MRMAP& m, Attribute<Vec3>* pos){
		
		MRMAP m2(m.m_);
		m2.current_level_ = m2.maximum_level_;
		MRMAP m3(m.m_);
		m3.current_level_ = std::max(0u,m3.maximum_level_-1);
		
		particule_volume_ =
			add_attribute<Particule_SPH_MR*, Volume>(m2, "SPH_particule_constraint_solver_particule_volume_" + id);
		
		initial_volume_ = add_attribute<double, Volume>(m2, "SPH_particule_constraint_solver_initial_volume_" + id);
		initial_centroid_volume_ =
			add_attribute<Vec3, Volume>(m2, "SPH_particule_constraint_solver_initial_centroid_volume_" + id);
		
		geometry::compute_centroid<Vec3, Volume>(m2, pos, initial_centroid_volume_.get());
		geometry::compute_volume(m2, pos, initial_volume_.get());
		
		foreach_cell(m2, [&](Volume v) -> bool {
			Particule_SPH_MR* p = new Particule_SPH_MR(value<Vec3>(m2, initial_centroid_volume_.get(), v), DENSITY_SPH*value<double>(m2, initial_volume_.get(), v));
			p->initial_volume_ = value<double>(m2, initial_volume_.get(), v);
			value<Particule_SPH_MR*>(m2, particule_volume_.get(), v) = p;
			return true;
		});
		CellMarker<MRMAP, Volume> marker(m2);
		while(m2.current_level_> 0){
			foreach_cell(m3, [&](Volume v) -> bool {
				if(m3.volume_level(v.dart) == m2.volume_level(v.dart))
					return true;
				Particule_SPH_MR* p = new Particule_SPH_MR(Vec3(0,0,0),0);
				value<Particule_SPH_MR*>(m3, particule_volume_.get(), v) = p;
				foreach_dart_of_orbit(m3,v,[&](Dart d)->bool{
					if(marker.is_marked(Volume(d)))
						return true;
					marker.mark(Volume(d));
					p->child_.push_front(value<Particule_SPH_MR*>(m2, particule_volume_.get(), Volume(d)));
					return true;
				});
				double masse =0;
				double volume = 0;
				Vec3 pos = Vec3::Zero();
				for(auto c : p->child_){
					masse+=c->masse_;
					volume+=c->initial_volume_;
					pos +=c->masse_ * c->initial_position_;
				}
				p->initial_position_ = pos/masse;
				p->current_position_ = pos/masse;
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
	
	void SPH_skinning(const MRMAP& m, Attribute<Vec3>* pos){
			parallel_foreach_cell(m,[&](Vertex v){
				Particule_SPH_MR* p = value<Particule_SPH_MR*>(m, particule_vertex_.get(), v);
				Vec3 new_pos = Vec3::Zero();
				for(Particule_SPH_MR* p2 : p->neighborhood_ ){
					new_pos += p2->initial_volume_
							*(p2->deformation_gradient_*(p->initial_position_-p2->initial_position_)+p2->current_position_)
							*Kernel_W((p->initial_position_ - p2->initial_position_).norm(), p->h_);
				}
				value<Vec3>(m, pos, v) = p->shepard_filter_*new_pos;
				return true;
			});
	}
	
	void particule_integration(double time_step){
		for (auto p : particules_)
		{
			p->RK_coeff[8] = p->force_;
		}
		compute_force_particules();
		for (auto p : particules_)
		{
			if(p->is_fixed){
				continue;
			}
			p->force_ -= 0.005*p->speed_;
			// k1
			p->RK_coeff[0] =
				(time_step * p->force_ / p->masse_);
			// j1
			p->RK_coeff[1] =
				time_step * (p->speed_ + p->RK_coeff[0]);
	
			p->current_position_ += p->RK_coeff[1] / 2.0f;
	
			p->force_ = p->RK_coeff[8];
		}
		compute_force_particules();
		for (auto p : particules_)
		{
			if(p->is_fixed){
				continue;
			}
			p->force_ -= 0.005*p->speed_;
			// k2
			p->RK_coeff[2] =
				(time_step * p->force_ / p->masse_);
			// j2
			p->RK_coeff[3] =
				time_step *
				(p->speed_ + p->RK_coeff[2] / 2.0f);
	
			p->current_position_ -= p->RK_coeff[1] / 2.0f;
	
			p->current_position_ += p->RK_coeff[3] / 2.0f;
	
			p->force_ = p->RK_coeff[8];
		}
		compute_force_particules();
		for (auto p : particules_)
		{
			if(p->is_fixed){
				continue;
			}
			p->force_ -= 0.005*p->speed_;
			// k3
			p->RK_coeff[4] =
				(time_step * p->force_ / p->masse_);
			// j3
			p->RK_coeff[5] =
				time_step *
				( p->speed_ + p->RK_coeff[4] / 2.0f);
	
			p->current_position_ -= p->RK_coeff[3] / 2.0f;
	
			p->current_position_ += p->RK_coeff[5];
	
			p->force_ = p->RK_coeff[8];
		}
		compute_force_particules();
		for (auto p : particules_)
		{
			if(p->is_fixed){
				continue;
			}
			p->force_ -= 0.005*p->speed_;
			// k4
			p->RK_coeff[6] =
				(time_step * p->force_ / p->masse_);
			// j4
			p->RK_coeff[7] =
				time_step * (p->speed_ + p->RK_coeff[6]);
	
			p->current_position_ -= p->RK_coeff[5];
	
			Vec3 k1 = p->RK_coeff[0];
			Vec3 k2 = p->RK_coeff[2];
			Vec3 k3 = p->RK_coeff[4];
			Vec3 k4 = p->RK_coeff[6];
	
			Vec3 j1 = p->RK_coeff[1];
			Vec3 j2 = p->RK_coeff[3];
			Vec3 j3 = p->RK_coeff[5];
			Vec3 j4 = p->RK_coeff[7];
	
			Vec3 diff_pos = 1. / 6. * (j1 + (2 * j2) + (2 * j3) + j4);
	
			Vec3 diff_speed = 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
	
			if (diff_pos.norm() < 1.0e-10)
				diff_pos = Vec3(0, 0, 0);
	
			if (diff_speed.norm() < 1.0e-10)
				diff_speed = Vec3(0, 0, 0);
	
			p->current_position_ += diff_pos;
	
			p->speed_ = p->speed_ + diff_speed;
			p->force_ = Vec3(0, 0, 0);
		}
	}
	
	template <typename FUNC>
	void set_particule_fixed(const FUNC& f){
		for(Particule_SPH_MR* p:particules_){
			p->is_fixed = f(*p);
		}
	}
	
	template <typename FUNC>
	void set_particule_forces(const FUNC& f){
		for(Particule_SPH_MR* p:particules_){
			p->force_ = f(*p);
		}
	}
	
	void test_particule_kernel(){
		
		for(Particule_SPH_MR* p:particules_){
			using Real = double;
			float eps = 1.0e-4f;
			Vec3 xi = p->initial_position_;
			Real sum = 0.0;
			Vec3 sumV = Vec3::Zero();
			bool positive = true;
			Real V =p->initial_volume_;
			std::cout << V << std::endl;
			for(Particule_SPH_MR* p2:p->neighborhood_){
				const Vec3 xj = p2->initial_position_;
				const Real W = Kernel_W((xi - xj).norm(), p->h_);
				sum += W * V;
				sumV += gradient(xi - xj, (xi - xj).norm(), p->h_) * V;
				if (W < -eps)
					positive = false;
			}
			if (fabs(sum - 1.0) < eps)
			{
				std::cout << "Kernel OK" << std::endl;
			}else{
				std::cout << sum << std::endl;
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
	}
 
	void solve_constraint(const MRMAP& m, Attribute<Vec3>* pos, Attribute<Vec3>* , double time_step) override
	{
		particule_integration(time_step);
		//SPH_skinning(m, pos);
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_SPH_PARTICULE_SPH_MULTIRESOLUTION_H_
