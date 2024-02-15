#ifndef CGOGN_SIMULATION_XPBD_XPBD_MULTIRESOLUTION_H_
#define CGOGN_SIMULATION_XPBD_XPBD_MULTIRESOLUTION_H_
#include <algorithm>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/cmap/EMR_Map3_Adaptative.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>

#define POISSON_RATIO 0.45
#define YOUNG_MODULUS 1e7

#define LAME_MU (YOUNG_MODULUS / (2.0 * (1.0 + POISSON_RATIO)))
#define LAME_LAMBDA ((YOUNG_MODULUS * POISSON_RATIO) / ((1.0 + POISSON_RATIO) * (1.0 - 2.0 * POISSON_RATIO)))

#define SHEAR_MODULUS (YOUNG_MODULUS / (2.0 * (1.0 + POISSON_RATIO)))
#define BULK_MODULUS (YOUNG_MODULUS / (3.0 * (1.0 - 2.0 * POISSON_RATIO)))

#define NUM_SUBSTEP 5
#define DENSITY_SPH 10
#define NORMALIZE_TERM (21 / (2 * M_PI))

#define EPS 1e-12

namespace cgogn
{
namespace simulation
{
class XPBD_SPH_Multiresolution
{
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	struct Particule_SPH_MR
	{

		Vec3 initial_position_;
		Vec3 current_position_;
		Vec3 previous_position_;
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
		double CH;
		double CD;
		Vec3 dFdX;

		Particule_SPH_MR(Vec3 pos, double masse)
			: initial_position_(pos), current_position_(pos), force_(0, 0, 0), speed_(0, 0, 0), masse_(masse),
			  is_fixed(false)
		{
			rotation_ = std::move(Mat3d::Identity());
			deformation_gradient_ = std::move(Mat3d::Identity());
		}
	};

	using Self = XPBD_SPH_Multiresolution;
	using MAP = EMR_Map3_Adaptative;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Volume = typename mesh_traits<MAP>::Volume;
	using Face = typename mesh_traits<MAP>::Face;

	enum tree_volume_node
	{
		ROOT,
		NONE,
		CURRENT,
		COARSE,
		TOPOLOGY
	};

	struct tree_volume
	{
		int id;
		tree_volume* fils;
		tree_volume* pere;
		tree_volume* frere;
		Dart volume_dart;
		Mat3d F_;
		Vec3 init_cm_;
		Vec3 cm_;
		Vec3 v_cm_;
		double error;

		bool is_topo;
		bool have_contact;

		tree_volume_node type;
		uint32 clock;
		tree_volume() : fils(nullptr), pere(nullptr), frere(nullptr), is_topo(false), have_contact(false), clock(0)
		{
		}

		template <typename FUNC>
		void for_each_child(const FUNC& fn)
		{
			static_assert(is_func_parameter_same<FUNC, tree_volume*>::value,
						  "Given function must receive tree_volume* as parameter");
			static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
			tree_volume* it = fils;
			while (it != nullptr)
			{
				if (!(fn(it)))
					return;
				it = it->frere;
			}
		}

		void print()
		{
			if (!fils)
			{
				return;
			}

			switch (type)
			{
			case ROOT:
				std::cout << "ROOT_" + std::to_string(id);
				break;
			case NONE:
				std::cout << "NONE_" + std::to_string(id);
				break;
			case CURRENT:
				std::cout << "CURRENT_" + std::to_string(id);
				break;
			case COARSE:
				std::cout << "COARSE_" + std::to_string(id);
				break;
			case TOPOLOGY:
				std::cout << "TOPOLOGY_" + std::to_string(id);
				break;
			}

			std::cout << "->";
			tree_volume* it = fils->frere;
			while (it != nullptr)
			{
				switch (it->type)
				{
				case NONE:
					std::cout << "NONE_";
					break;
				case CURRENT:
					std::cout << "CURRENT_";
					break;
				case COARSE:
					std::cout << "COARSE_";
					break;
				case TOPOLOGY:
					std::cout << "TOPOLOGY_";
					break;
				}
				std::cout << std::to_string(it->id) + ",";
				it = it->frere;
			}
			switch (fils->type)
			{
			case NONE:
				std::cout << "NONE_";
				break;
			case CURRENT:
				std::cout << "CURRENT_";
				break;
			case COARSE:
				std::cout << "COARSE_";
				break;
			case TOPOLOGY:
				std::cout << "TOPOLOGY_";
				break;
			}
			std::cout << std::to_string(fils->id) + ";" << std::endl;
			for_each_child([&](tree_volume* c) -> bool {
				c->print();
				return true;
			});
		}
	};

public:
	static inline int nb_solver = 0;
	int id;
	std::vector<Particule_SPH_MR*> particules_;

	// Stable Values
	std::shared_ptr<Attribute<Particule_SPH_MR*>> particule_vertex_;
	std::shared_ptr<Attribute<Particule_SPH_MR*>> particule_volume_;
	// Integration Values

	std::shared_ptr<Attribute<double>> initial_volume_;
	std::shared_ptr<Attribute<Vec3>> initial_centroid_volume_;

	tree_volume* hierarchy_;
	std::shared_ptr<Attribute<tree_volume*>> hierarchy_node_;
	uint32 nb_volume_current;
	uint32 nb_volume_init;

	enum type_particule
	{
		VOLUME_PARTICULE,
		VERTEX_PARTICULE
	};
	type_particule particule_type;

	XPBD_SPH_Multiresolution()
		: id(nb_solver++), particule_vertex_(nullptr), particule_volume_(nullptr), initial_volume_(nullptr),
		  initial_centroid_volume_(nullptr), hierarchy_(nullptr), hierarchy_node_(nullptr)
	{
	}

	void init_solver(MAP& m, std::shared_ptr<Attribute<Vec3>> pos);
	void init_particule_MR(MAP& m, Attribute<Vec3>* pos);
	void compute_neighborhood_Volume(MAP& m, Attribute<Vec3>* pos);

	double Kernel_W(double dist, double h) const;

	Vec3 gradient(Vec3 xij, double dist, double h) const;

	Mat3d Corrected_matrix(Particule_SPH_MR& p, double h) const;

	void compute_corrected_matrix();

	Vec3 Corrected_gradient(Particule_SPH_MR& vi, Particule_SPH_MR& vj, double h) const;

	void compute_deformation_gradient_particule(Particule_SPH_MR* p);
	void compute_deformation_gradient_particules();

	template <typename FUNC>
	void set_particule_fixed(const FUNC& f)
	{
		for (Particule_SPH_MR* p : particules_)
		{
			p->is_fixed = f(*p);
		}
	}

	template <typename FUNC>
	void set_particule_forces(const FUNC& f)
	{
		for (Particule_SPH_MR* p : particules_)
		{
			p->force_ = f(*p);
		}
	}

	void activate_volume(MAP& m, std::vector<Volume>& list_Volumes);
	void remove_volume(MAP& m, std::vector<Volume>& list_Volumes);
	void activate_remove_volume(MAP& m, std::vector<Volume>& list_activate, std::vector<Volume>& list_remove);
	void update_topo(MAP& m);

	void constraint_Neo_Hookean_H(MAP& m, Particule_SPH_MR& p, double h);
	void constraint_Neo_Hookean_D(MAP& m, Particule_SPH_MR& p, double h);
	void constraint_Zero_Energy(MAP& m, Particule_SPH_MR& p, double);
	void applyDamping(MAP& m, Volume v, double damping_coeff, double time_step);

	void compute_error(MAP& m, std::vector<Volume>& volume_activate, std::vector<Volume>& volume_disable);

	void solve_surface(MAP& m, MAP& geom, Volume v);

	void solver(MAP& m, MAP* geom, double timestep, bool allow_modif_topo = true);

	void compute_error_point(MAP& m, const Vec3& p, double quotat);

	template <typename FUNC>
	void compute_contact(MAP& m, MAP& geom, const FUNC& f_contact, uint32 max_level = 100u)
	{
		foreach_cell(m, [&](Volume v) -> bool {
			tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
			t->have_contact = false;
			if (m.volume_level(v.dart) >= max_level)
				return true;
			foreach_incident_vertex(geom, v, [&](Vertex w) -> bool {
				if (f_contact(w))
					t->have_contact = true;
				return !t->have_contact;
			});
			return true;
		});
		std::vector<Volume> volume_disable;
		std::vector<Volume> volume_activate;
		foreach_cell(m, [&](Volume v) -> bool {
			tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
			if (t->type == CURRENT && t->fils != nullptr)
			{
				if (t->have_contact)
				{
					t->have_contact = false;
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
				t->pere->have_contact = false;
				t->pere->for_each_child([&](tree_volume* c) -> bool {
					t->pere->have_contact = c->have_contact;
					return !c->have_contact;
				});
				if (!t->pere->have_contact)
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
		activate_remove_volume(m, volume_activate, volume_disable);
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_XPBD_XPBD_MULTIRESOLUTION_H_
