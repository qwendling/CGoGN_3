#ifndef CGOGN_SIMULATION_XPBD_XPBD_MULTIRESOLUTION_H_
#define CGOGN_SIMULATION_XPBD_XPBD_MULTIRESOLUTION_H_
#include <algorithm>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/cmap/EMR_Map3_Adaptative.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>

#define POISSON_RATIO 0.35
#define YOUNG_MODULUS 1e5

#define LAME_MU (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define LAME_LAMBDA ((YOUNG_MODULUS * POISSON_RATIO) / ((1 + POISSON_RATIO) * (1 - 2 * POISSON_RATIO)))

#define SHEAR_MODULUS (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define BULK_MODULUS (YOUNG_MODULUS / (3 * (1 - 2 * POISSON_RATIO)))

#define NUM_SUBSTEP 5
#define DENSITY 10

#define EPS 1e-12

namespace cgogn
{
namespace simulation
{
class XPBD_Multiresolution
{

	using Self = XPBD_Multiresolution;
	using MAP = EMR_Map3_Adaptative;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
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
	// Stable Values
	std::shared_ptr<Attribute<Vec3>> init_pos_;
	std::shared_ptr<Attribute<Vec3>> init_cm_;
	std::shared_ptr<Attribute<double>> init_volume_;
	std::shared_ptr<Attribute<double>> masse_;
	std::shared_ptr<Attribute<Mat3d>> inv_Q_;
	std::shared_ptr<Attribute<bool>> fixed_vertex;
	std::shared_ptr<Attribute<std::vector<Vertex>>> inc_vertices_;
	std::shared_ptr<Attribute<double>> Det_F_Volume_;
	std::shared_ptr<Attribute<Mat3d>> F_;
	std::shared_ptr<Attribute<double>> s_;
	// Integration Values

	std::shared_ptr<Attribute<Vec3>> centroid_;
	std::shared_ptr<Attribute<double>> distance_plan_;
	std::shared_ptr<Attribute<Vec3>> pos_;
	std::shared_ptr<Attribute<Vec3>> pos_prev_;
	std::shared_ptr<Attribute<Vec3>> speed_;
	std::shared_ptr<Attribute<Vec3>> f_ext_;
	std::shared_ptr<Attribute<Vec3>> Grad_C_i_;
	std::shared_ptr<Attribute<Vec3>> Grad_C2_i_;
	std::shared_ptr<Attribute<double>> error_volume_;

	tree_volume* hierarchy_;
	std::shared_ptr<Attribute<tree_volume*>> hierarchy_node_;
	uint32 nb_volume_current;
	uint32 nb_volume_init;

	XPBD_Multiresolution()
		: init_pos_(nullptr), init_cm_(nullptr), masse_(nullptr), inv_Q_(nullptr), inc_vertices_(nullptr),
		  Det_F_Volume_(nullptr), F_(nullptr), s_(nullptr), centroid_(nullptr), pos_(nullptr), pos_prev_(nullptr),
		  speed_(nullptr), f_ext_(nullptr), Grad_C_i_(nullptr), error_volume_(nullptr), hierarchy_(nullptr),
		  hierarchy_node_(nullptr)
	{
	}

	void init_solver(MAP& m, std::shared_ptr<Attribute<Vec3>> pos);

	void activate_volume(MAP& m, std::vector<Volume>& list_Volumes);
	void remove_volume(MAP& m, std::vector<Volume>& list_Volumes);
	void activate_remove_volume(MAP& m, std::vector<Volume>& list_activate, std::vector<Volume>& list_remove);
	void update_topo(MAP& m);

	void constraint_Neo_Hookean_H(MAP& m, Volume v, double h);
	void constraint_Neo_Hookean_D(MAP& m, Volume v, double h);
	void constraint_Zero_Energy(MAP& m, Volume v, double);
	void applyDamping(MAP& m, Volume v, double damping_coeff, double time_step);

	void compute_error(MAP& m, std::vector<Volume>& volume_activate, std::vector<Volume>& volume_disable);

	void solve_surface(MAP& m, MAP& geom, Volume v);

	void solver(MAP& m, MAP* geom, double timestep, bool allow_modif_topo = true);

	void compute_error_point(MAP& m, const Vec3& p, double quotat);

	template <typename FUNC>
	void compute_contact(MAP& m, const FUNC& f_contact)
	{
		foreach_cell(m, [&](Volume v) -> bool {
			tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
			t->have_contact = false;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
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

	template <typename FUNC>
	void apply_cut(MAP& m, Vec3 dir_plan, double w, const FUNC& callback_vertices);
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_XPBD_XPBD_MULTIRESOLUTION_H_
