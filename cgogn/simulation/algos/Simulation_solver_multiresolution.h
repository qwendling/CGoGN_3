#ifndef CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
#define CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>
#include <cgogn/simulation/algos/Simulation_solver.h>
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_constraint.h>
#include <forward_list>
#include <string>

#define QUOTA_VOLUME 1000
#define PROPORTION_QUOTA 0.1f
#if 0
#define MODIF_MAX PROPORTION_QUOTA* QUOTA_VOLUME
#else
#define MODIF_MAX 1
#endif

namespace cgogn
{
namespace simulation
{
template <typename MR_MAP>
class Simulation_solver_multiresolution : public Simulation_solver<MR_MAP>
{
	using Self = Simulation_solver_multiresolution<MR_MAP>;
	using Inherit = Simulation_solver<MR_MAP>;
	template <typename T>
	using Attribute = typename mesh_traits<MR_MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Vertex = typename mesh_traits<MR_MAP>::Vertex;
	using Face = typename mesh_traits<MR_MAP>::Face;
	using Volume = typename mesh_traits<MR_MAP>::Volume;

	struct Module
	{
		Simulation_solver_multiresolution* ssm_;
	};

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
		/*bool is_current;
		bool is_coarse;*/
		tree_volume_node type;
		int clock;
		tree_volume() : fils(nullptr), pere(nullptr), frere(nullptr), clock(0)
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
	Propagation_Constraint<MR_MAP>* pc_;
	Simulation_constraint<MR_MAP>* sc_;
	std::shared_ptr<Simulation_constraint<MR_MAP>> sc_fine_;
	std::shared_ptr<Simulation_constraint<MR_MAP>> sc_coarse_;
	std::shared_ptr<Attribute<std::array<Vertex, 4>>> parents_;
	std::shared_ptr<Attribute<Vec3>> relative_pos_;
	MR_MAP* mecanical_mesh_;
	MR_MAP* fine_meca_mesh_;
	MR_MAP* coarse_meca_mesh_;
	MR_MAP* topology_;

	tree_volume* hierarchy_;

	std::forward_list<tree_volume*> list_volume_coarse_;
	std::forward_list<tree_volume*> list_volume_current_;

	std::shared_ptr<Attribute<tree_volume*>> hierarchy_node_;

	std::shared_ptr<Attribute<Vec3>> pos_coarse_;
	std::shared_ptr<Attribute<Vec3>> pos_current_;
	std::shared_ptr<Attribute<Vec3>> forces_coarse_;
	std::shared_ptr<Attribute<Vec3>> forces_current_;

	std::shared_ptr<Attribute<double>> diff_volume_current_fine_;
	std::shared_ptr<Attribute<double>> diff_volume_coarse_current_;
	std::shared_ptr<Attribute<double>> volume_coarse_;

	std::shared_ptr<Attribute<double>> area_face_;
	Vec3 gravity_;

	CellCache<MR_MAP>* cache_current_vol_;
	int clock;

	Simulation_solver_multiresolution()
		: Simulation_solver<MR_MAP>(), pc_(nullptr), parents_(nullptr), relative_pos_(nullptr), hierarchy_(nullptr),
		  gravity_(0, 0, 0), cache_current_vol_(nullptr), clock(1)
	{
	}

	void init_cell_cache(CellCache<MR_MAP>& cc, std::forward_list<tree_volume*>& lv)
	{
		cc.template clear<Volume>();
		for (auto t : lv)
		{
			cc.add(Volume(t->volume_dart));
		}
	}

	void create_coarse_view()
	{
		list_volume_coarse_.clear();
		coarse_meca_mesh_ = mecanical_mesh_->get_copy();
		std::vector<Volume> volume_coarse;
		foreach_cell(*coarse_meca_mesh_, [&](Volume v) -> bool {
			tree_volume* t = value<tree_volume*>(*coarse_meca_mesh_, hierarchy_node_, v);
			if (t->pere->clock != clock)
			{
				t->pere->clock = clock;
				bool can_be_coarse = true;
				t->pere->for_each_child([&](tree_volume* c) -> bool { return can_be_coarse = (c->type == CURRENT); });
				if (can_be_coarse)
				{
					t->pere->type = COARSE;
					volume_coarse.push_back(Volume(t->pere->volume_dart));
					list_volume_coarse_.push_front(t->pere);
				}
			}
			return true;
		});
		clock++;
		for (Volume v : volume_coarse)
		{
			std::cout << "disable volume lvl : " << coarse_meca_mesh_->volume_level(v.dart) << std::endl;
			coarse_meca_mesh_->disable_volume_subdivision(v, true);
		}
	}

	void create_fine_view()
	{
		fine_meca_mesh_ = mecanical_mesh_->get_copy();
		foreach_cell(*mecanical_mesh_, [&](Volume v) -> bool {
			fine_meca_mesh_->activate_volume_subdivision(v);
			return true;
		});
	}

	void init_solver(MR_MAP& m, Simulation_constraint<MR_MAP>* sc, Attribute<Vec3>* pos,
					 Propagation_Constraint<MR_MAP>* pc = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& speed = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& forces = nullptr)
	{
		using MR_Base = typename MR_MAP::Inherit;
		Inherit::init_solver(m, sc, speed, forces);
		pc_ = pc;
		mecanical_mesh_ = &m;
		fine_meca_mesh_ = m.get_copy();
		coarse_meca_mesh_ = m.get_copy();
		topology_ = new MR_MAP(m.m_);

		std::vector<tree_volume*> tmp_watcher;
		std::forward_list<tree_volume*> list_volume_current_tmp;

		// Construction de la hierarchie de volume

		MR_Base tmp(m);
		tmp.current_level_ = m.current_level_;
		std::function<void(tree_volume*, MR_Base&)> progress_tree;

		int nb_node = 0;

		progress_tree = [&](tree_volume* p, MR_Base& cph) -> void {
			if (cph.current_level_ == cph.maximum_level_)
				return;
			foreach_incident_vertex(cph, Volume(p->volume_dart), [&](Vertex v) -> bool {
				cph.current_level_++;

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
			value<tree_volume*>(tmp, hierarchy_node_, v) = t;
			progress_tree(t, tmp2);
			return true;
		});
		// init_cell_cache(*cache_current_vol_, list_volume_current_);

		// Construction des vues

		list_volume_current_.clear();

		foreach_cell(m, [&](Volume v) -> bool {
			tree_volume* t = value<tree_volume*>(m, hierarchy_node_, v);
			t->type = CURRENT;
			std::cout << "current volume lvl : " << m.volume_level(v.dart) << std::endl;
			if (t->fils)
				list_volume_current_.push_front(t);
			return true;
		});

		create_coarse_view();
		create_fine_view();

		hierarchy_->print();

		fine_meca_mesh_->current_level_ = mecanical_mesh_->current_level_;
		reset_forces(*fine_meca_mesh_);

		sc_ = sc;
		sc_->init_solver(*mecanical_mesh_, pos);
		sc_fine_ = std::shared_ptr<Simulation_constraint<MR_MAP>>(sc->get_new_ptr());
		sc_fine_->init_solver(*fine_meca_mesh_, pos);
		sc_coarse_ = std::shared_ptr<Simulation_constraint<MR_MAP>>(sc->get_new_ptr());
		sc_coarse_->init_solver(*coarse_meca_mesh_, pos);

		pos_coarse_ = get_attribute<Vec3, Vertex>(m, "Solver_multiresolution_pos_coarse");
		if (pos_coarse_ == nullptr)
			pos_coarse_ = add_attribute<Vec3, Vertex>(m, "Solver_multiresolution_pos_coarse");

		pos_current_ = get_attribute<Vec3, Vertex>(m, "Solver_multiresolution_pos_current");
		if (pos_current_ == nullptr)
			pos_current_ = add_attribute<Vec3, Vertex>(m, "Solver_multiresolution_pos_current");

		forces_current_ = get_attribute<Vec3, Vertex>(m, "Solver_multiresolution_forces_current");
		if (forces_current_ == nullptr)
			forces_current_ = add_attribute<Vec3, Vertex>(m, "Solver_multiresolution_forces_current");

		forces_coarse_ = get_attribute<Vec3, Vertex>(m, "Solver_multiresolution_forces_coarse");
		if (forces_coarse_ == nullptr)
			forces_coarse_ = add_attribute<Vec3, Vertex>(m, "Solver_multiresolution_forces_coarse");

		diff_volume_current_fine_ = get_attribute<double, Volume>(m, "Solver_multiresolution_diff_volume_current_fine");
		if (diff_volume_current_fine_ == nullptr)
			diff_volume_current_fine_ =
				add_attribute<double, Volume>(m, "Solver_multiresolution_diff_volume_current_fine");
		volume_coarse_ = get_attribute<double, Volume>(m, "Solver_multiresolution_volume_coarse");
		if (volume_coarse_ == nullptr)
			volume_coarse_ = add_attribute<double, Volume>(m, "Solver_multiresolution_volume_coarse");
		diff_volume_coarse_current_ =
			get_attribute<double, Volume>(m, "Solver_multiresolution_diff_volume_coarse_current");
		if (diff_volume_coarse_current_ == nullptr)
			diff_volume_coarse_current_ =
				add_attribute<double, Volume>(m, "Solver_multiresolution_diff_volume_coarse_current");

		area_face_ = get_attribute<double, Face>(m, "Solver_multiresolution_area_face");
		if (area_face_ == nullptr)
			area_face_ = add_attribute<double, Face>(m, "Solver_multiresolution_area_face");
	}

	void reset_forces(MR_MAP& m)
	{
		parallel_foreach_cell(m.m_, [&](Vertex v) -> bool {
			value<Vec3>(m, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			value<Vec3>(m, this->speed_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
	}

	int get_size_list(const std::forward_list<tree_volume*>& l)
	{
		int i = 0;
		for (auto t : l)
			i++;
		return i;
	}

	bool update_topo(Attribute<Vec3>* vertex_position)
	{
		std::clock_t start;
		double duration;
		start = std::clock();
		if (list_volume_current_.empty() || list_volume_coarse_.empty())
		{
			return false;
		}
		list_volume_current_.sort([&](tree_volume* t1, tree_volume* t2) {
			double v1 = value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), Volume(t1->volume_dart));
			double v2 = value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), Volume(t2->volume_dart));
			return v1 > v2;
		});
		list_volume_coarse_.sort([&](tree_volume* t1, tree_volume* t2) {
			double v1 =
				value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), Volume(t1->volume_dart));
			double v2 =
				value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), Volume(t2->volume_dart));
			return v1 < v2;
		});

		std::forward_list<tree_volume*> list_new_current_;
		std::forward_list<tree_volume*> list_new_coarse_;
		std::vector<Volume> list_new_volume_coarse;
		std::vector<Volume> list_new_volume_current;
		std::vector<Volume> list_new_volume_fine;

		CellMarkerStore<MR_MAP, Vertex> marker(*fine_meca_mesh_);
		foreach_cell(*fine_meca_mesh_, [&](Vertex v) -> bool {
			marker.mark(v);
			return true;
		});
		tree_volume* min_coarse = list_volume_coarse_.front();
		tree_volume* max_fine = list_volume_current_.front();
		if (min_coarse == max_fine->pere)
		{
			list_volume_coarse_.pop_front();
			if (list_volume_coarse_.empty())
				return false;
			min_coarse = list_volume_coarse_.front();
		}
		double v1 =
			value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), Volume(max_fine->volume_dart));
		double v2 =
			value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), Volume(min_coarse->volume_dart));
		int nb_modif = 0;
		while (v2 / v1 < 0.7 && nb_modif < MODIF_MAX)
		{
			list_volume_current_.pop_front();
			list_volume_coarse_.pop_front();

			// Activation
			/*max_fine->is_current = false;
			max_fine->is_coarse = true;*/
			max_fine->type = COARSE;
			list_new_coarse_.push_front(max_fine);

			// if (max_fine->pere && max_fine->pere->is_coarse)
			if (max_fine->pere && max_fine->pere->type == COARSE)
			{
				list_volume_coarse_.remove(max_fine->pere);
				// max_fine->pere->is_coarse = false;
				max_fine->pere->type = NONE;
				coarse_meca_mesh_->activate_volume_subdivision(Volume(max_fine->volume_dart));
				max_fine->pere->for_each_child([&](tree_volume* c) -> bool {
					list_new_volume_coarse.push_back(Volume(c->volume_dart));
					return true;
				});
			}

			mecanical_mesh_->activate_volume_subdivision(Volume(max_fine->volume_dart));
			max_fine->for_each_child([&](tree_volume* c) -> bool {
				list_new_volume_current.push_back(Volume(c->volume_dart));
				if (c->fils)
				{
					c->for_each_child([&](tree_volume* cc) -> bool {
						list_new_volume_fine.push_back(Volume(cc->volume_dart));
						return true;
					});
					fine_meca_mesh_->activate_volume_subdivision(Volume(c->volume_dart));
					list_new_current_.push_front(c);
				}
				// c->is_current = true;
				c->type = CURRENT;
				return true;
			});
			// Disable
			// min_coarse->is_current = true;
			min_coarse->type = CURRENT;

			list_new_current_.push_front(min_coarse);

			mecanical_mesh_->disable_volume_subdivision(Volume(min_coarse->volume_dart), true);
			list_new_volume_current.push_back(Volume(min_coarse->volume_dart));
			min_coarse->for_each_child([&](tree_volume* c) -> bool {
				list_volume_current_.remove(c);
				// c->is_current = false;
				c->type = NONE;

				if (c->fils)
				{
					list_new_volume_fine.push_back(Volume(c->volume_dart));
					fine_meca_mesh_->disable_volume_subdivision(Volume(c->volume_dart), true);
				}
				return true;
			});

			if (min_coarse->pere)
			{

				bool can_be_coarse = true;
				min_coarse->pere->for_each_child([&](tree_volume* c) {
					// if (!c->is_current)
					if (c->type != CURRENT)
						can_be_coarse = false;
					return can_be_coarse;
				});
				if (can_be_coarse)
				{
					list_new_volume_coarse.push_back(Volume(min_coarse->volume_dart));
					coarse_meca_mesh_->disable_volume_subdivision(Volume(min_coarse->volume_dart), true);
					// min_coarse->pere->is_coarse = true;
					min_coarse->pere->type = COARSE;
					list_new_coarse_.push_front(min_coarse->pere);
				}
			}

			// end
			if (list_volume_coarse_.empty() || list_volume_current_.empty())
				break;
			min_coarse = list_volume_coarse_.front();
			max_fine = list_volume_current_.front();
			if (min_coarse == max_fine->pere)
			{
				list_volume_coarse_.pop_front();
				if (list_volume_coarse_.empty())
					break;
				min_coarse = list_volume_coarse_.front();
			}
			v1 = value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), Volume(max_fine->volume_dart));
			v2 = value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(),
							   Volume(min_coarse->volume_dart));
			nb_modif++;
		}
		for (auto t : list_new_coarse_)
			list_volume_coarse_.push_front(t);
		for (auto t : list_new_current_)
			list_volume_current_.push_front(t);

		if (!list_new_volume_current.empty() || !list_new_coarse_.empty())
		{
			CellMarker<MR_MAP, Vertex> coarse_marker(*coarse_meca_mesh_);
			std::vector<Vertex> list_new_vertices_coarse;
			uint tmp_cmp = 0;

			for (Volume v : list_new_volume_coarse)
			{
				foreach_incident_vertex(*coarse_meca_mesh_, v, [&](Vertex w) -> bool {
					if (!coarse_marker.is_marked(w))
					{
						coarse_marker.mark(w);
						list_new_vertices_coarse.push_back(w);
					}
					return true;
				});
				tmp_cmp++;
			}
			CellMarker<MR_MAP, Vertex> current_marker(*mecanical_mesh_);
			std::vector<Vertex> list_new_vertices_current;
			for (Volume v : list_new_volume_current)
			{
				foreach_incident_vertex(*mecanical_mesh_, v, [&](Vertex w) -> bool {
					if (!current_marker.is_marked(w))
					{
						current_marker.mark(w);
						list_new_vertices_current.push_back(w);
					}
					return true;
				});
			}
			CellMarker<MR_MAP, Vertex> fine_marker(*fine_meca_mesh_);
			std::vector<Vertex> list_new_vertices_fine;
			for (Volume v : list_new_volume_fine)
			{
				foreach_incident_vertex(*fine_meca_mesh_, v, [&](Vertex w) -> bool {
					if (!fine_marker.is_marked(w))
					{
						fine_marker.mark(w);
						list_new_vertices_fine.push_back(w);
					}
					return true;
				});
			}

			std::clock_t start_update;
			start_update = std::clock();

			if (list_new_vertices_coarse.size() > 0)
			{
				sc_coarse_->update_topo(*coarse_meca_mesh_, list_new_vertices_coarse);
			}

			start_update = std::clock();
			if (list_new_vertices_current.size() > 0)
			{
				sc_->update_topo(*mecanical_mesh_, list_new_vertices_current);
			}

			start_update = std::clock();
			if (list_new_vertices_fine.size() > 0)
			{
				sc_fine_->update_topo(*fine_meca_mesh_, list_new_vertices_fine);
			}

			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "time activation/disable + update topo simu : " << duration << std::endl;
			pc_->propagate(*mecanical_mesh_, *fine_meca_mesh_, vertex_position, this->forces_ext_.get(),
						   sc_fine_->masse_.get(), relative_pos_.get(), parents_.get(),
						   [&](Vertex v) -> bool { return !marker.is_marked(v); });
			pc_->propagate(*mecanical_mesh_, *fine_meca_mesh_, this->speed_.get(), this->forces_ext_.get(),
						   sc_fine_->masse_.get(), relative_pos_.get(), parents_.get(),
						   [&](Vertex v) -> bool { return !marker.is_marked(v); });
			return true;
		}
		return false;
	}

	void compute_error(Attribute<Vec3>* vertex_position, Attribute<double>* masse, double time_step)
	{

		if (pc_)
		{
			pc_->propagate(*coarse_meca_mesh_, *mecanical_mesh_, pos_coarse_.get(), this->forces_coarse_.get(), masse,
						   relative_pos_.get(), parents_.get(), time_step);
			pc_->propagate(*mecanical_mesh_, *fine_meca_mesh_, pos_current_.get(), this->forces_current_.get(),
						   sc_fine_->masse_.get(), relative_pos_.get(), parents_.get(), time_step);
		}

		std::clock_t start_volume;
		double duration_volume = 0;

		std::clock_t start;
		double duration;
		start = std::clock();

		bool volume_current_is_finish = false;
		auto volume_current = [&]() {
			auto fn = [&](tree_volume* t) -> bool {
				start_volume = std::clock();
				value<double>(*mecanical_mesh_, this->diff_volume_coarse_current_.get(), Volume(t->volume_dart)) =
					geometry::volume(*mecanical_mesh_, Volume(t->volume_dart), pos_coarse_.get());
				duration_volume += (std::clock() - start_volume) / (double)CLOCKS_PER_SEC;
				return true;
			};
			for (tree_volume* tp : list_volume_coarse_)
			{
				tp->for_each_child(fn);
			}
			volume_current_is_finish = true;
		};

		start = std::clock();
		volume_current();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "\033[1;32m time volume current : \033[0m" << duration << std::endl;

		bool volume_fine_is_finish = false;
		auto volume_fine = [&]() {
			auto fn = [&](tree_volume* t) -> bool {
				start_volume = std::clock();
				value<double>(*fine_meca_mesh_, this->diff_volume_current_fine_.get(), Volume(t->volume_dart)) =
					geometry::volume(*fine_meca_mesh_, Volume(t->volume_dart), pos_current_.get());
				duration_volume += (std::clock() - start_volume) / (double)CLOCKS_PER_SEC;
				return true;
			};
			for (tree_volume* tp : list_volume_current_)
			{
				tp->for_each_child(fn);
			}
			volume_fine_is_finish = true;
		};

		start = std::clock();
		volume_fine();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "\033[1;33m time volume fine : \033[0m" << duration << std::endl;

		bool error_coarse_current_is_finish = false;
		auto error_coarse_current = [&]() {
			auto fn = [&](tree_volume* t) -> bool {
				Volume v = Volume(t->volume_dart);
				start_volume = std::clock();
				double vol = geometry::volume(*mecanical_mesh_, v, vertex_position);
				duration_volume += (std::clock() - start_volume) / (double)CLOCKS_PER_SEC;
				double vol2 = value<double>(*mecanical_mesh_, this->diff_volume_coarse_current_.get(), v);
				double diff;
				diff = vol2 / vol;
				if (diff > 1)
					diff = 1 / diff;
				diff = 1 - diff;
				diff *= log(fabs(vol2 - vol) + 1);
				value<double>(*mecanical_mesh_, this->diff_volume_coarse_current_.get(), v) = diff;
				return true;
			};
			for (tree_volume* tp : list_volume_coarse_)
			{
				tp->for_each_child(fn);
			}
			for (tree_volume* t : list_volume_coarse_)
			{
				// if (t->is_current)
				if (t->type == CURRENT)
					continue;
				Volume v = Volume(t->volume_dart);
				value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), v) = 0;
				int i = 0;
				t->for_each_child([&](tree_volume* t2) -> bool {
					i++;
					Volume v2 = Volume(t2->volume_dart);
					value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), v) +=
						value<double>(*mecanical_mesh_, this->diff_volume_coarse_current_.get(), v2);
					return true;
				});
				if (i == 0)
				{
					value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), v) = 0;
					continue;
				}
				value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), v) /= i;
			}
			error_coarse_current_is_finish = true;
			// cv.notify_all();
		};

		start = std::clock();
		error_coarse_current();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "\033[1;34m time error coarse current : \033[0m" << duration << std::endl;

		bool error_current_fine_is_finish = false;
		auto error_current_fine = [&]() {
			auto fn2 = [&](tree_volume* t) -> bool {
				Volume v = Volume(t->volume_dart);
				start_volume = std::clock();
				double vol = geometry::volume(*fine_meca_mesh_, v, vertex_position);
				duration_volume += (std::clock() - start_volume) / (double)CLOCKS_PER_SEC;
				double vol2 = value<double>(*fine_meca_mesh_, this->diff_volume_current_fine_.get(), v);
				double diff;
				diff = vol2 / vol;
				if (diff > 1)
					diff = 1 / diff;
				diff = 1 - diff;
				diff *= log(fabs(vol2 - vol) + 1);
				value<double>(*fine_meca_mesh_, this->diff_volume_current_fine_.get(), v) = diff;
				return true;
			};
			for (tree_volume* tp : list_volume_current_)
			{
				tp->for_each_child(fn2);
			}

			auto fn = [&](tree_volume* t) -> bool {
				Volume v = Volume(t->volume_dart);
				value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), v) = 0;
				if (t->fils == nullptr)
					return true;
				int i = 0;
				t->for_each_child([&](tree_volume* t2) -> bool {
					i++;
					Volume v2 = Volume(t2->volume_dart);
					value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), v) +=
						value<double>(*fine_meca_mesh_, this->diff_volume_current_fine_.get(), v2);
					return true;
				});
				if (i == 0)
				{
					value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), v) = 0;
					return true;
				}
				value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), v) /= i;

				return true;
			};
			for (tree_volume* tp : list_volume_current_)
			{
				fn(tp);
			}
			error_current_fine_is_finish = true;
		};

		start = std::clock();
		error_current_fine();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "\033[1;35m time error current fine : \033[0m" << duration << std::endl;
		std::cout << "temps total compute volume : " << duration_volume << std::endl;
	}

	void add_new_vertices_in_simulation(MR_MAP& view, Attribute<Vec3>* vertex_position)
	{
		CellMarkerStore<MR_MAP, Vertex> marker(view);
		foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool { return true; });
		pc_->propagate(*mecanical_mesh_, view, vertex_position, nullptr, nullptr, relative_pos_.get(), parents_.get(),
					   [&](Vertex v) -> bool { return !marker.is_marked(v); });
		pc_->propagate(*mecanical_mesh_, *fine_meca_mesh_, this->speed_.get(), this->forces_ext_.get(), nullptr,
					   relative_pos_.get(), parents_.get(), [&](Vertex v) -> bool { return !marker.is_marked(v); });
	}

	void compute_time_step(MR_MAP& m_geom, Attribute<Vec3>* vertex_position, Attribute<double>* masse, double time_step,
						   bool& modif_topo)
	{
		static int nb_step = 0;
		std::cout << "step : " << nb_step++ << std::endl;
		if (!this->constraint_)
			return;

		std::clock_t start;
		double duration;
		start = std::clock();
		/////////////////////////////////////////////////
		////	         coarse_solve				/////
		/////////////////////////////////////////////////
		auto solve_coarse = [&]() {
			for (const Vertex& v : sc_->vertices_cache)
			{
				value<Vec3>(*mecanical_mesh_, this->forces_coarse_.get(), v) =
					value<Vec3>(*coarse_meca_mesh_, this->forces_ext_.get(), v);
			}

			sc_coarse_->solve_constraint(*coarse_meca_mesh_, vertex_position, this->forces_coarse_.get(), time_step);

			for (const Vertex& v : sc_coarse_->vertices_cache)
			{
				if (this->fixed_vertex && value<bool>(*coarse_meca_mesh_, this->fixed_vertex.get(), v))
				{
					value<Vec3>(*coarse_meca_mesh_, pos_coarse_, v) =
						value<Vec3>(*coarse_meca_mesh_, vertex_position, v);
					continue;
				}
				// compute speed
				Vec3 s = 0.995 * value<Vec3>(*coarse_meca_mesh_, this->speed_.get(), v) +
						 time_step * value<Vec3>(*coarse_meca_mesh_, this->forces_coarse_.get(), v) /
							 value<double>(*coarse_meca_mesh_, sc_coarse_->masse_, v);
				s += time_step * gravity_;
				value<Vec3>(*coarse_meca_mesh_, pos_coarse_, v) =
					value<Vec3>(*coarse_meca_mesh_, vertex_position, v) + time_step * s;
			}
		};
		solve_coarse();

		/////////////////////////////////////////////////
		////	         current_solve				/////
		/////////////////////////////////////////////////

		auto solve_current = [&]() {
			for (const Vertex& v : sc_fine_->vertices_cache)
			{
				value<Vec3>(*fine_meca_mesh_, this->forces_current_.get(), v) =
					value<Vec3>(*mecanical_mesh_, this->forces_ext_.get(), v);
			}
			sc_->solve_constraint(*mecanical_mesh_, vertex_position, this->forces_current_.get(), time_step);

			for (const Vertex& v : sc_->vertices_cache)
			{
				if (this->fixed_vertex && value<bool>(*mecanical_mesh_, this->fixed_vertex.get(), v))
				{
					value<Vec3>(*mecanical_mesh_, pos_current_, v) = value<Vec3>(*mecanical_mesh_, vertex_position, v);
					continue;
				}
				// compute speed
				Vec3 s = 0.995 * value<Vec3>(*mecanical_mesh_, this->speed_.get(), v) +
						 time_step * value<Vec3>(*mecanical_mesh_, this->forces_current_.get(), v) /
							 value<double>(*mecanical_mesh_, masse, v);
				s += time_step * gravity_;
				value<Vec3>(*mecanical_mesh_, pos_current_, v) =
					value<Vec3>(*mecanical_mesh_, vertex_position, v) + time_step * s;
			}
		};

		solve_current();

		/////////////////////////////////////////////////
		////	         fine_solve				/////
		/////////////////////////////////////////////////

		auto solve_fine = [&]() {
			sc_fine_->solve_constraint(*fine_meca_mesh_, vertex_position, this->forces_ext_.get(), time_step);

			for (const Vertex& v : sc_fine_->vertices_cache)
			{
				if (this->fixed_vertex && value<bool>(*fine_meca_mesh_, this->fixed_vertex.get(), v))
					continue;
				// compute speed
				value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) =
					0.995 * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) +
					time_step * value<Vec3>(*fine_meca_mesh_, this->forces_ext_.get(), v) /
						value<double>(*fine_meca_mesh_, sc_fine_->masse_, v);
				value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) += time_step * gravity_;
				value<Vec3>(*fine_meca_mesh_, vertex_position, v) =
					value<Vec3>(*fine_meca_mesh_, vertex_position, v) +
					time_step * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v);
			}
		};

		solve_fine();

		if (pc_)
		{
			pc_->propagate(*fine_meca_mesh_, m_geom, vertex_position, this->forces_ext_.get(), sc_fine_->masse_.get(),
						   relative_pos_.get(), parents_.get(), time_step);
		}

		foreach_cell(mecanical_mesh_->m_, [&](Vertex v) -> bool {
			value<Vec3>(m_geom, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			return true;
		});

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "\033[1;31m time step : \033[0m" << duration << std::endl;
		start = std::clock();
		compute_error(vertex_position, masse, time_step);
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "error : " << duration << std::endl;
		start = std::clock();
		modif_topo = update_topo(vertex_position) || modif_topo;
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "modif topo : " << duration << std::endl;
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
