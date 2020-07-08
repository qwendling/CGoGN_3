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

#define QUOTA_VOLUME 1000
#define PROPORTION_QUOTA 0.1f
#define MODIF_MAX PROPORTION_QUOTA* QUOTA_VOLUME

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

	struct tree_volume
	{
		tree_volume* fils;
		tree_volume* pere;
		tree_volume* frere;
		Dart volume_dart;
		bool is_current;
		int clock;
		tree_volume() : fils(nullptr), pere(nullptr), frere(nullptr), is_current(false), clock(0)
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
	};

public:
	Propagation_Constraint<MR_MAP>* pc_;
	Simulation_constraint<MR_MAP>* sc_;
	std::shared_ptr<Simulation_constraint<MR_MAP>> sc_fine_;
	std::shared_ptr<Simulation_constraint<MR_MAP>> sc_coarse_;
	std::shared_ptr<Attribute<std::array<Vertex, 3>>> parents_;
	std::shared_ptr<Attribute<Vec3>> relative_pos_;
	MR_MAP* mecanical_mesh_;
	MR_MAP* fine_meca_mesh_;
	MR_MAP* coarse_meca_mesh_;
	std::forward_list<tree_volume*> list_volume_coarse_;
	std::forward_list<tree_volume*> list_volume_current_;

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
		: Simulation_solver<MR_MAP>(), pc_(nullptr), parents_(nullptr), relative_pos_(nullptr), gravity_(0, 0, -9.81f),
		  cache_current_vol_(nullptr), clock(1)
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

	void init_solver(MR_MAP& m, Simulation_constraint<MR_MAP>* sc, Attribute<Vec3>* pos,
					 Propagation_Constraint<MR_MAP>* pc = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& speed = nullptr,
					 const std::shared_ptr<Attribute<Vec3>>& forces = nullptr)
	{
		Inherit::init_solver(m, sc, speed, forces);
		pc_ = pc;
		mecanical_mesh_ = &m;
		fine_meca_mesh_ = m.get_child();
		coarse_meca_mesh_ = m.get_copy();

		if (cache_current_vol_ != nullptr)
			delete cache_current_vol_;
		cache_current_vol_ = new CellCache<MR_MAP>(*mecanical_mesh_);

		std::vector<tree_volume*> tmp_watcher;

		list_volume_coarse_.clear();
		list_volume_current_.clear();
		int cmp_cur = 0;
		CPH3 tmp(m);
		tmp.current_level_ = m.current_level_;
		std::function<void(tree_volume*, CPH3&)> progress_tree;
		progress_tree = [&](tree_volume* p, CPH3& cph) {
			uint32 cph_level = cph.volume_level(p->volume_dart);
			uint32 l = 0;
			if (mecanical_mesh_->dart_is_visible(p->volume_dart))
				l = mecanical_mesh_->volume_level(p->volume_dart);
			if (l == cph_level)
			{
				list_volume_current_.push_front(p);
				tmp_watcher.push_back(p);
				cmp_cur++;
				p->is_current = true;
			}
			if (!cph.volume_is_subdivided(p->volume_dart))
			{
				return;
			}

			std::vector<Volume> sub_volume;
			foreach_incident_vertex(cph, Volume(p->volume_dart), [&](Vertex v) -> bool {
				sub_volume.push_back(Volume(v.dart));
				return true;
			});
			cph.current_level_++;
			for (auto v : sub_volume)
			{
				tree_volume* t = new tree_volume();
				t->volume_dart = v.dart;
				t->frere = p->fils;
				p->fils = t;
				t->pere = p;

				progress_tree(t, cph);
			}
			cph.current_level_--;
		};
		foreach_cell(tmp, [&](CPH3::Volume v) -> bool {
			CPH3 tmp2(tmp);
			tree_volume* t = new tree_volume();
			t->volume_dart = v.dart;
			progress_tree(t, tmp2);
			return true;
		});
		// init_cell_cache(*cache_current_vol_, list_volume_current_);

		for (tree_volume* t : list_volume_current_)
		{
			if (t->pere && t->pere->clock != clock)
			{
				t->pere->clock = clock;
				bool can_be_coarse = true;
				t->pere->for_each_child([&](tree_volume* c) {
					if (!c->is_current)
						can_be_coarse = false;
					return can_be_coarse;
				});
				if (can_be_coarse)
				{
					coarse_meca_mesh_->disable_volume_subdivision(Volume(t->pere->volume_dart), true);
					list_volume_coarse_.push_front(t->pere);
				}
			}
			fine_meca_mesh_->activate_volume_subdivision(Volume(t->volume_dart));
		}

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

	void update_topo()
	{
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

		tree_volume* min_coarse = list_volume_coarse_.front();
		tree_volume* max_fine = list_volume_current_.front();
		double v1 =
			value<double>(*mecanical_mesh_, this->diff_volume_current_fine_.get(), Volume(max_fine->volume_dart));
		double v2 =
			value<double>(*coarse_meca_mesh_, this->diff_volume_coarse_current_.get(), Volume(min_coarse->volume_dart));
		int nb_modif = 0;
		if (min_coarse->is_current)
			return;
		std::cout << v1 << std::endl;
		std::cout << v2 << std::endl;
		while (v2 - v1 < -0.00001 && nb_modif < MODIF_MAX)
		{
			list_volume_current_.pop_front();
			list_volume_coarse_.pop_front();
			// Activation
			if (max_fine->pere)
			{
				coarse_meca_mesh_->activate_volume_subdivision(Volume(max_fine->pere->volume_dart));
				list_volume_coarse_.remove(max_fine->pere);
				max_fine->pere->for_each_child([&](tree_volume* c) -> bool {
					list_volume_current_.remove(c);
					return true;
				});
			}
			list_new_coarse_.push_front(max_fine);
			max_fine->is_current = false;
			mecanical_mesh_->activate_volume_subdivision(Volume(max_fine->volume_dart));
			max_fine->for_each_child([&](tree_volume* c) -> bool {
				fine_meca_mesh_->activate_volume_subdivision(Volume(c->volume_dart));
				list_new_current_.push_front(c);
				return true;
			});

			// Disable
			min_coarse->is_current = true;
			if (min_coarse->pere)
			{

				bool can_be_coarse = true;
				min_coarse->pere->for_each_child([&](tree_volume* c) {
					if (!c->is_current)
						can_be_coarse = false;
					return can_be_coarse;
				});
				if (can_be_coarse)
				{
					coarse_meca_mesh_->disable_volume_subdivision(Volume(min_coarse->volume_dart), true);
					min_coarse->pere->for_each_child([&](tree_volume* t) -> bool {
						list_new_current_.push_front(t);
						return true;
					});
					list_new_coarse_.push_front(min_coarse->pere);
				}
			}
			else
			{
				list_new_current_.push_front(min_coarse);
			}

			min_coarse->for_each_child([&](tree_volume* t) -> bool {
				list_volume_current_.remove(t);
				return true;
			});

			mecanical_mesh_->disable_volume_subdivision(Volume(min_coarse->volume_dart), true);
			min_coarse->for_each_child([&](tree_volume* c) -> bool {
				if (c->fils)
					fine_meca_mesh_->disable_volume_subdivision(Volume(c->volume_dart), true);
				return true;
			});

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
		if (!list_new_current_.empty())
		{
			std::cout << "hello" << std::endl;
			sc_coarse_->update_topo(*coarse_meca_mesh_, {});
			sc_->update_topo(*mecanical_mesh_, {});
			sc_fine_->update_topo(*fine_meca_mesh_, {});
		}
		else
		{
			std::cout << "bye" << std::endl;
		}
	}

	void compute_time_step(MR_MAP& m_geom, Attribute<Vec3>* vertex_position, Attribute<double>* masse, double time_step)
	{
		if (!this->constraint_)
			return;

		ThreadPool* pool = thread_pool();
		std::clock_t start;
		double duration;
		start = std::clock();
		/////////////////////////////////////////////////
		////	         coarse_solve				/////
		/////////////////////////////////////////////////

		auto solve_coarse = [&]() {
			parallel_foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool {
				value<Vec3>(*mecanical_mesh_, this->forces_coarse_.get(), v) =
					value<Vec3>(*coarse_meca_mesh_, this->forces_ext_.get(), v);
				return true;
			});
			sc_coarse_->solve_constraint(*coarse_meca_mesh_, vertex_position, this->forces_coarse_.get(), time_step);
			parallel_foreach_cell(*coarse_meca_mesh_, [&](Vertex v) -> bool {
				// compute speed
				Vec3 s = 0.995 * value<Vec3>(*coarse_meca_mesh_, this->speed_.get(), v) +
						 time_step * value<Vec3>(*coarse_meca_mesh_, this->forces_coarse_.get(), v) /
							 value<double>(*coarse_meca_mesh_, sc_coarse_->masse_, v);
				s += time_step * gravity_;
				value<Vec3>(*coarse_meca_mesh_, pos_coarse_, v) =
					value<Vec3>(*coarse_meca_mesh_, vertex_position, v) + time_step * s;
				return true;
			});
			if (!pc_)
				return;
			pc_->propagate(*coarse_meca_mesh_, *mecanical_mesh_, pos_coarse_.get(), nullptr, this->forces_coarse_.get(),
						   masse, relative_pos_.get(), parents_.get(), time_step);
		};

		auto future_solve_coarse = pool->enqueue(solve_coarse);

		/////////////////////////////////////////////////
		////	         current_solve				/////
		/////////////////////////////////////////////////

		auto solve_current = [&]() {
			parallel_foreach_cell(*fine_meca_mesh_, [&](Vertex v) -> bool {
				value<Vec3>(*fine_meca_mesh_, this->forces_current_.get(), v) =
					value<Vec3>(*mecanical_mesh_, this->forces_ext_.get(), v);
				return true;
			});

			sc_->solve_constraint(*mecanical_mesh_, vertex_position, this->forces_current_.get(), time_step);
			parallel_foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool {
				// compute speed
				Vec3 s = 0.995 * value<Vec3>(*mecanical_mesh_, this->speed_.get(), v) +
						 time_step * value<Vec3>(*mecanical_mesh_, this->forces_current_.get(), v) /
							 value<double>(*mecanical_mesh_, masse, v);
				s += time_step * gravity_;
				value<Vec3>(*mecanical_mesh_, pos_current_, v) =
					value<Vec3>(*mecanical_mesh_, vertex_position, v) + time_step * s;
				return true;
			});

			if (!pc_)
				return;
			pc_->propagate(*mecanical_mesh_, *fine_meca_mesh_, pos_current_.get(), nullptr, this->forces_current_.get(),
						   sc_fine_->masse_.get(), relative_pos_.get(), parents_.get(), time_step);
		};

		auto future_solve_current = pool->enqueue(solve_current);

		/////////////////////////////////////////////////
		////	         fine_solve				/////
		/////////////////////////////////////////////////

		auto solve_fine = [&]() {
			sc_fine_->solve_constraint(*fine_meca_mesh_, vertex_position, this->forces_ext_.get(), time_step);
			parallel_foreach_cell(*fine_meca_mesh_, [&](Vertex v) -> bool {
				// compute speed
				value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) =
					0.995 * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) +
					time_step * value<Vec3>(*fine_meca_mesh_, this->forces_ext_.get(), v) /
						value<double>(*fine_meca_mesh_, sc_fine_->masse_, v);
				value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) += time_step * gravity_;
				value<Vec3>(*fine_meca_mesh_, vertex_position, v) =
					value<Vec3>(*fine_meca_mesh_, vertex_position, v) +
					time_step * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v);
				return true;
			});

			if (!pc_)
				return;
			pc_->propagate(*fine_meca_mesh_, m_geom, vertex_position, this->speed_.get(), this->forces_ext_.get(),
						   sc_fine_->masse_.get(), relative_pos_.get(), parents_.get(), time_step);
		};

		future_solve_coarse.wait();
		future_solve_current.wait();
		auto future_solve_fine = pool->enqueue(solve_fine);

		/////////////////////////////////////////////////
		////	         error_criteria				/////
		/////////////////////////////////////////////////

		auto volume_current = [&]() {
			future_solve_coarse.wait();
			auto fn = [&](tree_volume* t) -> bool {
				value<double>(*mecanical_mesh_, this->diff_volume_coarse_current_.get(), Volume(t->volume_dart)) =
					geometry::volume(*mecanical_mesh_, Volume(t->volume_dart), pos_coarse_.get());
				return true;
			};
			for (tree_volume* tp : list_volume_coarse_)
			{
				tp->for_each_child(fn);
			}
		};
		auto future_volume_current = pool->enqueue(volume_current);

		auto volume_fine = [&]() {
			future_solve_current.wait();
			auto fn = [&](tree_volume* t) -> bool {
				value<double>(*fine_meca_mesh_, this->diff_volume_current_fine_.get(), Volume(t->volume_dart)) =
					geometry::volume(*fine_meca_mesh_, Volume(t->volume_dart), pos_current_.get());
				return true;
			};
			for (tree_volume* tp : list_volume_current_)
			{
				tp->for_each_child(fn);
			}
		};
		auto future_volume_fine = pool->enqueue(volume_fine);

		auto error_coarse_current = [&]() {
			future_volume_current.wait();
			future_solve_current.wait();
			auto fn = [&](tree_volume* t) -> bool {
				Volume v = Volume(t->volume_dart);
				double vol = geometry::volume(*mecanical_mesh_, v, vertex_position);
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
				if (t->is_current)
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
		};

		auto future_error_cc = pool->enqueue(error_coarse_current);

		auto error_current_fine = [&]() {
			future_volume_fine.wait();
			future_solve_fine.wait();
			auto fn2 = [&](tree_volume* t) -> bool {
				Volume v = Volume(t->volume_dart);
				double vol = geometry::volume(*fine_meca_mesh_, v, vertex_position);
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
		};

		auto future_error_cf = pool->enqueue(error_current_fine);

		parallel_foreach_cell(mecanical_mesh_->m_, [&](Vertex v) -> bool {
			value<Vec3>(m_geom, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
		future_error_cc.wait();
		future_error_cf.wait();
		update_topo();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "time step : " << duration << std::endl;
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
