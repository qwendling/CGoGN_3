#ifndef CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
#define CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/algos/Simulation_constraint.h>
#include <cgogn/simulation/algos/Simulation_solver.h>
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_constraint.h>
#include <forward_list>

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
	using Volume = typename mesh_traits<MR_MAP>::Volume;

	struct tree_volume
	{
		tree_volume* fils;
		tree_volume* pere;
		tree_volume* frere;
		Dart volume_dart;
		tree_volume() : fils(nullptr), pere(nullptr), frere(nullptr)
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
				if (!fn(it))
					return;
				it = fils->frere;
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
	std::forward_list<tree_volume*> list_volume_current_;
	std::forward_list<tree_volume*> list_volume_fine_;
	std::forward_list<tree_volume*> list_volume_coarse_;

	std::shared_ptr<Attribute<Vec3>> pos_tmp_;
	std::shared_ptr<Attribute<Vec3>> forces_tmp_;

	std::shared_ptr<Attribute<double>> volume_fine_;
	std::shared_ptr<Attribute<double>> volume_current_;
	std::shared_ptr<Attribute<double>> volume_coarse_;

	Simulation_solver_multiresolution()
		: Simulation_solver<MR_MAP>(), pc_(nullptr), parents_(nullptr), relative_pos_(nullptr)
	{
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

		std::vector<tree_volume*> tmp_watcher;

		list_volume_coarse_.clear();
		list_volume_fine_.clear();
		list_volume_current_.clear();
		int cmp_cur = 0, cmp_fine = 0, cmp_coarse = 0;
		CPH3 tmp(m);
		tmp.current_level_ = 0;
		std::function<void(tree_volume*, CPH3&, uint32)> progress_tree;
		progress_tree = [&](tree_volume* p, CPH3& cph, uint32 l) {
			uint32 cph_level = cph.volume_level(p->volume_dart);
			if (l == cph_level)
			{
				list_volume_current_.push_front(p);
				tmp_watcher.push_back(p);
				cmp_cur++;
				fine_meca_mesh_->activate_volume_subdivision(Volume(p->volume_dart));
				coarse_meca_mesh_->disable_volume_subdivision(Volume(p->volume_dart), true);
			}
			if (l == cph_level + 1)
			{
				list_volume_coarse_.push_front(p);
				cmp_coarse++;
			}
			if (l == cph_level - 1)
			{
				list_volume_fine_.push_front(p);
				cmp_fine++;
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

				progress_tree(t, cph, l);
			}
			cph.current_level_--;
		};
		foreach_cell(tmp, [&](CPH3::Volume v) -> bool {
			if (tmp.volume_is_subdivided(v.dart))
			{

				CPH3 tmp2(tmp);
				tree_volume* t = new tree_volume();
				t->volume_dart = v.dart;
				uint32 l = mecanical_mesh_->volume_level(v.dart);
				if (l == 0)
				{
					list_volume_coarse_.push_front(t);
					cmp_coarse++;
				}
				progress_tree(t, tmp2, l);
			}
			return true;
		});
		fine_meca_mesh_->current_level_ = mecanical_mesh_->current_level_ + 1;
		reset_forces(*fine_meca_mesh_);

		sc_ = sc;
		sc_->init_solver(*mecanical_mesh_, pos);
		sc_fine_ = std::shared_ptr<Simulation_constraint<MR_MAP>>(sc->get_new_ptr());
		sc_fine_->init_solver(*fine_meca_mesh_, pos);
		sc_coarse_ = std::shared_ptr<Simulation_constraint<MR_MAP>>(sc->get_new_ptr());
		sc_coarse_->init_solver(*coarse_meca_mesh_, pos);

		pos_tmp_ = get_attribute<Vec3, Vertex>(m, "Solver_multiresolution_pos_tmp");
		if (pos_tmp_ == nullptr)
			pos_tmp_ = add_attribute<Vec3, Vertex>(m, "Solver_multiresolution_pos_tmp");
		forces_tmp_ = get_attribute<Vec3, Vertex>(m, "Solver_multiresolution_forces_tmp");
		if (forces_tmp_ == nullptr)
			forces_tmp_ = add_attribute<Vec3, Vertex>(m, "Solver_multiresolution_forces_tmp");

		volume_fine_ = get_attribute<double, Volume>(m, "Solver_multiresolution_volume_fine");
		if (volume_fine_ == nullptr)
			volume_fine_ = add_attribute<double, Volume>(m, "Solver_multiresolution_volume_fine");
		volume_coarse_ = get_attribute<double, Volume>(m, "Solver_multiresolution_volume_coarse");
		if (volume_coarse_ == nullptr)
			volume_coarse_ = add_attribute<double, Volume>(m, "Solver_multiresolution_volume_coarse");
		volume_current_ = get_attribute<double, Volume>(m, "Solver_multiresolution_volume_current");
		if (volume_current_ == nullptr)
			volume_current_ = add_attribute<double, Volume>(m, "Solver_multiresolution_volume_current");
	}

	void reset_forces(MR_MAP& m)
	{
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			value<Vec3>(m, this->speed_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
	}

	void compute_time_step(MR_MAP& m_geom, Attribute<Vec3>* vertex_position, Attribute<double>* masse, double time_step)
	{
		if (!this->constraint_)
			return;

		/////////////////////////////////////////////////
		////	         coarse_solve				/////
		/////////////////////////////////////////////////
		parallel_foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool {
			value<Vec3>(*mecanical_mesh_, this->forces_tmp_.get(), v) =
				value<Vec3>(*coarse_meca_mesh_, this->forces_ext_.get(), v);
			return true;
		});
		sc_coarse_->solve_constraint(*coarse_meca_mesh_, vertex_position, this->forces_tmp_.get(), time_step);
		parallel_foreach_cell(*coarse_meca_mesh_, [&](Vertex v) -> bool {
			// compute speed
			auto s = 0.995 * value<Vec3>(*coarse_meca_mesh_, this->speed_.get(), v) +
					 time_step * value<Vec3>(*coarse_meca_mesh_, this->forces_tmp_.get(), v) /
						 value<double>(*coarse_meca_mesh_, sc_coarse_->masse_, v);
			value<Vec3>(*coarse_meca_mesh_, pos_tmp_, v) =
				value<Vec3>(*coarse_meca_mesh_, vertex_position, v) + time_step * s;
			return true;
		});

		if (!pc_)
			return;
		pc_->propagate(
			*coarse_meca_mesh_, *mecanical_mesh_, pos_tmp_.get(), this->forces_tmp_.get(), relative_pos_.get(),
			parents_.get(),
			[&](Vertex v) {
				auto s = 0.995 * value<Vec3>(*mecanical_mesh_, this->speed_.get(), v) +
						 time_step * value<Vec3>(*mecanical_mesh_, this->forces_tmp_.get(), v) /
							 value<double>(*mecanical_mesh_, sc_coarse_->masse_, v);
				value<Vec3>(*mecanical_mesh_, pos_tmp_.get(), v) =
					value<Vec3>(*mecanical_mesh_, vertex_position, v) + time_step * s;
			},
			time_step);

		parallel_foreach_cell(*mecanical_mesh_, [&](Volume v) -> bool {
			value<double>(*mecanical_mesh_, this->volume_current_.get(), v) =
				geometry::volume(*mecanical_mesh_, v, pos_tmp_.get());
			return true;
		});

		/////////////////////////////////////////////////
		////	         current_solve				/////
		/////////////////////////////////////////////////

		parallel_foreach_cell(*fine_meca_mesh_, [&](Vertex v) -> bool {
			value<Vec3>(*fine_meca_mesh_, this->forces_tmp_.get(), v) =
				value<Vec3>(*mecanical_mesh_, this->forces_ext_.get(), v);
			return true;
		});
		sc_->solve_constraint(*mecanical_mesh_, vertex_position, this->forces_tmp_.get(), time_step);
		parallel_foreach_cell(*mecanical_mesh_, [&](Vertex v) -> bool {
			// compute speed
			auto s = 0.995 * value<Vec3>(*mecanical_mesh_, this->speed_.get(), v) +
					 time_step * value<Vec3>(*mecanical_mesh_, this->forces_tmp_.get(), v) /
						 value<double>(*mecanical_mesh_, masse, v);
			value<Vec3>(*mecanical_mesh_, pos_tmp_, v) =
				value<Vec3>(*mecanical_mesh_, vertex_position, v) + time_step * s;
			return true;
		});

		if (!pc_)
			return;
		pc_->propagate(
			*mecanical_mesh_, *fine_meca_mesh_, pos_tmp_.get(), this->forces_tmp_.get(), relative_pos_.get(),
			parents_.get(),
			[&](Vertex v) {
				auto s = 0.995 * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) +
						 time_step * value<Vec3>(*fine_meca_mesh_, this->forces_tmp_.get(), v) /
							 value<double>(*fine_meca_mesh_, masse, v);
				value<Vec3>(*fine_meca_mesh_, pos_tmp_.get(), v) =
					value<Vec3>(*fine_meca_mesh_, vertex_position, v) + time_step * s;
			},
			time_step);

		parallel_foreach_cell(*mecanical_mesh_, [&](Volume v) -> bool {
			value<double>(*mecanical_mesh_, this->volume_current_.get(), v) =
				fabs(value<double>(*mecanical_mesh_, this->volume_current_.get(), v) -
					 geometry::volume(*mecanical_mesh_, v, pos_tmp_.get()));
			return true;
		});

		parallel_foreach_cell(*fine_meca_mesh_, [&](Volume v) -> bool {
			value<double>(*fine_meca_mesh_, this->volume_fine_.get(), v) =
				geometry::volume(*fine_meca_mesh_, v, pos_tmp_.get());
			return true;
		});

		/////////////////////////////////////////////////
		////	         fine_solve				/////
		/////////////////////////////////////////////////
		sc_fine_->solve_constraint(*fine_meca_mesh_, vertex_position, this->forces_ext_.get(), time_step);
		parallel_foreach_cell(*fine_meca_mesh_, [&](Vertex v) -> bool {
			std::cout << value<Vec3>(*fine_meca_mesh_, this->forces_ext_.get(), v) << std::endl;
			// compute speed
			value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) =
				0.995 * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v) +
				time_step * value<Vec3>(*fine_meca_mesh_, this->forces_ext_.get(), v) /
					value<double>(*fine_meca_mesh_, sc_fine_->masse_, v);
			value<Vec3>(*fine_meca_mesh_, vertex_position, v) =
				value<Vec3>(*fine_meca_mesh_, vertex_position, v) +
				time_step * value<Vec3>(*fine_meca_mesh_, this->speed_.get(), v);
			return true;
		});

		if (!pc_)
			return;
		pc_->propagate(
			*fine_meca_mesh_, m_geom, vertex_position, this->forces_tmp_.get(), relative_pos_.get(), parents_.get(),
			[&](Vertex v) {
				value<Vec3>(m_geom, this->speed_.get(), v) = 0.995 * value<Vec3>(m_geom, this->speed_.get(), v) +
															 time_step *
																 value<Vec3>(m_geom, this->forces_ext_.get(), v) /
																 value<double>(m_geom, sc_fine_->masse_, v);
				value<Vec3>(m_geom, vertex_position, v) =
					value<Vec3>(m_geom, vertex_position, v) + time_step * value<Vec3>(m_geom, this->speed_.get(), v);
			},
			time_step);

		parallel_foreach_cell(*fine_meca_mesh_, [&](Volume v) -> bool {
			value<double>(*fine_meca_mesh_, this->volume_fine_.get(), v) =
				fabs(value<double>(*fine_meca_mesh_, this->volume_fine_.get(), v) -
					 geometry::volume(*fine_meca_mesh_, v, vertex_position));
			return true;
		});

		parallel_foreach_cell(m_geom, [&](Vertex v) -> bool {
			value<Vec3>(m_geom, this->forces_ext_.get(), v) = Vec3(0, 0, 0);
			return true;
		});
	}
};
} // namespace simulation
} // namespace cgogn

#endif // CGOGN_SIMULATION_SIMULATION_SOLVER_MULTIRESOLUTION_H
