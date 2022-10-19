/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>

#include <GLFW/glfw3.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/core/types/cmap/phi.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_emr_modeling/volume_emr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

#include <cgogn/simulation/algos/Simulation_solver_multiresolution.h>
#include <cgogn/simulation/algos/lattice_shape_matching/lattice_shape_matching.h>
#include <cgogn/simulation/algos/multiresolution_propagation/propagation_plastique.h>

#include <cgogn/core/utils/thread.h>

using MRMesh = cgogn::EMR_Map3_Adaptative;
using Mesh = MRMesh::BASE;
using EMR_Map3 = cgogn::EMR_Map3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<MRMesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<MRMesh>::Vertex;
using Edge = typename cgogn::mesh_traits<MRMesh>::Edge;
using Face = typename cgogn::mesh_traits<MRMesh>::Face;
using Volume = typename cgogn::mesh_traits<MRMesh>::Volume;
using Dart = cgogn::Dart;

using Vec3 = cgogn::geometry::Vec3;

void test_foreach_dart(MRMesh* mrm, Volume v)
{
	foreach_dart_of_orbit(*mrm, v, [&](Dart) -> bool { return true; });
};

void foreachcell_benchmark(MRMesh* mrm)
{
	std::clock_t start;
	std::clock_t start2;
	double duration;
	double duration2;

	start = std::clock();
	duration2 = 0;
	cgogn::CellMarker<MRMesh, Vertex> cm(*mrm);
	for (Dart d = mrm->begin(), end = mrm->end(); d != end; d = mrm->next(d))
	{
		const Vertex c(d);
		start2 = std::clock();
		if (!is_boundary(*mrm, d) && !cm.is_marked(c))
		{
			cm.mark(c);
		}
		duration2 += (std::clock() - start2) / (double)CLOCKS_PER_SEC;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours vertex : " << duration << std::endl;
	std::cout << "temps it parcours vertex : " << duration2 << std::endl;

	start = std::clock();
	duration2 = 0;
	cgogn::CellMarker<MRMesh, Edge> cm2(*mrm);
	for (Dart d = mrm->begin(), end = mrm->end(); d != end; d = mrm->next(d))
	{
		const Edge c(d);
		start2 = std::clock();
		if (!is_boundary(*mrm, d) && !cm2.is_marked(c))
		{
			cm2.mark(c);
		}
		duration2 += (std::clock() - start2) / (double)CLOCKS_PER_SEC;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours edge : " << duration << std::endl;
	std::cout << "temps it parcours edge : " << duration2 << std::endl;

	start = std::clock();
	duration2 = 0;
	cgogn::CellMarker<MRMesh, Face> cm3(*mrm);
	for (Dart d = mrm->begin(), end = mrm->end(); d != end; d = mrm->next(d))
	{
		const Face c(d);
		start2 = std::clock();
		if (!is_boundary(*mrm, d) && !cm3.is_marked(c))
		{
			cm3.mark(c);
		}
		duration2 += (std::clock() - start2) / (double)CLOCKS_PER_SEC;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours face : " << duration << std::endl;
	std::cout << "temps it parcours face : " << duration2 << std::endl;

	start = std::clock();
	duration2 = 0;
	cgogn::CellMarker<MRMesh, Volume> cm4(*mrm);
	for (Dart d = mrm->begin(), end = mrm->end(); d != end; d = mrm->next(d))
	{
		const Volume c(d);
		start2 = std::clock();
		if (!is_boundary(*mrm, d) && !cm4.is_marked(c))
		{
			cm4.mark(c);
		}
		duration2 += (std::clock() - start2) / (double)CLOCKS_PER_SEC;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours volume : " << duration << std::endl;
	std::cout << "temps it parcours volume : " << duration2 << std::endl;
}

void test(MRMesh* mrm, Attribute<Vec3>* attr)
{
	std::clock_t start;
	double duration;

	// foreachcell_benchmark( mrm);

	start = std::clock();

	cgogn::foreach_cell(*mrm, [&](Vertex) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours vertex : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Edge) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours edge : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Face) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours face : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours volume : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours volume 2 : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
		cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
			cgogn::foreach_incident_vertex(*mrm, w, [&](Vertex) -> bool { return true; });
			return true;
		});
		return true;
	});

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorA : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		cgogn::foreach_incident_vertex(*mrm, v, [&](Vertex w) -> bool {
			cgogn::foreach_incident_volume(*mrm, w, [&](Volume) -> bool { return true; });
			return true;
		});
		return true;
	});

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorB : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		cgogn::geometry::centroid<Vec3>(*mrm, v, attr);
		return true;
	});

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps Centroid : " << duration << std::endl;
};

void test2(cgogn::CMap3* mrm, Attribute<Vec3>* attr)
{
	std::clock_t start;
	double duration;

	start = std::clock();

	cgogn::foreach_cell(*mrm, [&](Vertex) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours vertex : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Edge) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours edge : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Face) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours face : " << duration << std::endl;
	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume) -> bool { return true; });
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps parcours volume : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
		cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
			cgogn::foreach_incident_vertex(*mrm, w, [&](Vertex) -> bool { return true; });
			return true;
		});
		return true;
	});

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorA : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		cgogn::foreach_incident_vertex(*mrm, v, [&](Vertex w) -> bool {
			cgogn::foreach_incident_volume(*mrm, w, [&](Volume) -> bool { return true; });
			return true;
		});
		return true;
	});

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorB : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		Vec3 c = cgogn::geometry::centroid<Vec3>(*mrm, v, attr);
		return true;
	});

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps Centroid : " << duration << std::endl;
};

void activate_all_volume(MRMesh* m)
{
	std::vector<Volume> vec_vol;
	cgogn::foreach_cell(*m, [&](Volume v) -> bool {
		vec_vol.push_back(Volume(m->volume_youngest_dart(v.dart)));
		return true;
	});
	for (auto v : vec_vol)
	{
		m->activate_volume_subdivision_fast(v);
	}
}

int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " volume_mesh_file" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	std::string ext = cgogn::extension(filename);
	Mesh* m = new Mesh();
	cgogn::CMap3* m2 = new cgogn::CMap3();
	bool imported = false;

	if (ext.compare("tet") == 0)
	{
		imported = cgogn::io::import_TET(*m, filename);
		imported = cgogn::io::import_TET(*m2, filename);
	}
	else if (ext.compare("mesh") == 0 || ext.compare("meshb") == 0)
	{
		imported = cgogn::io::import_MESHB(*m, filename);
		imported = cgogn::io::import_MESHB(*m2, filename);
	}
	else
		imported = false;

	if (!imported)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	MRMesh* mrm = new MRMesh(*m);

	cgogn::index_cells<Mesh::Face>(*mrm);
	cgogn::index_cells<Mesh::Volume>(*mrm);
	cgogn::index_cells<Mesh::Edge>(*mrm);

	std::shared_ptr<Attribute<Vec3>> positionCmap = cgogn::get_attribute<Vec3, Vertex>(*m2, "position");
	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");
	std::shared_ptr<Attribute<Vec3>> force = cgogn::add_attribute<Vec3, Vertex>(*mrm, "force");
	std::shared_ptr<Attribute<double>> volume = cgogn::add_attribute<double, Volume>(*mrm, "volume");
	std::shared_ptr<Attribute<std::array<Vertex, 4>>> parent =
		cgogn::add_attribute<std::array<Vertex, 4>, Vertex>(*m, "parents");
	std::shared_ptr<Attribute<Vec3>> relative_pos = cgogn::add_attribute<Vec3, Vertex>(*m, "relative_position");

	std::clock_t start;
	double duration;

	m->add_resolution();
	mrm->change_resolution_level(1);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, parent.get(), relative_pos.get());

	m->add_resolution();
	mrm->change_resolution_level(2);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, parent.get(), relative_pos.get());

	m->add_resolution();
	mrm->change_resolution_level(3);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, parent.get(), relative_pos.get());

	mrm->change_resolution_level(0);

	std::cout << "Resolution 0 : " << std::endl;
	test(mrm, position.get());
	test2(m2, position.get());

	start = std::clock();
	activate_all_volume(mrm);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps activate mrmap : " << duration << std::endl;
	start = std::clock();
	cgogn::modeling::butterflySubdivisionVolumeRegular(*m2, 0.0f, {position.get()});
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps subdivise cmap : " << duration << std::endl;

	std::cout << "Resolution 1 : " << std::endl;
	test(mrm, position.get());
	test2(m2, position.get());

	start = std::clock();
	activate_all_volume(mrm);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps activate mrmap : " << duration << std::endl;
	start = std::clock();
	cgogn::modeling::butterflySubdivisionVolumeRegular(*m2, 0.0f, {position.get()});
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps subdivise cmap : " << duration << std::endl;

	std::cout << "Resolution 2 : " << std::endl;
	test(mrm, position.get());
	test2(m2, position.get());

	start = std::clock();
	activate_all_volume(mrm);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps activate mrmap : " << duration << std::endl;
	start = std::clock();
	cgogn::modeling::butterflySubdivisionVolumeRegular(*m2, 0.0f, {position.get()});
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps subdivise cmap : " << duration << std::endl;

	std::cout << "Resolution 3 : " << std::endl;
	test(mrm, position.get());
	test2(m2, position.get());

	mrm->change_resolution_level(1);
	// std::srand(std::time(nullptr));
	std::srand(2124512438);

#define CELL_RANDOM Volume
#if 0

	auto nb_cells = [&]() -> int {
		int result = 0;
		cgogn::foreach_cell(*mrm, [&](CELL_RANDOM) -> bool {
			result++;
			return true;
		});
		return result;
	};

#define N 10

	for (int i = 0; i < N; i++)
	{
		std::vector<CELL_RANDOM> vec_volume;

		cgogn::foreach_cell(*mrm, [&](CELL_RANDOM v) -> bool {
			if ((rand() / (double)RAND_MAX) * 100 < 10)
			{
				vec_volume.push_back(v);
			}
			return true;
		});
		std::clock_t start;
		double duration;
		double diff_volume = nb_cells();

		start = std::clock();
		for (auto v : vec_volume)
		{
			mrm->activate_volume_subdivision(v);
			assert(mrm->check_integrity());
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

		std::cout << "temps activate " << nb_cells() - diff_volume << " volume : " << duration << std::endl;
	}
#endif

	/*cgogn::simulation::Simulation_solver_multiresolution<MRMesh> simu_solver;
	cgogn::simulation::Propagation_Plastique<MRMesh> ps_;

	cgogn::simulation::lattice_shape_matching_constraint_solver<MRMesh> sm_solver_(0.9);
	sm_solver_.init_solver(*mrm, position.get());
	sm_solver_.solve_constraint(*mrm, position.get(), force.get(), 0.005f);

	cgogn::simulation::lattice_shape_matching_constraint_solver<cgogn::CMap3> sm_solver2_(0.9);
	sm_solver2_.init_solver(*mrm->get_map(), position.get());
	sm_solver2_.solve_constraint(*mrm->get_map(), position.get(), force.get(), 0.005f);

	simu_solver.init_solver(*mrm, &sm_solver_, position.get(), &ps_);
	simu_solver.parents_ = parent;
	simu_solver.relative_pos_ = relative_pos;

	bool tmp;
	for (int i = 0; i < 100; i++)
	{
		simu_solver.compute_time_step(*mrm, position.get(), sm_solver_.masse_.get(), 0.005f, tmp);
	}*/

	/*cgogn::launch_thread([&]() {
		bool tmp;
		for (int i = 0; i < 100; i++)
		{
			simu_solver.compute_time_step(*mrm, position.get(), sm_solver_.masse_.get(), 0.005f, tmp);
		}
	});
	cgogn::launch_thread([&]() {
		bool tmp;
		for (int i = 0; i < 100; i++)
		{
			simu_solver.compute_time_step(*mrm, position.get(), sm_solver_.masse_.get(), 0.005f, tmp);
		}
	});
	cgogn::launch_thread([&]() {
		bool tmp;
		for (int i = 0; i < 100; i++)
		{
			simu_solver.compute_time_step(*mrm, position.get(), sm_solver_.masse_.get(), 0.005f, tmp);
		}
	});*/

	/*std::clock_t start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		cgogn::value<double>(*mrm, volume.get(), v) = cgogn::geometry::volume(*mrm, v, position.get());
		return true;
	});
	double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "temps volume MR : " << duration << std::endl;

	start = std::clock();
	cgogn::foreach_cell(*map, [&](Volume v) -> bool {
		cgogn::value<double>(*map, volume.get(), v) = cgogn::geometry::volume(*map, v, position.get());
		return true;
	});
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "temps volume CMAP : " << duration << std::endl;

	mrm->change_resolution_level(0);
	std::vector<Volume> vect_vol;
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		vect_vol.push_back(v);
		return true;
	});

	for (Volume v : vect_vol)
	{
		mrm->activate_volume_subdivision(v);
	}

	test(mrm);

	sm_solver_.init_solver(*mrm, position.get());
	sm_solver_.solve_constraint(*mrm, position.get(), force.get(), 0.005f);

	start = std::clock();
	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		cgogn::value<double>(*mrm, volume.get(), v) = cgogn::geometry::volume(*mrm, v, position.get());
		return true;
	});
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "temps volume MR view : " << duration << std::endl;*/
	/*bool tmp;
	for (int i = 0; i < 1000; i++)
	{
		simu_solver.compute_time_step(*mrm, position.get(), sm_solver_.masse_.get(), 0.005f, tmp);
	}*/
	delete m2;
	delete mrm;
	delete m;

	return 0;
}
