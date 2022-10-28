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
using MAP = cgogn::CPH3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<MRMesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<MRMesh>::Vertex;
using Edge = typename cgogn::mesh_traits<MRMesh>::Edge;
using Face = typename cgogn::mesh_traits<MRMesh>::Face;
using Volume = typename cgogn::mesh_traits<MRMesh>::Volume;
using Dart = cgogn::Dart;

using Vec3 = cgogn::geometry::Vec3;

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

	std::cout << "subdivision mr" << std::endl;
	m->add_resolution();
	mrm->change_resolution_level(1);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, parent.get(), relative_pos.get());

	m->add_resolution();
	mrm->change_resolution_level(2);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, parent.get(), relative_pos.get());

	m->add_resolution();
	mrm->change_resolution_level(3);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, parent.get(), relative_pos.get());

	mrm->change_resolution_level(2);

	std::cout << "subdivision mono" << std::endl;
	cgogn::modeling::butterflySubdivisionVolumeRegular(*m2, 0.0f, {position.get()});
	cgogn::modeling::butterflySubdivisionVolumeRegular(*m2, 0.0f, {position.get()});
	MAP* cph = new MAP(*m2);

	std::srand(2124512438);

	std::list<Volume> volume_to_subdivided;
	std::list<Volume> volume_to_simplified;

	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		volume_to_subdivided.push_back(Volume(mrm->volume_youngest_dart(v.dart)));
		return true;
	});

	std::cout << "debut benchmark" << std::endl;

	std::cout << "nb_volume;subdivide mr;subdivide mono;simplified mr;simplified mono" << std::endl;

	std::vector<Volume> choix_volume;
	int nb_cell = cgogn::nb_cells<Volume>(*mrm);

	for (int i = 0; i < 1000; i++)
	{
		std::cout << nb_cell << ";";
		for (auto it = volume_to_subdivided.begin(); it != volume_to_subdivided.end();)
		{
			if (rand() % 100 < 10)
			{
				Volume v = *it;
				choix_volume.push_back(v);
				volume_to_simplified.push_back(v);
				it = volume_to_subdivided.erase(it);
				nb_cell += 7;
			}
			else
			{
				++it;
			}
		}
		start = std::clock();
		for (Volume v : choix_volume)
		{
			mrm->activate_volume_subdivision_fast(v);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";
		// std::cout << "temps subdivided 10% volume mr : " << duration << std::endl;

		start = std::clock();
		for (Volume v : choix_volume)
		{
			auto fn = [](Vertex) {};
			cgogn::modeling::butterflySubdivisionVolume(*cph, 0.0f, {position.get()}, v, fn, fn, fn);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";
		// std::cout << "temps subdivided 10% volume mono : " << duration << std::endl;

		choix_volume.clear();

		for (auto it = volume_to_simplified.begin(); it != volume_to_simplified.end();)
		{
			if (rand() % 100 < 10)
			{
				Volume v = *it;
				choix_volume.push_back(v);
				volume_to_subdivided.push_back(v);
				it = volume_to_simplified.erase(it);
				nb_cell -= 7;
			}
			else
			{
				++it;
			}
		}
		start = std::clock();
		for (Volume v : choix_volume)
		{
			mrm->disable_volume_subdivision_fast(Volume(cgogn::phi1(*mrm, v.dart)), true);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";
		// std::cout << "temps simplified 10% volume mr : " << duration << std::endl;

		start = std::clock();
		for (Volume v : choix_volume)
		{
			cph->disable_volume_subdivision(v, true);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << std::endl;
		// std::cout << "temps simplified 10% volume mono : " << duration << std::endl;

		choix_volume.clear();
	}

	delete m2;
	delete mrm;
	delete m;

	return 0;
}
