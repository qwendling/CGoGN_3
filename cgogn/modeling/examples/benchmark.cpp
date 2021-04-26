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
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_emr_modeling/volume_emr_modeling.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/volume_selection/volume_selection.h>

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

void test(MRMesh* mrm)
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
	duration = 0;
	for (cgogn::Dart d = mrm->begin(), e = mrm->end(); d != e; d = mrm->next(d))
	{
		start = std::clock();

		if (mrm->edge_level(d) != 0)
		{
			foreach_dart_of_orbit(*mrm, Face(d), [&](Dart) -> bool { return true; });
			duration += (std::clock() - start) / (double)CLOCKS_PER_SEC;
		}
	}
	std::cout << "temps face foreach dart : " << duration << std::endl;
	start = std::clock();
	duration = 0;
	for (cgogn::Dart d = mrm->begin(), e = mrm->end(); d != e; d = mrm->next(d))
	{
		start = std::clock();

		if (mrm->edge_level(d) != 0)
		{
			test_foreach_dart(mrm, Volume(d));
			duration += (std::clock() - start) / (double)CLOCKS_PER_SEC;
		}
	}
	std::cout << "temps volume foreach dart : " << duration << std::endl;
	start = std::clock();
	for (cgogn::Dart d = mrm->begin(), e = mrm->end(); d != e; d = mrm->next(d))
	{
		cgogn::phi1(*mrm, d);
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps phi1 : " << duration << std::endl;
	start = std::clock();
	for (cgogn::Dart d = mrm->begin(), e = mrm->end(); d != e; d = mrm->next(d))
	{
		phi2(*mrm, d);
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps phi2 : " << duration << std::endl;
	start = std::clock();
	for (cgogn::Dart d = mrm->begin(), e = mrm->end(); d != e; d = mrm->next(d))
	{
		phi3(*mrm, d);
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps phi3 : " << duration << std::endl;
};

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
	bool imported = false;

	if (ext.compare("tet") == 0)
		imported = cgogn::io::import_TET(*m, filename);
	else if (ext.compare("mesh") == 0 || ext.compare("meshb") == 0)
		imported = cgogn::io::import_MESHB(*m, filename);
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

	std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*mrm, "position");

	m->add_resolution();
	mrm->change_resolution_level(1);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, nullptr, nullptr);
	m->add_resolution();
	mrm->change_resolution_level(2);
	cgogn::modeling::butterflyMultiresolution(*mrm, 0.34f, {position.get()}, nullptr, nullptr);

	mrm->change_resolution_level(0);
	// std::srand(std::time(nullptr));
	std::srand(2124512438);

#define CELL_RANDOM Volume

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
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		// std::cout << "temps activate " << nb_cells() - diff_volume << " volume : " << duration << std::endl;
	}

	test(mrm);

	return 0;
}
