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
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/ui_modules/volume_selection.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/ui_modules/volume_emr_modeling.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/simulation/ui_modules/animation_multiresolution.h>

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

	start = std::clock();
	for (int i = 0; i < 10; i++)
	{
		cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
			cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
				cgogn::foreach_incident_vertex(*mrm, w, [&](Vertex) -> bool { return true; });
				return true;
			});
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorA : " << duration / 10.0f << std::endl;

	start = std::clock();

	for (int i = 0; i < 10; i++)
	{
		cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
			cgogn::foreach_incident_vertex(*mrm, v, [&](Vertex w) -> bool {
				cgogn::foreach_incident_volume(*mrm, w, [&](Volume) -> bool { return true; });
				return true;
			});
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorB : " << duration / 10.0f << std::endl;

	start = std::clock();
	for (int i = 0; i < 10; i++)
	{
		cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
			cgogn::geometry::centroid<Vec3>(*mrm, v, attr);
			return true;
		});
		cgogn::foreach_cell(*mrm, [&](cgogn::CMap3::Vertex v) -> bool {
			cgogn::CellMarkerStore<MRMesh, cgogn::CMap3::Vertex> mv(*mrm);
			Vec3 cm(0, 0, 0);
			cgogn::foreach_incident_volume(*mrm, v, [&](cgogn::CMap3::Volume w) -> bool {
				cgogn::foreach_incident_vertex(*mrm, w, [&](cgogn::CMap3::Vertex v2) -> bool {
					if (!mv.is_marked(v2))
					{
						cm += cgogn::value<Vec3>(*mrm, attr, v2);
						mv.mark(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps centroid + smooth : " << duration / 10.0f << std::endl;
};

void test2(cgogn::CMap3* mrm, Attribute<Vec3>* attr)
{
	std::clock_t start;
	double duration;

	start = std::clock();
	for (int i = 0; i < 10; i++)
	{
		cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
			cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
				cgogn::foreach_incident_vertex(*mrm, w, [&](Vertex) -> bool { return true; });
				return true;
			});
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorA : " << duration / 10.0f << std::endl;

	start = std::clock();

	for (int i = 0; i < 10; i++)
	{
		cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
			cgogn::foreach_incident_vertex(*mrm, v, [&](Vertex w) -> bool {
				cgogn::foreach_incident_volume(*mrm, w, [&](Volume) -> bool { return true; });
				return true;
			});
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps CirculatorB : " << duration / 10.0f << std::endl;

	start = std::clock();
	for (int i = 0; i < 10; i++)
	{
		cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
			cgogn::geometry::centroid<Vec3>(*mrm, v, attr);
			return true;
		});
		cgogn::foreach_cell(*mrm, [&](cgogn::CMap3::Vertex v) -> bool {
			cgogn::CellMarkerStore<cgogn::CMap3, cgogn::CMap3::Vertex> mv(*mrm);
			Vec3 cm(0, 0, 0);
			cgogn::foreach_incident_volume(*mrm, v, [&](cgogn::CMap3::Volume w) -> bool {
				cgogn::foreach_incident_vertex(*mrm, w, [&](cgogn::CMap3::Vertex v2) -> bool {
					if (!mv.is_marked(v2))
					{
						cm += cgogn::value<Vec3>(*mrm, attr, v2);
						mv.mark(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "temps centroid + smooth : " << duration / 10.0f << std::endl;
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

	delete m2;
	delete mrm;
	delete m;

	return 0;
}
