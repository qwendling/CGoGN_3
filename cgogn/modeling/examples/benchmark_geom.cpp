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

#include <algorithm>
#include <random>

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
using uint32 = cgogn::numerics::uint32;

using Vec3 = cgogn::geometry::Vec3;

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

std::shared_ptr<Attribute<Vec3>> extract_cmap(cgogn::CMap3* m, EMR_Map3& mrm, Attribute<Vec3>* pos)
{
	std::unordered_map<Dart, Dart> um;
	// Dart creation
	for (Dart d = mrm.begin(), end = mrm.end(); d != end; d = mrm.next(d))
	{
		Dart d2 = cgogn::add_dart(*m);
		um.insert({d, d2});
	}
	// Topology
	for (Dart d = mrm.begin(), end = mrm.end(); d != end; d = mrm.next(d))
	{
		Dart d2 = um[d];
		(*(m->phi1_))[d2.index] = cgogn::phi1(mrm, d);
		(*(m->phi_1_))[d2.index] = cgogn::phi_1(mrm, d);
		(*(m->phi2_))[d2.index] = cgogn::phi2(mrm, d);
		(*(m->phi3_))[d2.index] = cgogn::phi3(mrm, d);
	}
	cgogn::index_cells<cgogn::CMap3::Vertex>(*m);
	cgogn::index_cells<cgogn::CMap3::Volume>(*m);
	std::shared_ptr<Attribute<Vec3>> pos_cmap = cgogn::add_attribute<Vec3, Vertex>(*m, "pos_cmap");
	// geometry
	cgogn::foreach_cell(mrm, [&](EMR_Map3::Vertex v) -> bool {
		cgogn::value<Vec3>(*m, pos_cmap, cgogn::CMap3::Vertex(um[v.dart])) = cgogn::value<Vec3>(mrm, pos, v);
		return true;
	});
	return pos_cmap;
}

#define PERCENT_MODIF 10

int main(int argc, char** argv)
{
	int percent_nb_modif = 10;
	std::string filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " volume_mesh_file [Percent modif topo]" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);
	if (argc == 3)
		percent_nb_modif = atoi(argv[2]);

	cgogn::thread_start();

	std::string ext = cgogn::extension(filename);
	Mesh* m = new Mesh();
	bool imported = false;

	if (ext.compare("tet") == 0)
	{
		imported = cgogn::io::import_TET(*m, filename);
	}
	else if (ext.compare("mesh") == 0 || ext.compare("meshb") == 0)
	{
		imported = cgogn::io::import_MESHB(*m, filename);
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

	std::srand(2124512438);

	std::vector<Volume> volume_to_subdivided;
	std::vector<Volume> volume_to_simplified;

	cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
		volume_to_subdivided.push_back(Volume(mrm->volume_youngest_dart(v.dart)));
		return true;
	});

	std::cout << "debut benchmark" << std::endl;

	std::vector<Volume> choix_volume;
	std::vector<Volume> vect_vol;
	std::vector<Vertex> vect_vertex;
	std::vector<Volume> vect_vol2;
	std::vector<Vertex> vect_vertex2;
	int nb_cell = cgogn::nb_cells<Volume>(*mrm);
	bool test_dis = false;
	std::cout << "nb_volume;nb subdivision;subdivide mr;nb simplification;simplified mr;centroid+smoothing MR;centroid "
				 "+ smoothing mono"
			  << std::endl;
	for (int i = 0; i < 100; i++)
	{
		// std::cout << "Volume mono :" << cgogn::nb_cells<Volume>(*cph) << "Volume MR :" <<
		// cgogn::nb_cells<Volume>(*mrm) << std::endl; std::cout << "Face mono :" << cgogn::nb_cells<Face>(*cph) <<
		// "Face MR :" << cgogn::nb_cells<Face>(*mrm) << std::endl; std::cout << "Edge mono :" <<
		// cgogn::nb_cells<Edge>(*cph) << "Edge MR :" << cgogn::nb_cells<Edge>(*mrm) << std::endl; std::cout << "Vertex
		// mono :" << cgogn::nb_cells<Vertex>(*cph) << "Vertex MR :" << cgogn::nb_cells<Vertex>(*mrm) << std::endl;
		std::cout << nb_cell << ";";
		for (auto it = volume_to_subdivided.begin(); it != volume_to_subdivided.end();)
		{
			if (rand() % 100 < percent_nb_modif)
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
		std::cout << choix_volume.size() << ";";
		start = std::clock();
		for (Volume v : choix_volume)
		{
			mrm->activate_volume_subdivision_fast(v);
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";

		choix_volume.clear();

		for (auto it = volume_to_simplified.begin(); it != volume_to_simplified.end();)
		{
			if (test_dis)
				break;
			if (rand() % 100 < percent_nb_modif)
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
		std::cout << choix_volume.size() << ";";
		start = std::clock();
		for (Volume v : choix_volume)
		{
			mrm->disable_volume_subdivision(Volume(cgogn::phi1(*mrm, v.dart)), true);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration;
		// std::cout << "temps simplified 10% volume mr : " << duration << std::endl;

		choix_volume.clear();
		start = std::clock();

		cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
			cgogn::geometry::centroid<Vec3>(*mrm, v, position.get());
			return true;
		});
		cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
			cgogn::CellMarkerStore<MRMesh, Vertex> mv(*mrm);
			Vec3 cm(0, 0, 0);
			cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
				cgogn::foreach_incident_vertex(*mrm, w, [&](Vertex v2) -> bool {
					if (!mv.is_marked(v2))
					{
						cm += cgogn::value<Vec3>(*mrm, position.get(), v2);
						mv.mark(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";

		cgogn::CMap3* tmp_map = new cgogn::CMap3();
		std::shared_ptr<Attribute<Vec3>> pos_cmap = extract_cmap(tmp_map, *mrm, position.get());

		start = std::clock();

		cgogn::foreach_cell(*tmp_map, [&](Volume v) -> bool {
			cgogn::geometry::centroid<Vec3>(*tmp_map, v, pos_cmap.get());
			return true;
		});
		cgogn::foreach_cell(*tmp_map, [&](cgogn::CMap3::Vertex v) -> bool {
			cgogn::CellMarkerStore<cgogn::CMap3, cgogn::CMap3::Vertex> mv(*tmp_map);
			Vec3 cm(0, 0, 0);
			cgogn::foreach_incident_volume(*tmp_map, v, [&](cgogn::CMap3::Volume w) -> bool {
				cgogn::foreach_incident_vertex(*tmp_map, w, [&](cgogn::CMap3::Vertex v2) -> bool {
					if (!mv.is_marked(v2))
					{
						cm += cgogn::value<Vec3>(*tmp_map, pos_cmap.get(), v2);
						mv.mark(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration;

		std::cout << std::endl;
		delete tmp_map;
	}

	std::cout << "bench volume constant" << std::endl;
	std::cout << "nb_volume;nb subdivision;subdivide mr;nb simplification;simplified mr;centroid+smoothing MR;centroid "
				 "+ smoothing mono"
			  << std::endl;

	int nb_modif = volume_to_subdivided.size() * 0.1;
	std::random_device rd;
	std::mt19937 g(19111996);
	for (int i = 0; i < 100; i++)
	{
		// std::clog << cph->current_level_ << ";" << cph->maximum_level_ << std::endl;
		std::cout << nb_cell << ";";

		int j;
		std::shuffle(volume_to_subdivided.begin(), volume_to_subdivided.end(), g);
		for (j = 0; j < nb_modif; ++j)
		{
			Volume v = volume_to_subdivided.back();
			choix_volume.push_back(v);
			volume_to_simplified.push_back(v);
			volume_to_subdivided.pop_back();
			nb_cell += 7;
		}
		std::cout << choix_volume.size() << ";";
		start = std::clock();
		for (Volume v : choix_volume)
		{
			mrm->activate_volume_subdivision_fast(v);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";
		// std::cout << "temps subdivided 10% volume mr : " << duration << std::endl;

		choix_volume.clear();
		// std::shuffle(volume_to_simplified.begin(), volume_to_simplified.end(), g);

		for (j = 0; j < nb_modif; ++j)
		{
			Volume v = volume_to_simplified.back();
			choix_volume.push_back(v);
			volume_to_subdivided.push_back(v);
			volume_to_simplified.pop_back();
			nb_cell -= 7;
		}
		std::cout << choix_volume.size() << ";";
		start = std::clock();
		for (Volume v : choix_volume)
		{
			mrm->disable_volume_subdivision(v, true);
		}
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration;
		// std::cout << "temps simplified 10% volume mr : " << duration << std::endl;

		choix_volume.clear();
		start = std::clock();

		cgogn::foreach_cell(*mrm, [&](Volume v) -> bool {
			cgogn::geometry::centroid<Vec3>(*mrm, v, position.get());
			return true;
		});
		cgogn::foreach_cell(*mrm, [&](Vertex v) -> bool {
			cgogn::CellMarkerStore<MRMesh, Vertex> mv(*mrm);
			Vec3 cm(0, 0, 0);
			cgogn::foreach_incident_volume(*mrm, v, [&](Volume w) -> bool {
				cgogn::foreach_incident_vertex(*mrm, w, [&](Vertex v2) -> bool {
					if (!mv.is_marked(v2))
					{
						cm += cgogn::value<Vec3>(*mrm, position.get(), v2);
						mv.mark(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration << ";";

		cgogn::CMap3* tmp_map = new cgogn::CMap3();
		std::shared_ptr<Attribute<Vec3>> pos_cmap = extract_cmap(tmp_map, *mrm, position.get());

		start = std::clock();

		cgogn::foreach_cell(*tmp_map, [&](Volume v) -> bool {
			cgogn::geometry::centroid<Vec3>(*tmp_map, v, pos_cmap.get());
			return true;
		});
		cgogn::foreach_cell(*tmp_map, [&](cgogn::CMap3::Vertex v) -> bool {
			cgogn::CellMarkerStore<cgogn::CMap3, cgogn::CMap3::Vertex> mv(*tmp_map);
			Vec3 cm(0, 0, 0);
			cgogn::foreach_incident_volume(*tmp_map, v, [&](cgogn::CMap3::Volume w) -> bool {
				cgogn::foreach_incident_vertex(*tmp_map, w, [&](cgogn::CMap3::Vertex v2) -> bool {
					if (!mv.is_marked(v2))
					{
						cm += cgogn::value<Vec3>(*tmp_map, pos_cmap.get(), v2);
						mv.mark(v2);
					}
					return true;
				});
				return true;
			});
			return true;
		});

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << duration;

		std::cout << std::endl;
		delete tmp_map;
	}

	delete mrm;
	delete m;

	return 0;
}
