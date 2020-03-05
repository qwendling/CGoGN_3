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

#include <cgogn/ui/modules/graph_render/graph_render.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>

#include <cgogn/modeling/algos/graph_to_hex.h>

int main(int argc, char** argv)
{
	using Graph = cgogn::Graph;
	using Surface = cgogn::CMap2;
	using Volume = cgogn::CMap3;

	std::string filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " graph_file" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Tubular mesh");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Graph> mpg(app);
	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Volume> mpv(app);

	cgogn::ui::GraphRender<Graph> gr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::SurfaceRender<Volume> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mpg);
	v1->link_module(&mps);
	v1->link_module(&mpv);

	v1->link_module(&gr);
	v1->link_module(&sr);
	v1->link_module(&vr);

	Graph* g = mpg.load_graph_from_file(filename);
	if (!g)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	Surface* s = mps.add_mesh("contact");
	Volume* v = mpv.add_mesh("tubes");

	if (cgogn::modeling::graph_to_hex(*g, *s, *v))
	{
		mps.emit_connectivity_changed(s);
		mpv.emit_connectivity_changed(v);
	}

	return app.launch();
}
