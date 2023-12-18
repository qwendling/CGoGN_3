/**************	*****************************************************************
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

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/topo_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

// using Mesh = cgogn::IncidenceGraph;
using Mesh = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

class LocalInterface : public cgogn::ui::ViewModule
{

public:
	LocalInterface(const cgogn::ui::App& app)
		: cgogn::ui::ViewModule(app, "LocalInterface"), mesh_(nullptr), vertex_position_(nullptr),
		  mesh_provider_(nullptr), surface_render_(nullptr), topo_render_(nullptr),
		  moving_color_(1.0f, 0.0f, 1.0f, 1.0f)
	{
		view_ = app.current_view();
	}

	~LocalInterface()
	{
	}
	void force_update()
	{
		for (cgogn::ui::View* v : linked_views_)
			v->request_update();
	}

	void init() override
	{
		mesh_provider_ = static_cast<cgogn::ui::MeshProvider<Mesh>*>(
			app_.module("MeshProvider (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));

		surface_render_ = static_cast<cgogn::ui::SurfaceRender<Mesh>*>(
			app_.module("SurfaceRender (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));
		topo_render_ = static_cast<cgogn::ui::TopoRender<Mesh>*>(
			app_.module("TopoRender (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));
	}

	void left_panel() override
	{

		if (ImGui::Button("init moving"))
		{
			if (!moving_dart_.is_nil())
				topo_render_->set_dart_color(moving_dart_, moving_color_);
			srand(time(NULL));
			moving_dart_ = cgogn::Dart(933);
			/*cgogn::foreach_cell(*mesh_, [&](Face f) -> bool {
				int r = rand() % 1000;
				if (r < 10)
				{
					moving_dart_ = f.dart;
					return false;
				}
				return true;
			});*/
			cgogn::foreach_incident_edge(*mesh_, Vertex(moving_dart_), [&](Edge e) -> bool { return true; });
			topo_render_->set_dart_color(moving_dart_, moving_color_);
			force_update();
		}

		if (ImGui::Button("ph1"))
		{
			cgogn::Dart new_moving_dart_ = cgogn::phi1(*mesh_, moving_dart_);
			topo_render_->set_dart_color(new_moving_dart_, moving_color_);
			topo_render_->reset_dart_color(moving_dart_);
			moving_dart_ = new_moving_dart_;
			force_update();
		}

		if (ImGui::Button("phi2"))
		{
			cgogn::Dart new_moving_dart_ = cgogn::phi2(*mesh_, moving_dart_);
			topo_render_->set_dart_color(new_moving_dart_, moving_color_);
			topo_render_->reset_dart_color(moving_dart_);
			moving_dart_ = new_moving_dart_;
			force_update();
		}
		if (!moving_dart_.is_nil())
			ImGui::Text("Dart index : %d", moving_dart_.index);
	}
	Mesh* mesh_;
	cgogn::ui::View* view_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	cgogn::ui::MeshProvider<Mesh>* mesh_provider_;
	cgogn::ui::SurfaceRender<Mesh>* surface_render_;
	cgogn::ui::TopoRender<Mesh>* topo_render_;
	cgogn::Dart d_hexa_;
	cgogn::Dart d_pyra_;
	cgogn::Dart moving_dart_;
	Eigen::Vector4f moving_color_;
	float expl_vol_;
};

int main(int argc, char** argv)
{
	std::string filename;
	if (argc <= 1)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple surface viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);
	cgogn::ui::TopoRender<Mesh> tpr(app);
	LocalInterface interf(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&tpr);
	v1->link_module(&interf);

	if (filename.length() > 0)
	{
		Mesh* m = mp.load_surface_from_file(filename);
		if (!m)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}

		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
		std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");

		cgogn::foreach_cell(*m, [&](Vertex v) -> bool {
			cgogn::value<Vec3>(*m, vertex_position, v) += Vec3(1, 0, 0);
			return true;
		});
		cgogn::value<Vec3>(*m, vertex_position, Vertex(cgogn::Dart(0)));

		interf.mesh_ = m;
		interf.vertex_position_ = vertex_position;

		// std::shared_ptr<Attribute<Vec3>> face_color = cgogn::add_attribute<Vec3, Face>(*m, "color");
		// std::shared_ptr<Attribute<Scalar>> face_weight = cgogn::add_attribute<Scalar, Face>(*m, "weight");

		// cgogn::foreach_cell(*m, [&](Face f) -> bool {
		// 	Vec3 c(0, 0, 0);
		// 	c[rand() % 3] = 1;
		// 	cgogn::value<Vec3>(*m, face_color, f) = c;
		// 	cgogn::value<Scalar>(*m, face_weight, f) = double(rand()) / RAND_MAX;
		// 	return true;
		// });

		// mp.set_mesh_bb_vertex_position(*m, vertex_position);

		sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());

		sr.set_vertex_position(*v1, *m, vertex_position);
		sr.set_vertex_normal(*v1, *m, vertex_normal);
		tpr.set_vertex_position(*v1, *m, vertex_position);
	}

	return app.launch();
}
