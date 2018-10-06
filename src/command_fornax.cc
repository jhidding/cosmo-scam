#include "base/common.hh"
#include "base/argv.hh"
#include "base/format.hh"
#include "ply/ply.hh"
#include "ply/read.hh"
#include "geometry/geometry.hh"
#include "render/render.hh"
#include "render/map_projection.hh"
#include "two_mass.hh"
#include "material/colour.hh"
#include "tex/tex.hh"
#include "pnm.hh"

#include "support.hh"
#include "cosmic.hh"
#include "materials.hh"
#include "coma_filter.hh"
#include "make_coords.hh"

#include "base/mvector.hh"

using namespace System;
using namespace Scam;
using namespace TwoMass;

inline std::function<void (Context)> prepare_context_extract(int w, int h, double L, bool rv)
{
	return [w,h,L,rv] (Context cx)
	{
		cx->translate(w/2, h/2);
		cx->scale(h/L, h/L);
	};
}

void command_fornax(int argc_, char **argv_)
{
	Argv argv = read_arguments(argc_, argv_,
		Option({0, "h", "help", "false",
			"print this help."}),
		Option({0, "2m", "2mass", "false",
			"include 2Mass into rendering."}),
		Option({0, "ac", "abell", "false",
			"include Abell cluster catalog."}),
		Option({0, "r", "reverse", "false",
			"reverse colours."}),
		Option({0, "v", "velocity", "false",
			"include velocities."}),
		Option({0, "sky", "sky", "false",
			"sky plot"}),
		Option({Option::VALUED | Option::CHECK, "rx", "resx", "1920",
			"image size X"}),
		Option({Option::VALUED | Option::CHECK, "ry", "resy", "1080",
			"image size Y"}),
		Option({Option::VALUED | Option::CHECK, "vf", "vel-factor", "0.01",
			"velocity vector factor."}),
		Option({Option::VALUED | Option::CHECK, "vw", "vel-width", "1",
			"velocity vector width."}),
		Option({Option::VALUED | Option::CHECK, "nf", "frames", "999",
			"number of frames in the movie."}),
		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),
		Option({Option::VALUED | Option::CHECK, "L", "size", "180",
			"size of the box."}),
		Option({Option::VALUED | Option::CHECK, "t", "time", "1.0",
			"growing mode parameter."}),
		Option({Option::VALUED | Option::CHECK, "f", "fila-lim", "1e6",
			"lower limit of filament density to show."}),
		Option({Option::VALUED | Option::CHECK, "w", "wall-lim", "20.0",
			"lower limit of wall density to show."}) );

	if (argv.get<bool>("help"))
	{
		std::cerr << "Scam -- 3D vector graphics for science.\n"
			"Copyright Johan Hidding, June 2014 - licence: GPL3.\n\n";
		argv.print(std::cerr);
		exit(0);
	}

	double t = argv.get<double>("time");
	double L = argv.get<double>("size");
	double wall_lim = argv.get<double>("wall-lim");
	double fila_lim = argv.get<double>("fila-lim");

	std::string fn_i_wall = Misc::format(argv["id"], ".walls.", time_string(t), ".ply");
	std::string fn_i_fila = Misc::format(argv["id"], ".filam.", time_string(t), ".ply");

	std::cerr << "Reading " << fn_i_wall << " ..." << std::endl;
	auto ply = make_ptr<PLY::PLY>();

	std::ifstream fi(fn_i_wall);
	PLY::Format format = PLY::read_header(fi, ply);
	ply->print_header(std::cout, PLY::BINARY);
	auto v = read_vertices(fi, ply, format);
	std::cout << "read " << v.size() << " vertices.\n";
	auto polygons = read_polygons(fi, ply, v, format, argv.get<double>("wall-lim"));
	std::cout << "read " << polygons.size() << " polygons.\n";

	std::cerr << "Reading " << fn_i_fila << " ..." << std::endl;
	v->clear(); ply = make_ptr<PLY::PLY>(); fi.close();

	std::ifstream fi2(fn_i_fila);
	format = PLY::read_header(fi2, ply);
	ply->print_header(std::cout, PLY::BINARY);
	v = read_vertices(fi2, ply, format);
	std::cout << "read " << v.size() << " vertices.\n";
	auto segments = read_segments(fi2, ply, v, format, argv.get<double>("fila-lim"));
	std::cout << "read " << segments.size() << " segments.\n";
	v->clear(); ply.reset(); fi2.close();

	Maybe<Array<Vertex>> galaxies;
	if (argv.get<bool>("2mass"))
		galaxies = Just(read_2mass("2mrs_1175_done.dat"));
	else
		galaxies = Nothing;

	bool rv = argv.get<bool>("reverse");
	Maybe<Array<Vertex>> clusters;
	if (argv.get<bool>("abell"))
		clusters = Just(read_abell("abell_catalog.tsv", true, rv));
	else
		clusters = Nothing;

	std::string fn_i_vel = Misc::format(argv["id"], ".nodes.", time_string(t), ".conan");
	Maybe<Array<Vertex>> velocities;
	if (argv.get<bool>("velocity"))
		velocities = Just(read_velocities(fn_i_vel));
	else
		velocities = Nothing;

	auto fornax_galaxies = read_fornax("fornax.txt");
	// A1367	234.81	 73.03	0.0215 Leo
	// A1656	 58.09	 87.96	0.0232 Coma
	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();
	unsigned N_imag = argv.get<unsigned>("frames");

	double rad = 20;
	double mat_scale = 1.0;
	auto cluster_label = scale_material(make_cluster_label_material(rv), mat_scale);
	// auto wall_material = scale_material(make_rainbow_material(rv, 10, 20), mat_scale);
	auto wall_material = scale_material(make_wall_material(rv),mat_scale);
	auto filament_material = scale_material(make_filament_material(rv), mat_scale);
	auto galaxy_material1 = scale_material(make_galaxy_material(rv), mat_scale / 25);
	auto galaxy_material2 = scale_material(make_galaxy_material(not rv), mat_scale / 5);
	auto velocity_material = scale_material(make_vel_material(rv, argv.get<double>("vel-width")), mat_scale);

	Fornax_filter cf(rad);
	Shell_filter S(10, 20);

	Array<ptr<RenderObject>> scene;
	scene.push_back(ptr<RenderObject>(new PolygonObject(
		S(cf(polygons)), wall_material)));
	scene.push_back(ptr<RenderObject>(new SegmentObject(
		S(cf(segments)), filament_material)));
	if (galaxies) {
	scene.push_back(ptr<RenderObject>(new VertexObject(
		fornax_galaxies, galaxy_material1)));
	scene.push_back(ptr<RenderObject>(new VertexObject(
		S(cf(*galaxies)), galaxy_material2)));
	}
	if (clusters)
	scene.push_back(ptr<RenderObject>(new VertexObject(
		S(cf(*clusters)), cluster_label)));
	if (velocities)
	scene.push_back(ptr<RenderObject>(new VectorObject(
		S(cf(*velocities)), velocity_material, argv.get<double>("vel-factor"))));

	if (argv.get<bool>("sky"))
	{
		Spherical_rotation side = eq_to_sg; //(0,-M_PI/2-0.0001,M_PI);
		Spherical_rotation sph_id(0, M_PI/2-1e-5, 0);

		/* add coordinate lines */
		add_coordinates(scene, side);
		scene.push_back(ptr<RenderObject>(new SegmentObject(
			make_parallel(side*ga_to_eq, 5, 1e5), [] (Info info, Context cx)
		{
					cx->set_source_rgb(0.5, 0.5, 0.5);
					cx->set_line_width(0.002);
					cx->stroke();
		})));
		scene.push_back(ptr<RenderObject>(new SegmentObject(
			make_parallel(side*ga_to_eq, -5, 1e5), [] (Info info, Context cx)
		{
					cx->set_source_rgb(0.5, 0.5, 0.5);
					cx->set_line_width(0.002);
					cx->stroke();
		})));

		Point pointing = cf.center();
		Vector shub = -cf.shub();

		auto C = make_ptr<Map_projection_camera>(
			Point(L/2,L/2,L/2), pointing, shub,
			//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
				Map_projection(Aitoff_Hammer));

		unsigned rx = argv.get<unsigned>("resx"),
			 ry = argv.get<unsigned>("resy");
		auto R = Renderer::Image(rx, ry);
		R->apply([rx, ry] (Context cx)
		{
			double f = 5.0;
			cx->scale(ry*0.5*f, ry*0.5*f);
			cx->translate((1.0/f * rx) / ry, 1.0/f);
		});
		R->render(scene, C);
		R->write_to_png("fornax-sky.png");
		R->finish();
		return;
	} else {
		double N_frames = argv.get<double>("frames");
		Array<double> th;
		for (double t = 0; t <= 2*M_PI; t += 2*M_PI/N_frames) th.push_back(t);

		for (unsigned i = 0; i < th.size(); ++i)
		{
			Point p = cf.center();
			Vector dp = Vector(cos(th[i]), sin(th[i]), sin(th[i])/5).normalize();
			auto C = make_ptr<Camera>(p + dp * (rad * 1.2), cf.center(), Vector(0,0,1),
				dome_projection);
			auto R = Renderer::Image(4096, 4096);
			R->apply(prepare_context_extract(4096, 4096, 3.1415927, rv));
			R->render(scene, C);
			R->write_to_png(Misc::format("dome", time_string(t), "-",
				std::setfill('0'), std::setw(4), i, ".png"));
			R->finish();
		}
	}
}


#include "base/global.hh"
System::Global<Command> _COMMAND_FORNAX("fornax", command_fornax);
