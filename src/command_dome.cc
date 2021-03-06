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

#include "support.hh"
#include "cosmic.hh"
#include "materials.hh"
#include "coma_filter.hh"
#include "make_coords.hh"

#include <iomanip>

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

inline Material scale_material(Material M, double scale)
{
	return [M, scale] (Info info, Context cx)
	{
		cx->save();
		cx->scale(scale, scale);
		M(info, cx);
		cx->restore();
	};
}

void command_dome(int argc_, char **argv_)
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

	// A1367	234.81	 73.03	0.0215 Leo
	// A1656	 58.09	 87.96	0.0232 Coma
	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();
	unsigned N_imag = argv.get<unsigned>("frames");

	double rad = 17;
	double mat_scale = 1.0;
	auto cluster_label = scale_material(make_cluster_label_material(rv), mat_scale);
	auto wall_material = scale_material(make_wall_material(rv), mat_scale);
	auto filament_material = scale_material(make_filament_material(rv), mat_scale);
	auto galaxy_material = scale_material(make_galaxy_material(rv), mat_scale*1.5);
	auto velocity_material = scale_material(make_vel_material(rv, argv.get<double>("vel-width")), mat_scale);
	
	Home_filter cf(rad);

	Array<ptr<RenderObject>> scene;
	scene.push_back(ptr<RenderObject>(new PolygonObject(
		cf(polygons), wall_material)));
	scene.push_back(ptr<RenderObject>(new SegmentObject(
		cf(segments), filament_material)));
	scene.push_back(ptr<RenderObject>(new SegmentObject(
		cf.fog(), [] (Info info, Context cx)
	{
		cx->set_line_width(0.1);
		cx->set_source_rgb(0.5,0.7,1.0);
		cx->stroke();
	})));
	if (galaxies)
	scene.push_back(ptr<RenderObject>(new VertexObject(
		cf(*galaxies), galaxy_material)));
	if (clusters)
	scene.push_back(ptr<RenderObject>(new VertexObject(
		cf(*clusters), cluster_label)));
	if (velocities)
	scene.push_back(ptr<RenderObject>(new VectorObject(
		cf(*velocities), velocity_material, argv.get<double>("vel-factor"))));

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


#include "base/global.hh"
System::Global<Command> _COMMAND_DOME("dome", command_dome);

