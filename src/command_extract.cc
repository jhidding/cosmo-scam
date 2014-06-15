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

#include <iomanip>

using namespace System;
using namespace Scam;
using namespace TwoMass;

std::function<void (Context)> prepare_context_extract(int w, int h, double L, bool rv)
{
	return [w,h,L,rv] (Context cx)
	{
		cx->translate(w/2, h/2);
		cx->scale(h/(1.1*L), h/(1.1*L));
		//cx->translate(-L/2, -L/2);
		//cx->rectangle(-L/2,-L/2, L,L);
		cx->arc(0,0,L/2,0,2*M_PI);
		if (rv) cx->set_source_rgba(1,1,1,0.5);
		else cx->set_source_rgba(0,0,0,0.5);
		cx->set_line_width(1);
		cx->stroke();
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

void command_tully(int argc_, char **argv_)
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

	// A1367	234.81	 73.03	0.0215 Leo
	// A1656	 58.09	 87.96	0.0232 Coma
	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();

	double rad = 40;
	auto cluster_label = scale_material(make_cluster_label_material(rv), rad/(sqrt(2)));
	auto wall_material = scale_material(make_wall_material(rv), rad/(sqrt(2)));
	auto filament_material = scale_material(make_filament_material(rv), rad/(sqrt(2)));
	auto galaxy_material = scale_material(make_galaxy_material(rv), rad*1.5/(sqrt(2)));
	
	Tully_filter cf(rad);

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

	Array<double> th;
	for (double t = 0; t <= 4*M_PI; t += 4*M_PI/40.) th.push_back(t);

	for (unsigned i = 0; i < th.size(); ++i)
	{	
		Point p = cf.center();
		Vector dp(cos(th[i]), sin(th[i]), sin(3./2*th[i]));
		auto C = make_ptr<Camera>(p + dp * 50, cf.center(), Vector(0,0,1),
			parallel_projection);
		auto R = Renderer::Image(1920, 1080);
		R->apply(prepare_context_extract(1920, 1080, rad*2, rv));
		R->render(scene, C);
		R->write_to_png(Misc::format("tully-b", time_string(t), "-", 
			std::setfill('0'), std::setw(3), i, ".png"));
		R->finish();
	}
}

void command_perseus(int argc_, char **argv_)
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

	std::string fn_i_wall = Misc::format(argv["id"], ".", time_string(t), ".walls.ply");
	std::string fn_i_fila = Misc::format(argv["id"], ".", time_string(t), ".filam.ply");

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
		clusters = Just(read_abell("perseus_catalog.tsv", true, rv));
	else
		clusters = Nothing;

	// A1367	234.81	 73.03	0.0215 Leo
	// A1656	 58.09	 87.96	0.0232 Coma
	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();

	double rad = 30;
	auto cluster_label = scale_material(make_cluster_label_material(rv), rad/(sqrt(2)));
	auto wall_material = scale_material(make_wall_material(rv), rad/(sqrt(2)));
	auto filament_material = scale_material(make_filament_material(rv), rad/(sqrt(2)));
	auto galaxy_material = scale_material(make_galaxy_material(rv), rad*1.5/(sqrt(2)));
	
	Perseus_filter cf(rad);

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

	Array<double> th;
	for (double t = 0; t <= 4*M_PI; t += 4*M_PI/1000.) th.push_back(t);

	for (unsigned i = 0; i < th.size(); ++i)
	{	
		Point p = cf.center();
		Vector dp(cos(th[i]), sin(th[i]), sin(3./2*th[i]));
		auto C = make_ptr<Camera>(p + dp * 50, cf.center(), Vector(0,0,1),
			parallel_projection);
		auto R = Renderer::Image(1920, 1080);
		R->apply(prepare_context_extract(1920, 1080, rad*2, rv));
		R->render(scene, C);
		R->write_to_png(Misc::format("perseus-b-", std::setfill('0'), std::setw(3), i, ".png"));
		R->finish();
	}
}

void command_coma(int argc_, char **argv_)
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

	std::string fn_i_wall = Misc::format(argv["id"], ".", time_string(t), ".walls.ply");
	std::string fn_i_fila = Misc::format(argv["id"], ".", time_string(t), ".filam.ply");

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
		clusters = Just(read_abell("coma_catalog.tsv", true, rv));
	else
		clusters = Nothing;

	// A1367	234.81	 73.03	0.0215 Leo
	// A1656	 58.09	 87.96	0.0232 Coma
	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();

	auto cluster_label = scale_material(make_cluster_label_material(rv), 100./(2*sqrt(2)));
	auto wall_material = scale_material(make_wall_material(rv), 100./(2*sqrt(2)));
	auto filament_material = scale_material(make_filament_material(rv), 100./(2*sqrt(2)));
	auto galaxy_material = scale_material(make_galaxy_material(rv), 140./(2*sqrt(2)));
	
	Coma_filter cf(30, 50);

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

	Array<double> th;
	for (double t = 0; t <= 2*M_PI; t += 2*M_PI/250.) th.push_back(t);

	for (unsigned i = 0; i < th.size(); ++i)
	{	
		Point p = cf.center();
		Vector dp(cos(th[i]), sin(th[i]), 0);
		auto C = make_ptr<Camera>(p + dp * 50, cf.center(), -cf.shub(),
			parallel_projection);
		auto R = Renderer::Image(1920, 1080);
		R->apply(prepare_context_extract(1920, 1080, 100, rv));
		R->render(scene, C);
		R->write_to_png(Misc::format("coma-b-", std::setfill('0'), std::setw(3), i, ".png"));
		R->finish();
	}
}

#include "base/global.hh"
System::Global<Command> _COMMAND_COMA("coma", command_coma);
System::Global<Command> _COMMAND_PERSEUS("perseus", command_perseus);
System::Global<Command> _COMMAND_TULLY("tully", command_tully);

