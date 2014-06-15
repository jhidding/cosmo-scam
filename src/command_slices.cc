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

using namespace System;
using namespace Scam;
using namespace TwoMass;

std::function<void (Context)> prepare_context_slice(int w, int h, double L, bool rv)
{
	return [w,h,L,rv] (Context cx)
	{
		cx->translate(w/2, h/2);
		cx->scale(h/(1.1*L), h/(1.1*L));
		//cx->translate(-L/2, -L/2);
		cx->rectangle(-L/2,-L/2, L,L);
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

void command_slices(int argc_, char **argv_)
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
		Option({Option::VALUED | Option::CHECK, "fof", "haloes", "none",
			"include Halo catalog."}),
		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),
		Option({Option::VALUED | Option::CHECK, "L", "size", "100",
			"size of the box."}),
		Option({Option::VALUED | Option::CHECK, "t", "time", "1.0",
			"growing mode parameter."}),
		Option({Option::VALUED | Option::CHECK, "rs", "rstep", "10",
			"selection radius."}),
		Option({Option::VALUED | Option::CHECK, "dr", "dr", "10",
			"selection radius."}),
		Option({Option::VALUED | Option::CHECK, "f", "fila-lim", "200.0",
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
	auto polygons = read_polygons(fi, ply, v, format);
	std::cout << "read " << polygons.size() << " polygons.\n";

	std::cerr << "Reading " << fn_i_fila << " ..." << std::endl;
	v->clear(); ply = make_ptr<PLY::PLY>(); fi.close();
	
	std::ifstream fi2(fn_i_fila);
	format = PLY::read_header(fi2, ply);
	ply->print_header(std::cout, PLY::BINARY);
	v = read_vertices(fi2, ply, format);
	std::cout << "read " << v.size() << " vertices.\n";
	auto segments = read_segments(fi2, ply, v, format);
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

	Maybe<Array<Vertex>> haloes;
	if (argv["haloes"] != "none")
		haloes = Just(read_haloes(argv["haloes"]));
	else
		haloes = Nothing;

	double step = argv.get<double>("rstep");
	double dr = argv.get<double>("dr");

	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();

	auto cluster_label = scale_material(make_cluster_label_material(rv), L/(2*sqrt(2)));
	auto wall_material = make_wall_material(rv);
	auto filament_material = make_filament_material(rv);
	auto galaxy_material = scale_material(make_galaxy_material(rv), L/(2*sqrt(2)));
	auto halo_material = make_halo_material(rv);
	
	std::vector<double> slices;
	for (double r = 0.0; r+dr <=L; r +=step) slices.push_back(r);
		
	for (unsigned k = 0; k < 3; ++k) {
	//#pragma omp parallel for
	for (unsigned i = 0; i < slices.size(); ++i)
	{
		double r = slices[i];
		std::string dim_name = std::vector<std::string>({"x", "y", "z"})[k];
		std::string fn_output_png = timed_filename(argv["id"], 
			Misc::format("slice-", dim_name, "-", r, "-", r+dr), t, "png");

		Plane P1, P2; Point cam_pos, cam_tgt; Vector cam_shub;
		switch (k)
		{
			case 0: P1 = Plane(Point(   r, 0, 0), Vector( 1, 0, 0));
				P2 = Plane(Point(r+dr, 0, 0), Vector( 1, 0, 0)); 
				cam_pos = Point(L, L/2, L/2);
				cam_tgt = Point(L/2, L/2, L/2);
				cam_shub = Vector(0, 0, 1); 
				break;
			case 1: P1 = Plane(Point(0,    r, 0), Vector(0,  1, 0));
				P2 = Plane(Point(0, r+dr, 0), Vector(0,  1, 0));
				cam_pos = Point(L/2, L, L/2);
				cam_tgt = Point(L/2, L/2, L/2);
				cam_shub = Vector(0, 0, 1); 
				break;
			case 2: P1 = Plane(Point(0, 0,    r), Vector(0, 0,  1));
				P2 = Plane(Point(0, 0, r+dr), Vector(0, 0,  1));
				cam_pos = Point(L/2, L/2, L);
				cam_tgt = Point(L/2-0.001, L/2+0.001, L/2);
				cam_shub = Vector(-1, 0, 0); 
				break;
		}

		std::cout << "filtering polygons slice x_" << k << " = " << r << "\n";

		Array<Polygon> filtered_polygons;
		Array<Vertex>  filtered_galaxies;
		Array<Vertex>  filtered_haloes;
		Array<Segment> filtered_segments;
		Array<Vertex>  filtered_clusters;

		if (clusters)
		for (Vertex const &c : *clusters)
		{
			auto r_ = c.get_info<double>("r");
			double cr = 0.0; if (r_) cr = *r_;

			if ((not P1.is_below(c)) and P2.is_below(c))
				filtered_clusters.push_back(c);
		}

		for (Polygon const &p : polygons)
		{
			auto d = p.get_info<double>("density");
			if ((not d) or (*d < wall_lim)) continue;

			auto A = P1.split_polygon(p);
			if (not A.second) continue;
			auto B = P2.split_polygon(*A.second);
			if (not B.first) continue;

			filtered_polygons.push_back(*B.first);
		}

		if (galaxies)
		for (Vertex const &v : *galaxies)
		{
			if ((not P1.is_below(v)) and P2.is_below(v))
				filtered_galaxies.push_back(v);
		}
		if (haloes)
		for (Vertex const &v : *haloes)
		{
			if ((not P1.is_below(v)) and P2.is_below(v))
				filtered_haloes.push_back(v);
		}
		std::cerr << "including " << filtered_haloes.size() << " haloes\n";

		for (Segment const &s : segments)
		{
			auto d = s.get_info<double>("density");
			if ((not d) or (*d < fila_lim)) continue;
			auto A = P1.split_segment(s);
			if (not A.second) continue;
			auto B = P2.split_segment(*A.second);
			if (not B.first) continue;

			filtered_segments.push_back(*B.first);
		}
		std::cerr << "#segments: " << filtered_segments.size() << std::endl;

		Array<ptr<RenderObject>> scene;
		scene.push_back(ptr<RenderObject>(new PolygonObject(
			filtered_polygons, wall_material)));
		
		scene.push_back(ptr<RenderObject>(new SegmentObject(
			filtered_segments, filament_material)));

		if (galaxies)
		scene.push_back(ptr<RenderObject>(new VertexObject(
			filtered_galaxies, galaxy_material)));

		if (haloes)
		scene.push_back(ptr<RenderObject>(new VertexObject(
			filtered_haloes, halo_material)));

		scene.push_back(ptr<RenderObject>(new VertexObject(
			filtered_clusters, cluster_label)));
		
		auto C = make_ptr<Camera>(cam_pos, cam_tgt, cam_shub,
			parallel_projection);

		auto R = Renderer::Image(1920, 1080);
		R->apply(prepare_context_slice(1920, 1080, argv.get<double>("size"), rv));

		R->render(scene, C);
		R->write_to_png(fn_output_png);
		R->finish();
	} }
}

#include "base/global.hh"
System::Global<Command> _COMMAND_SLICES("slice", command_slices);

