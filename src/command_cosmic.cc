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

void command_cosmic(int argc_, char **argv_)
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

	/*
	Array<Cluster> cluster_catalog = {
		{ 186.75,  12.72,  10.00, "Virgo cluster",    0.0,  0.3 },

		{  49.50,  41.50,  53.66, "Perseus cluster", -0.2, -0.2 },
		{  28.20,  36.15,  45.30, "A262",             0.1,  0.1 },
		{  36.45,  41.87,  51.60, "A347",             0.2,  -0.03 },

		{ 243.89, -60.91,  47.07, "Norma cluster",    0.2, -0.2 },

		{ 194.95,  27.98,  69.25, "Coma cluster",     0.0,  0.3 }, 
		{ 176.12,  19.84,  65.95, "Leo cluster",      0.0, -0.3 } };
	*/

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

	Maybe<Array<Vertex>> clusters;
	if (argv.get<bool>("abell"))
		clusters = Just(read_abell("abell_catalog.tsv", false, argv.get<bool>("reverse")));
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

	bool rv = argv.get<bool>("reverse");
	auto cluster_label = make_cluster_label_material(rv);
	auto wall_material = make_wall_material(rv);
	auto filament_material = make_filament_material(rv);
	auto galaxy_material = make_galaxy_material(rv);
	auto halo_material = make_halo_material(rv);

	for (double r = 0.0; r + dr < L/2; r += step)
	{
		//std::string fn_output = timed_filename(argv["id"], Misc::format("slice-", r, "-", r+dr), t, "pdf");
		std::string fn_output_png = timed_filename(argv["id"], 
			Misc::format("slice-", r, "-", r+dr), t, "png");
		Sphere S1(Point(L/2,L/2,L/2), r), S2(Point(L/2,L/2,L/2), r + dr);
		std::cout << "filtering polygons ... radius " << r << "\n";
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

			if (cr < r-3) continue;
			if (cr > r+dr+3) continue;

			filtered_clusters.push_back(c);
		}

		for (Polygon const &p : polygons)
		{
			auto d = p.get_info<double>("density");
			if ((not d) or (*d < wall_lim)) continue;
			auto A = S1.split_polygon(p);
			if (not A.second) continue;
			auto B = S2.split_polygon(*A.second);
			if (not B.first) continue;

			filtered_polygons.push_back(*B.first);
		}

		if (galaxies)
		for (Vertex const &v : *galaxies)
		{
			if ((not S1.is_below(v)) and S2.is_below(v))
				filtered_galaxies.push_back(v);
		}
		if (haloes)
		for (Vertex const &v : *haloes)
		{
			if ((not S1.is_below(v)) and S2.is_below(v))
				filtered_haloes.push_back(v);
		}
		std::cerr << "including " << filtered_haloes.size() << " haloes\n";

		for (Segment const &s : segments)
		{
			auto d = s.get_info<double>("density");
			if ((not d) or (*d < fila_lim)) continue;
			auto A = S1.split_segment(s);
			if (not A.second) continue;
			auto B = S2.split_segment(*A.second);
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
		
		auto C = make_ptr<Map_projection_camera>(
			Point(L/2,L/2,L/2), Point(L/2, L/2+0.001, L), Vector(1, 0, 0),
			//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
				Map_projection(Aitoff_Hammer));

		auto R = Renderer::Image(1920, 1080);
		if (argv.get<bool>("reverse"))
			R->apply(prepare_context_rv(1920, 1080, r+dr/2));
		else
			R->apply(prepare_context(1920, 1080, r+dr/2));

		R->render(scene, C);
		R->write_to_png(fn_output_png);
		R->finish();
	}
}

#include "base/global.hh"
System::Global<Command> _COMMAND_COSMIC("cosmic", command_cosmic);

