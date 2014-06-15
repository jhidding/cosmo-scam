#include "base/common.hh"
#include "base/argv.hh"
#include "base/format.hh"
#include "ply/ply.hh"
#include "ply/read.hh"
#include "geometry/geometry.hh"

#include "support.hh"
#include "cosmic.hh"

using namespace System;
using namespace Scam;

void write_walls_to_ply(std::string const &fn, Array<Polygon> walls)
{
	PLY::PLY ply;
	ply.add_comment("walls");
	ply.add_element("vertex", 
		PLY::scalar_type<float>("x"), 
		PLY::scalar_type<float>("y"), 
		PLY::scalar_type<float>("z"));

	std::map<size_t,unsigned> vertex_ids; unsigned vi = 0;
	for (Polygon const &p : walls)
	{
		for (Vertex const &v : p)
		{
			size_t vid = std::hash<Vertex>()(v);
			if (vertex_ids.count(vid) == 0)
			{
				ply.put_data(
					PLY::scalar<float>(v.x()),
					PLY::scalar<float>(v.y()),
					PLY::scalar<float>(v.z()));

				vertex_ids[vid] = vi;
				++vi;
			}
		}
	}

	ply.add_element("face",
		PLY::list_type<unsigned>("vertex_indices"),
		PLY::scalar_type<float>("density"));

	for (Polygon const &p : walls)
	{
		std::vector<unsigned> verts;
		double density = *p.get_info<double>("density");

		for (Vertex const &v : p)
		{
			size_t vid = std::hash<Vertex>()(v);
			verts.push_back(vertex_ids[vid]);
		}

		ply.put_data(
			PLY::list<unsigned>(verts),
			PLY::scalar<float>(density));
	}

	ply.write(fn, PLY::BINARY);
}

void write_filam_to_ply(std::string const &fn, Array<Segment> filam)
{
	PLY::PLY ply;
	ply.add_comment("filaments");
	ply.add_element("vertex", 
		PLY::scalar_type<float>("x"), 
		PLY::scalar_type<float>("y"), 
		PLY::scalar_type<float>("z"));

	std::map<size_t,unsigned> vertex_ids; unsigned vi = 0;
	for (Segment const &p : filam)
	{
		auto f = [&] (Vertex const &v)
		{
			size_t vid = std::hash<Vertex>()(v);
			if (vertex_ids.count(vid) == 0)
			{
				ply.put_data(
					PLY::scalar<float>(v.x()),
					PLY::scalar<float>(v.y()),
					PLY::scalar<float>(v.z()));

				vertex_ids[vid] = vi;
				++vi;
			}
		};

		f(p.first()); f(p.second());
	}

	ply.add_element("edge",
		PLY::scalar_type<unsigned>("vertex1"),
		PLY::scalar_type<unsigned>("vertex2"),
		PLY::scalar_type<float>("density"));

	for (Segment const &p : filam)
	{
		double density = *p.get_info<double>("density");

		size_t vid1 = std::hash<Vertex>()(p.first()),
		       vid2 = std::hash<Vertex>()(p.second());

		ply.put_data(
			PLY::scalar<unsigned>(vertex_ids[vid1]),
			PLY::scalar<unsigned>(vertex_ids[vid2]),
			PLY::scalar<float>(density));
	}

	ply.write(fn, PLY::BINARY);
}

void command_split(int argc_, char **argv_)
{
	Argv argv = read_arguments(argc_, argv_,
		Option({0, "h", "help", "false", 
			"print this help."}),
		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),
		Option({Option::VALUED | Option::CHECK, "L", "size", "180",
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

	double step = argv.get<double>("rstep");
	double dr = argv.get<double>("dr");

	std::string fn_i_wall = Misc::format(argv["id"], ".", time_string(t), ".walls.ply");
	std::string fn_i_fila = Misc::format(argv["id"], ".", time_string(t), ".filam.ply");

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

	unsigned i = 0;
	for (double r = 0.0; r + dr <= L/2; r += step)
	{
		std::string fn_output_png = timed_filename(argv["id"], 
			Misc::format("slice-", r, "-", r+dr), t, "png");
		Sphere S1(Point(L/2,L/2,L/2), r), S2(Point(L/2,L/2,L/2), r + dr);
		std::cout << "filtering polygons ... radius " << r << "\n";
		Array<Polygon> filtered_polygons;
		Array<Segment> filtered_segments;

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
		std::cerr << "#polygons: " << filtered_polygons.size() << std::endl;

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

		std::string fn_o_wall = Misc::format(
			argv["id"], ".", time_string(t), ".walls.", i ,".ply");
		std::string fn_o_fila = Misc::format(
			argv["id"], ".", time_string(t), ".filam.", i, ".ply");

		write_filam_to_ply(fn_o_fila, filtered_segments);
		write_walls_to_ply(fn_o_wall, filtered_polygons);

		++i;
	}
}

#include "base/global.hh"
System::Global<Command> _COMMAND_SPLIT("split", command_split);

