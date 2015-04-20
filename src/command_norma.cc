#include "pnm.hh"
#include <iomanip>

using namespace System;
using namespace Scam;
using namespace TwoMass;

std::function<void (Context)> prepare_context_extract_norma(int w, int h, double L, bool rv)
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

class Norma_filter
{
	// A3627	325.34	 -7.26	0.0143 Norma
	// A1060	269.64	 26.51	0.0114 Hydra
	// A3526	302.42	 21.56	0.0110 Centaurus

	Plane 	P, Top, Bottom;
	Sphere 	S;
	Vector  hub, m_normal;

	public:                 
		Norma_filter(double width, double radius, double offset = 0.0)
		{
			double ra, dec;
			//ga_to_sg(radians(191.07), radians(44.41), ra, dec);
			//Point A779  = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0266 * 2997.92458);
			ga_to_sg(radians(325.34), radians(-7.26), ra, dec);
			Point A3627 = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0143 * 2997.92458);
			ga_to_sg(radians(269.64), radians(26.51), ra, dec);
			Point A1060 = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0114 * 2997.92458);
			ga_to_sg( radians(302.42), radians(21.56), ra, dec);
			Point A3526 = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0110 * 2997.92458);
			Point earth = Point(90,90,90);
			
			P = Plane(A3627 + (A3526-A3627)*0.5, 
				Vector::cross(A3526 - A3627, earth - A3627).normalize());
			Top = Plane(P.origin() + P.normal() * (width/2 + offset), P.normal());
			Bottom = Plane(P.origin() - P.normal() * (width/2 - offset), -P.normal());
			S = Sphere(P.origin(), radius);

			hub = P.normal();
			m_normal = P.normal();
		}

		Point center() const { return P.origin(); }
		Vector shub() const { return hub; } 
		Vector normal() const { return m_normal; }

		Array<Vertex> operator()(Array<Vertex> A) const
		{
			Array<Vertex> B;
			for (Vertex const &a : A)
			{
				if (S.is_below(a) and Top.is_below(a) and Bottom.is_below(a))
					B.push_back(a);
			}
			return B;
		}

		Array<Polygon> operator()(Array<Polygon> A) const
		{
			Array<Polygon> B;

			for (Polygon const &a : A)
			{
				auto T = S.split_polygon(a);
					if (not T.first) continue;
				auto U = Top.split_polygon(*T.first);
					if (not U.first) continue;
				auto V = Bottom.split_polygon(*U.first);
					if (not V.first) continue;
				B.push_back(*V.first);
			}

			return B;
		}

		Array<Segment> operator()(Array<Segment> A) const
		{
			Array<Segment> B;

			for (Segment const &a : A)
			{
				auto T = S.split_segment(a);
					if (not T.first) continue;
				auto U = Top.split_segment(*T.first);
					if (not U.first) continue;
				auto V = Bottom.split_segment(*U.first);
					if (not V.first) continue;
				B.push_back(*V.first);
			}

			return B;
		}
};

void command_norma(int argc_, char **argv_)
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
		Option({0, "side", "side", "false",
			"side plot"}),
		Option({0, "slices", "slices", "false",
			"sliced plot"}),
		Option({Option::VALUED | Option::CHECK, "off", "offset", "0",
			"slice offset"}),

		Option({Option::VALUED | Option::CHECK, "rx", "resx", "1920",
			"image size X"}),
		Option({Option::VALUED | Option::CHECK, "ry", "resy", "1080",
			"image size Y"}),
		Option({Option::VALUED | Option::CHECK, "nf", "frames", "999",
			"number of frames in the movie."}),
		Option({Option::VALUED | Option::CHECK, "vf", "vel-factor", "0.01",
			"velocity vector factor."}),
		Option({Option::VALUED | Option::CHECK, "vw", "vel-width", "1",
			"velocity vector width."}),
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

	TeX label_height; label_height << "lg";
	auto svg_height_ = label_height.svg();
	double sh = svg_height_->height();

	double rad = 30;
	auto cluster_label = make_cluster_label_material(rv);
	auto wall_material = make_wall_material(rv);
	auto filament_material = make_filament_material(rv);
	auto galaxy_material = make_galaxy_material(rv);
	auto velocity_material = scale_material(make_vel_material(rv, argv.get<double>("vel-width")), rad/(sqrt(2)));

	if (argv.get<bool>("slices"))
	{
		Norma_filter cf(140, 80);
		auto pp_walls = cf(polygons);
		auto pp_fills = cf(segments);
		
			Shell_filter S(20, 80);
			Array<ptr<RenderObject>> scene;

			scene.push_back(ptr<RenderObject>(new PolygonObject(
				S(pp_walls), make_rainbow_material(rv, 20, 50))));
			scene.push_back(ptr<RenderObject>(new SegmentObject(
				S(pp_fills), make_filament_material(false))));
			if (galaxies)
			scene.push_back(ptr<RenderObject>(new VertexObject(
				S(cf(*galaxies)), galaxy_material)));
			if (clusters)
			scene.push_back(ptr<RenderObject>(new VertexObject(
				S(cf(*clusters)), cluster_label)));
			if (velocities)
			scene.push_back(ptr<RenderObject>(new VectorObject(
				S(cf(*velocities)), velocity_material, argv.get<double>("vel-factor"))));

			Spherical_rotation side = eq_to_sg; //(0,-M_PI/2-0.0001,M_PI);
			Spherical_rotation sph_id(0, M_PI/2-1e-5, 0);

			add_coordinates(scene, side, 1e6);
			scene.push_back(ptr<RenderObject>(new SegmentObject(
				make_parallel(side*ga_to_eq, 5), [] (Info info, Context cx)
			{
						cx->set_source_rgb(0.5, 0.5, 0.5);
						cx->set_line_width(0.002);
						cx->stroke();
			})));
			scene.push_back(ptr<RenderObject>(new SegmentObject(
				make_parallel(side*ga_to_eq, -5), [] (Info info, Context cx)
			{
						cx->set_source_rgb(0.5, 0.5, 0.5);
						cx->set_line_width(0.002);
						cx->stroke();
			})));

			Point pointing = cf.center();
			double ra, dec; side(0, M_PI/2, ra, dec);
			Vector shub = -cf.shub();//-spherical_to_cartesian(ra, dec, L);

			side(6./12. * M_PI, 0, ra, dec);
			Vector angle = shub*0.1 + spherical_to_cartesian(ra, dec, L)*0.2;
					
			auto C = make_ptr<Map_projection_camera>(
				Point(L/2,L/2,L/2), pointing, shub,
			//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
				Map_projection(Aitoff_Hammer));
			/*
			auto C = make_ptr<Camera>(
				Point(90,90,90) + angle, 
				pointing, shub,
				//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
					
					scaled_parallel_projection(0.015));*/

			unsigned rx = argv.get<unsigned>("resx"),
				 ry = argv.get<unsigned>("resy");
			auto R = Renderer::Image(rx, ry);
			R->apply([rx, ry] (Context cx)
			{
				cx->scale(ry*0.5, ry*0.5);
				cx->translate((1.0 * rx)/ry, 1.0);
			});
			R->render(scene, C);
			R->write_to_png(Misc::format("norma-skun.png"));
			R->finish();

		return;
	}

	double soff = argv.get<double>("offset");
	Norma_filter cf(15, 180, soff);
	Array<ptr<RenderObject>> scene;
	scene.push_back(ptr<RenderObject>(new PolygonObject(
		cf(polygons), make_rainbow_material(rv, 1.3/1.5, 1.8/1.5))));
	scene.push_back(ptr<RenderObject>(new SegmentObject(
		cf(segments), make_filament_material(false))));
	if (galaxies)
	scene.push_back(ptr<RenderObject>(new VertexObject(
		cf(*galaxies), galaxy_material)));
	if (clusters)
	scene.push_back(ptr<RenderObject>(new VertexObject(
		cf(*clusters), cluster_label)));
	if (velocities)
	scene.push_back(ptr<RenderObject>(new VectorObject(
		cf(*velocities), velocity_material, argv.get<double>("vel-factor"))));

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
			cx->scale(ry*0.5, ry*0.5);
			cx->translate(1.0, 1.0);
		});
		R->render(scene, C);
		R->write_to_png("norma-sky.png");
		R->finish();
		return;
	}

	if (argv.get<bool>("side"))
	{
		// z
		Point pointing = Point(90,90,90);
		{		
		auto C = make_ptr<Camera>(
			pointing + cf.normal()*50, pointing, -(cf.center() - Point(L/2, L/2, L/2)),
			//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
				scaled_parallel_projection(0.02));

		unsigned rx = argv.get<unsigned>("resx"),
			 ry = argv.get<unsigned>("resy");
		auto R = Renderer::Image(rx, ry);
		R->apply([rx, ry] (Context cx)
		{
			cx->scale(ry*0.25, ry*0.25);
			cx->translate(2.0, 2.0);
			cx->set_line_join(Cairo::LINE_JOIN_ROUND);
		});
		R->render(scene, C);
		R->write_to_png(Misc::format("norma-z-", soff, ".png"));
		R->finish();
		}

	/*	//y
		{		
		auto C = make_ptr<Camera>(
			pointing + Vector(0.001, L, 0.001), pointing, -(pointing - Point(L/2, L/2, L/2)),
			//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
				scaled_parallel_projection(0.03));

		unsigned rx = argv.get<unsigned>("resx"),
			 ry = argv.get<unsigned>("resy");
		auto R = Renderer::Image(rx, ry);
		R->apply([rx, ry] (Context cx)
		{
			cx->scale(ry*0.5, ry*0.5);
			cx->translate(1.0, 1.0);
			cx->set_line_join(Cairo::LINE_JOIN_ROUND);
		});
		R->render(scene, C);
		R->write_to_png("pp-y.png");
		R->finish();
		}

		//xyz
		{		
		auto C = make_ptr<Camera>(
			pointing + Vector(L, L, L), pointing, -(pointing - Point(L/2, L/2, L/2)),
			//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
				scaled_parallel_projection(0.03));

		unsigned rx = argv.get<unsigned>("resx"),
			 ry = argv.get<unsigned>("resy");
		auto R = Renderer::Image(rx, ry);
		R->apply([rx, ry] (Context cx)
		{
			cx->scale(ry*0.5, ry*0.5);
			cx->translate(1.0, 1.0);
			cx->set_line_join(Cairo::LINE_JOIN_ROUND);
		});
		R->render(scene, C);
		R->write_to_png("pp-xyz.png");
		R->finish();
		}*/
		return;
	}

	Array<double> th;
	double N_frames = argv.get<double>("frames");
	for (double t = 0; t <= 4*M_PI; t += 4*M_PI/N_frames) th.push_back(t);

	unsigned rx = argv.get<unsigned>("resx"),
		 ry = argv.get<unsigned>("resy");
	for (unsigned i = 0; i < th.size(); ++i)
	{	
		Point p = cf.center();
		Vector dp(cos(th[i]), sin(th[i]), sin(3./2*th[i]));
		auto C = make_ptr<Camera>(p + dp * 50, cf.center(), Vector(0,0,1),
			scaled_parallel_projection(0.015));
		auto R = Renderer::Image(rx, ry);
		R->apply(prepare_context_extract_norma(rx, ry, rad*2, rv));
		R->render(scene, C);
		R->write_to_png(Misc::format("perseus-b-", std::setfill('0'), std::setw(3), i, ".png"));
		R->finish();
	}
}


#include "base/global.hh"
System::Global<Command> _COMMAND_NORMA("norma", command_norma);

