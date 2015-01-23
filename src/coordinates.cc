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
#include "make_coords.hh"

using namespace Scam;
using namespace TwoMass;
using namespace System;

std::function<void (Context)> prepare_context(unsigned w, unsigned h)
{
	return [w, h] (Context cx)
	{
		double L = 6;
		double margin = (h * (L/w) - L/2)/2;
		cx->scale(w/L,w/L);
		cx->translate(L/2, L/4 + margin);
		cx->set_source_rgb(0.5,0.5,0.5);
		cx->set_line_width(0.01);

		cx->save();
		cx->scale(2*sqrt(2), sqrt(2));
		cx->arc(0,0,1.0,0,6.283184);
		cx->restore();
		cx->stroke();

		cx->set_line_join(Cairo::LINE_JOIN_ROUND);
	};
}

std::function<void (Context)> prepare_context_rv(unsigned w, unsigned h)
{
	return [w, h] (Context cx)
	{
		double L = 6;
		double margin = (h * (L/w) - L/2)/2;
		cx->scale(w/L,w/L);
		cx->translate(L/2, L/4 + margin);
		cx->set_source_rgb(0.5,0.5,0.5);
		cx->set_line_width(0.01);

		cx->save();
		cx->scale(2*sqrt(2), sqrt(2));
		cx->arc(0,0,1.0,0,6.283184);
		cx->restore();
		cx->stroke();

		cx->set_line_join(Cairo::LINE_JOIN_ROUND);
		cx->set_line_cap(Cairo::LINE_CAP_ROUND);
	};
}

void command_coordinates(int argc_, char **argv_)
{
	System::Argv argv = System::read_arguments(argc_, argv_,
		Option({0, "h", "help", "false", 
			"print this help."}),
		Option({0, "2m", "2mass", "false",
			"include 2Mass into rendering."}),
		Option({0, "r", "reverse", "false",
			"reverse colours."}),
		Option({0, "lv", "local-void", "false",
			"point the camera into the local void"}),
		Option({0, "norma", "norma", "false",
			"point the camera to Norma"}),
		Option({0, "svg", "svg", "false",
			"produce SVG output"}),
		Option({Option::VALUED | Option::CHECK, "rx", "resx", "1920",
			"image size X"}),
		Option({Option::VALUED | Option::CHECK, "ry", "resy", "1080",
			"image size Y"}),
		Option({Option::VALUED | Option::CHECK, "i", "id", "test",
			"identifier for filenames."}),
		Option({Option::VALUED | Option::CHECK, "L", "size", "100",
			"size of the box."}) );

	if (argv.get<bool>("help"))
	{
		std::cerr << "Scam -- 3D vector graphics for science.\n"
			"Copyright Johan Hidding, June 2014 - licence: GPL3.\n\n";
		argv.print(std::cerr);
		exit(0);
	}

	/* How to draw a galaxy */
	auto galaxy_material = (argv.get<bool>("reverse") ?
		[] (Info info, Context cx)
		{
			auto m_ = info.get<double>("magnitude");
			double m = 1.0; if (m_) m = *m_;
			double s = (12.5 - m)/2.5;
			cx->set_source_rgba(0.0,0.4,0.8,0.5);
			cx->rel_move_to(-0.01*s,0);
			cx->rel_line_to(0.01*s, -0.01*s);
			cx->rel_line_to(0.01*s, 0.01*s);
			cx->rel_line_to(-0.01*s, 0.01*s);
			cx->rel_line_to(-0.01*s, -0.01*s);
			cx->close_path();

			cx->fill_preserve();
			cx->set_line_width(0.001);
			cx->set_source_rgba(1,1,1,0.8);
			cx->stroke();
		} :
		[] (Info info, Context cx)
		{
			auto m_ = info.get<double>("magnitude");
			double m = 1.0; if (m_) m = *m_;
			double s = (12.5 - m)/2.5;
			cx->set_source_rgba(0.8,0.4,0.0,0.5);
			cx->rel_move_to(-0.01*s,0);
			cx->rel_line_to(0.01*s, -0.01*s);
			cx->rel_line_to(0.01*s, 0.01*s);
			cx->rel_line_to(-0.01*s, 0.01*s);
			cx->rel_line_to(-0.01*s, -0.01*s);
			cx->close_path();

			cx->fill_preserve();
			cx->set_line_width(0.001);
			cx->set_source_rgba(0,0,0,0.8);
			cx->stroke();
		} );

	auto hour_label = (argv.get<bool>("reverse") ?
		[] (Info info, Context cx)
		{
			Cairo::TextExtents E;
			auto s_ = info.get_str("text");
			cx->select_font_face("Bitstream Charter", Cairo::FONT_SLANT_NORMAL,
				Cairo::FONT_WEIGHT_NORMAL);
			std::string s; if (s_) s = *s_;

			cx->set_font_size(0.05);
			cx->get_text_extents(s, E);

			double x, y; cx->get_current_point(x, y);

			cx->begin_new_path();
			cx->rectangle(x, y, E.width+0.02, E.height+0.02);
			cx->set_source_rgba(0,0,0,0.7);
			cx->fill();

			cx->arc(x, y, 0.010, 0, 2*M_PI);
			cx->set_source_rgb(1,1,1);
			cx->fill();

			cx->move_to(x+0.01,y+0.01+E.height);
			cx->show_text(s);
		} :
		[] (Info info, Context cx)
		{
			Cairo::TextExtents E;
			auto s_ = info.get_str("text");
			cx->select_font_face("Bitstream Charter", Cairo::FONT_SLANT_NORMAL, 
				Cairo::FONT_WEIGHT_NORMAL);
			std::string s; if (s_) s = *s_;

			cx->set_font_size(0.05);
			cx->get_text_extents(s, E);

			double x, y; cx->get_current_point(x, y);

			cx->begin_new_path();
			cx->rectangle(x, y, E.width+0.02, E.height+0.02);
			cx->set_source_rgba(1,1,1,0.7);
			cx->fill();

			cx->arc(x, y, 0.010, 0, 2*M_PI);
			cx->set_source_rgb(0,0,0);
			cx->fill();

			cx->move_to(x+0.01,y+0.01+E.height);
			cx->show_text(s);
		});

	auto circled_label = (argv.get<bool>("reverse") ?
		[] (Info info, Context cx)
		{
			Cairo::TextExtents E;

			auto s_ = info.get_str("text");
			cx->select_font_face("Bitstream Charter", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
			std::string s; if (s_) s = *s_;
			cx->set_source_rgb(1,1,1);
			cx->set_font_size(0.05);
			cx->get_text_extents(s, E);

			double x, y; cx->get_current_point(x, y);
			cx->begin_new_path();
			cx->arc(x, y, E.height*0.8, 0, M_PI*2);
			cx->set_line_width(0.005);
			cx->set_source_rgba(0,0,0,0.7);
			cx->fill_preserve();
			cx->set_source_rgb(1,1,1);
			cx->stroke();

			cx->move_to(x - E.width/2, y + E.height/2);
			cx->show_text(s);
		} :
		[] (Info info, Context cx)
		{
			Cairo::TextExtents E;

			auto s_ = info.get_str("text");
			cx->select_font_face("Bitstream Charter", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
			std::string s; if (s_) s = *s_;
			cx->set_source_rgb(0,0,0);
			cx->set_font_size(0.05);
			cx->get_text_extents(s, E);

			double x, y; cx->get_current_point(x, y);
			cx->begin_new_path();
			cx->arc(x, y, E.height*0.8, 0, M_PI*2);
			cx->set_line_width(0.005);
			cx->set_source_rgba(1,1,1,0.7);
			cx->fill_preserve();
			cx->set_source_rgb(0,0,0);
			cx->stroke();

			cx->move_to(x - E.width/2, y + E.height/2);
			cx->show_text(s);
		});


	auto make_label = [] (Spherical_rotation const &sph, double ra, double dec, std::string const &s)
	{
		double l, b; sph(radians(ra), radians(dec), l, b);
		Vertex v(Point(0,0,0) + spherical_to_cartesian(l, b, 1.0));
		v.set_info("text", s);
		return v;
	};

	Point pointing;
	if (argv.get<bool>("local-void"))
		pointing = Point(-0.0001,0.0001, 1);
	else if (argv.get<bool>("norma"))
		pointing = Point(-1, 0, 0);
	else
		pointing = Point(1, 0, 0);
	Vector shub = (argv.get<bool>("local-void") ?
		Vector(1, 0, 0) :
		Vector(0, 0, 1));
			
	auto C = make_ptr<Map_projection_camera>(
		Point(0,0,0), pointing, shub,
		//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
			Map_projection(Aitoff_Hammer));
//	auto C = make_ptr<Map_projection_camera>(
//		Point(0,0,0), Point(1, 0, 0), Vector(0, 0, 1),
//		//Point(-1.0, 0.5, 1.0), centre, Vector(0, -1, 0),
//			Map_projection(Aitoff_Hammer));

	Array<ptr<RenderObject>> scene;
	Spherical_rotation side = eq_to_sg; //(0,-M_PI/2-0.0001,M_PI);
	Spherical_rotation sph_id(0, M_PI/2-1e-5, 0);

	/* add labels */
	Array<Vertex> hour_labels = {
		make_label(side, 0, 0, "0h"),
		make_label(side, 90, 0, "6h"),
		make_label(side, 180, 0, "12h"),
		make_label(side, 270, 0, "18h") };

	Array<Vertex> circ_labels = {
		make_label(side, 0, 90, "N"),
		make_label(side, 0, -90, "S") };

	scene.push_back(ptr<RenderObject>(new VertexObject(
		hour_labels, hour_label)));	
	scene.push_back(ptr<RenderObject>(new VertexObject(
		circ_labels, circled_label)));	

	/* add galaxies */
	Maybe<Array<Vertex>> galaxies;
	if (argv.get<bool>("2mass"))
		galaxies = Just(read_2mass("2mrs_1175_done.dat", ga_to_sg, false));
	else
		galaxies = Nothing;

	if (galaxies)
	scene.push_back(ptr<RenderObject>(new VertexObject(
		*galaxies, galaxy_material)));

	/* add coordinate lines */
	add_coordinates(scene, side);
	scene.push_back(ptr<RenderObject>(new SegmentObject(
		make_parallel(side*ga_to_eq, 5, 1e5), [] (Info info, Context cx)
	{
				cx->set_source_rgb(0.5, 0.5, 0.5);
				cx->set_line_width(0.008);
				cx->stroke();
	})));
	scene.push_back(ptr<RenderObject>(new SegmentObject(
		make_parallel(side*ga_to_eq, -5, 1e5), [] (Info info, Context cx)
	{
				cx->set_source_rgb(0.5, 0.5, 0.5);
				cx->set_line_width(0.008);
				cx->stroke();
	})));

	/* draw the thing */
	std::string fn_output_png = Misc::format(argv["id"], ".png");
	std::string fn_output_svg = Misc::format(argv["id"], ".svg");
	unsigned rx = argv.get<unsigned>("resx"),
		 ry = argv.get<unsigned>("resy");

	ptr<Renderer> R;
	if (argv.get<bool>("svg"))
	{
		R = Renderer::SVG(rx, ry, fn_output_svg);
	}
	else
	{
		R = Renderer::Image(rx, ry);
	}

	if (argv.get<bool>("reverse"))
		R->apply(prepare_context_rv(rx, ry));
	else
		R->apply(prepare_context(rx, ry));
	R->render(scene, C);
/*	R->apply([] (Context cx)
	{
		cx->set_source_rgb(0.8,0.8,0.8);
		cx->set_line_width(0.01);

		cx->save();
		cx->scale(2*sqrt(2), sqrt(2));
		cx->arc(0,0,1.0,0,6.283184);
		cx->restore();
		cx->stroke();
	});*/
	R->write_to_png(fn_output_png);
	R->finish();
}

#include "base/global.hh"
System::Global<Command> _COMMAND_COORD("coordinates", command_coordinates);

