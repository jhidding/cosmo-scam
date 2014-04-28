#include "cosmic.hh"
#include "tex/tex.hh"

using namespace Scam;

std::function<void (Context)> Scam::prepare_context(unsigned w, unsigned h, double r)
{
	return [w, h, r] (Context cx)
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
		cx->stroke();
		cx->restore();

		cx->save();
		cx->translate(-L/2, L/4 - 2*margin);
		cx->scale(margin*3, margin*3);
		/*cx->rectangle(0,0, 2,1);
		cx->set_line_width(0.001);
		cx->set_source_rgb(1,0,1);
		cx->stroke();*/

		cx->move_to(0.9,0.24);
		cx->rel_curve_to(0.01, 0.03, 0.02, 0.04, 0.05, 0.05);
		cx->rel_line_to(1.0, 0.0);
		cx->rel_curve_to(0.03, -0.01, 0.04, -0.02, 0.05, -0.05);
		cx->rel_line_to(0.0, 0.12);
		cx->rel_curve_to(-0.01, -0.03, -0.02, -0.04, -0.05, -0.05);
		cx->rel_line_to(-1.0, 0.0);
		cx->rel_curve_to(-0.03, 0.01, -0.04, 0.02, -0.05, 0.05);
		cx->rel_line_to(0.0, -0.12);
		cx->close_path();
		cx->set_source_rgb(0,0,0);
		cx->fill();

		TeX text_ruler; 
		text_ruler << "$\\approx " << int(r) << "\\ {\\rm Mpc\\ h^{-1}}$\n";
		auto svg = text_ruler.svg();
		double u = svg->width(), v = svg->height(), H = 0.15;
		cx->save();
		cx->translate(0.9 + 0.55 - u/2*H/v, 0.24 + 0.03 - H);
		cx->scale(H/v, H/v);
		svg->render(cx);
		cx->restore();
		cx->restore();

		cx->save();
		cx->translate(L/2-margin*6, L/4 - 2*margin);
		cx->scale(margin*3, margin*3);
		/*cx->rectangle(0,0, 2,1);
		cx->set_line_width(0.001);
		cx->set_source_rgb(1,0,1);
		cx->stroke();*/
		cx->restore();

		cx->set_line_join(Cairo::LINE_JOIN_ROUND);
	};
}

std::function<void (Context)> Scam::prepare_context_rv(unsigned w, unsigned h, double r)
{
	return [w, h, r] (Context cx)
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

		cx->save();
		cx->translate(-L/2, L/4 - 2*margin);
		cx->scale(margin*3, margin*3);
		/*cx->rectangle(0,0, 2,1);
		cx->set_line_width(0.001);
		cx->set_source_rgb(1,0,1);
		cx->stroke();*/

		cx->move_to(0.9,0.24);
		cx->rel_curve_to(0.01, 0.03, 0.02, 0.04, 0.05, 0.05);
		cx->rel_line_to(1.0, 0.0);
		cx->rel_curve_to(0.03, -0.01, 0.04, -0.02, 0.05, -0.05);
		cx->rel_line_to(0.0, 0.12);
		cx->rel_curve_to(-0.01, -0.03, -0.02, -0.04, -0.05, -0.05);
		cx->rel_line_to(-1.0, 0.0);
		cx->rel_curve_to(-0.03, 0.01, -0.04, 0.02, -0.05, 0.05);
		cx->rel_line_to(0.0, -0.12);
		cx->close_path();
		cx->set_source_rgb(1,1,1);
		cx->fill();

		TeX text_ruler; 
		text_ruler << "$\\color{White}\\approx " << int(r) << "\\ {\\rm Mpc\\ h^{-1}}$\n";
		auto svg = text_ruler.svg();
		double u = svg->width(), v = svg->height(), H = 0.15;
		cx->save();
		cx->translate(0.9 + 0.55 - u/2*H/v, 0.24 + 0.03 - H);
		cx->scale(H/v, H/v);
		svg->render(cx);
		cx->restore();
		cx->restore();

		/*cx->save();
		cx->translate(L/2-margin*6, L/4 - 2*margin);
		cx->scale(margin*3, margin*3);
		cx->rectangle(0,0, 2,1);
		cx->set_line_width(0.001);
		cx->set_source_rgb(1,0,1);
		cx->stroke();
		cx->restore();*/

		cx->set_line_join(Cairo::LINE_JOIN_ROUND);
	};
}

