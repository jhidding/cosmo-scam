#include "materials.hh"
#include "tex/tex.hh"
#include "material/colour.hh"

using namespace Scam;

Material Scam::make_cluster_label_material(bool rv)
{
	return [rv] (Info info, Context cx)
	{
		auto ox_ = info.get<double>("ox"), oy_ = info.get<double>("oy");
		//auto name_ = info.get_str("name");
		auto svg_fn_ = info.get_str("svg");
		double ox = 0.15; if (ox_) ox = *ox_;
		double oy = 0.15; if (oy_) oy = *oy_;

		if (not svg_fn_) return;

		SvgFile svg(*svg_fn_);
		double s = 0.05 / 14;
		double w = svg.width() * s;
		double x, y;

		cx->rel_line_to(ox, oy);
		cx->get_current_point(x, y);
		cx->set_line_width(0.007);
		cx->set_source_rgb(0,0,0);
		cx->stroke_preserve();
		cx->set_line_width(0.003);
		if (rv)
		{
			cx->set_source_rgb(0,0.5,1);
		}
		else
		{	
			cx->set_source_rgb(1,0.9,0.5);
		}
		cx->stroke();

		cx->move_to(x, y);
		cx->rel_move_to(-w/2 - 0.05/4, (oy < 0 ? -0.05 - 0.05/4 : 0.05/4));
		cx->rel_curve_to(0.01/4, -0.03/4, 0.02/4, -0.04/4, 0.05/4, -0.05/4);
		cx->rel_line_to(w, 0);
		cx->rel_curve_to(0.03/4, 0.01/4, 0.04/4, 0.02/4, 0.05/4, 0.05/4);
		cx->rel_line_to(0, 0.05);
		cx->rel_curve_to(-0.01/4, 0.03/4, -0.02/4, 0.04/4, -0.05/4, 0.05/4);
		cx->rel_line_to(-w, 0);
		cx->rel_curve_to(-0.03/4, -0.01/4, -0.04/4, -0.02/4, -0.05/4, -0.05/4);
		cx->rel_line_to(0, -0.05);
		cx->close_path();

		if (not rv)
		{
			cx->set_source_rgba(1,0.5,0,0.7);
			cx->fill_preserve();
			cx->set_source_rgb(1,0.5,0);
			cx->stroke();
		}
		else
		{
			cx->set_source_rgba(0,0.5,1,0.5);
			cx->fill_preserve();
			cx->set_source_rgb(0,0.5,1);
			cx->stroke();
		}

		cx->save();
		cx->translate(x - w/2, (oy < 0 ? y - 0.05 - 0.05/4 : y + 0.05/4));
		cx->scale(s, s);
		svg.render(cx);
		cx->restore();
	};
}

Material Scam::make_wall_material(bool rv)
{
	return [rv] (Info info, Context cx)
	{
		auto s_ = info.get<double>("incidence");
		auto d_ = info.get<double>("density");
		double s, d;
		if (s_) s = fabs(*s_);
		if (d_) d = sqrt(*d_);

		d = std::max(5., d); d = std::min(10., d);
		d = (d - 5.) / 5.;


		auto Z = (rv ?
			Colour::HSVA(d*0.1667, 0.9 - d*0.35, 1.0-s/3, 0.5 + d/2 - s/2) :
			Colour::HSVA(0.1667 + d/3, 0.5 + s/5 + d*0.3, 1.0 - d*0.5, 0.5 + d/2 - s/2));

		cx->set_source_rgba(Z.r(), Z.g(), Z.b(), Z.a());
		cx->fill_preserve();
		
		if (rv) cx->set_source_rgba(0,0,0,0.2);
		else cx->set_source_rgba(1,1,1,0.18);

		cx->set_line_width(0.001);
		cx->stroke();
	};
}

Material Scam::make_filament_material(bool rv)
{
	return [rv] (Info info, Context cx)
	{
		auto d_ = info.get<double>("density");
		double d;
		if (d_) d = sqrt(*d_);
		d = std::max(10., d); d = std::min(30., d);
		d = (d - 10.) / 20.;
		cx->set_line_cap(Cairo::LINE_CAP_ROUND);

		if (rv)
		{
			cx->set_line_width(0.023 * d);
			cx->set_source_rgba(0,0,0,0.5);
			cx->stroke_preserve();
			cx->set_line_width(0.02 * d);
			cx->set_source_rgba(1,1,1,0.7);
			cx->stroke();
		}
		else 
		{
			cx->set_line_width(0.023 * d);
			cx->set_source_rgba(1,1,1,0.5);
			cx->stroke_preserve();
			cx->set_line_width(0.02 * d);
			cx->set_source_rgba(0,0,0,0.7);
			cx->stroke();
		}
	};
}

Material Scam::make_galaxy_material(bool rv)
{
	return [rv] (Info info, Context cx)
	{
		auto m_ = info.get<double>("magnitude");
		double m = 1.0; if (m_) m = *m_;
		double s = (12.5 - m)/2.5;
		if (rv) cx->set_source_rgba(0.0,0.4,0.8,0.5);
		else cx->set_source_rgba(1,0,0,0.5);
		cx->rel_move_to(-0.01*s,0);
		cx->rel_line_to(0.01*s, -0.01*s);
		cx->rel_line_to(0.01*s, 0.01*s);
		cx->rel_line_to(-0.01*s, 0.01*s);
		cx->rel_line_to(-0.01*s, -0.01*s);
		//cx->set_line_width(0.01);
		//cx->rel_line_to(0,0);
		//cx->stroke();
		cx->close_path();
		cx->fill_preserve();
		cx->set_line_width(0.001);
		if (rv) cx->set_source_rgba(1,1,1,0.8);
		else cx->set_source_rgba(0,0,0,0.8);
		cx->stroke();
	};
}

Material Scam::make_halo_material(bool rv)
{
	return [rv] (Info info, Context cx)
	{
		auto m_ = info.get<double>("mass");
		double m = 1.0; if (m_) m = *m_;
		double s = log(m)/20.;
		if (rv) cx->set_source_rgba(0.0,0.4,0.8,0.5);
		else cx->set_source_rgba(1,0,0,0.5);
		cx->rel_move_to(-0.01*s,0);
		cx->rel_line_to(0.01*s, -0.01*s);
		cx->rel_line_to(0.01*s, 0.01*s);
		cx->rel_line_to(-0.01*s, 0.01*s);
		cx->rel_line_to(-0.01*s, -0.01*s);
		//cx->set_line_width(0.01);
		//cx->rel_line_to(0,0);
		//cx->stroke();
		cx->close_path();
		cx->fill_preserve();
		cx->set_line_width(0.001);
		if (rv) cx->set_source_rgba(1,1,1,0.8);
		else cx->set_source_rgba(0,0,0,0.8);
		cx->stroke();
	};
}

Material Scam::make_vel_material(bool rv, double w)
{
	return [rv, w] (Info info, Context cx)
	{
		auto m_ = info.get<double>("mass");
		auto vx_ = info.get<double>("ux");
		auto vy_ = info.get<double>("uy");
		auto vz_ = info.get<double>("uz");

		double m = 1.0; if (m_) m = *m_;
		double s = sqrt(m) * w;

		double vx = *vx_, vy = *vy_, vz = *vz_;

		// we have a circle of radius r, and put a cone on it
		// out to radius R > r, what is the angle upto which to draw
		// the arc of the circle? r**2 + (tan(th)*r)**2 = R**2
		// 1 + tan(th) = (R/r)**2 -> th = atan((R/r)**2 - 1)

		double r = 0.01*s;
		double R = sqrt(vx*vx + vy*vy);
		double phi = atan2(vy, vx);
		double theta = (R > r ? atan((R*R)/(r*r) - 1) : 0);

		double x, y;
		cx->get_current_point(x, y);
		cx->begin_new_path();
		cx->arc(x, y, r, phi + theta, phi - theta);
		cx->line_to(x + vx, y + vy);
		cx->close_path();

		if (rv) cx->set_source_rgba(0.6,0.6,0.95,0.2);
		else cx->set_source_rgba(0,0,0,0.5);
		cx->fill_preserve();

		if (rv) cx->set_source_rgba(0,0,0,0.2);
		else cx->set_source_rgba(1,1,1,0.5);
		cx->set_line_width(0.001);
		cx->stroke();
	};
}

Material Scam::make_rainbow_material(bool rv, double a, double b)
{
	return [rv, b, a] (Info info, Context cx)
	{
		auto s_ = info.get<double>("incidence");
		auto d_ = info.get<double>("density");
		auto z_ = info.get<double>("z");
		double s, d, z;
		if (s_) s = fabs(*s_);
		if (d_) d = sqrt(*d_);
		if (z_) z = *z_;

		d = std::max(5., d); d = std::min(10., d);
		d = (d - 5.) / 5.;

		z = std::max(a, z); z = std::min(b, z);
		z = (z-a) / (b-a);

		auto Z = (rv ?
			Colour::HSVA(z*2./3, 0.9 - d*0.35, 1.0-s/3, 0.5 + d/2 - s/2) :
			//Colour::HSVA(z*2./3, s, 1.0 - d*0.5, 0.5 + d/2 - s/2));
			Colour::HSVA(z*2./3, 0.5 + s/5 + d*0.3, 1.0 - d*0.5, 0.5 + d/2 - s/2));

		cx->set_source_rgba(Z.r(), Z.g(), Z.b(), Z.a());
		cx->fill_preserve();
		
		if (rv) cx->set_source_rgba(0,0,0,0.2);
		else cx->set_source_rgba(1,1,1,0.18);

		cx->set_line_width(0.001);
		cx->stroke();
	};
}

Material Scam::make_mistwall_material(bool rv, double a, double b)
{
	return [rv, b, a] (Info info, Context cx)
	{
		auto s_ = info.get<double>("incidence");
		auto d_ = info.get<double>("density");
		auto z_ = info.get<double>("z");
		double s, d, z;
		if (s_) s = fabs(*s_);
		if (d_) d = sqrt(*d_);
		if (z_) z = *z_;

		d = std::max(5., d); d = std::min(10., d);
		d = (d - 5.) / 5.;

		z = std::max(a, z); z = std::min(b, z);
		z = (z-a) / (b-a);

		auto Z = (rv ?
			Colour::HSVA(z*2./3, 0.9 - d*0.35, 1.0-s/3, 0.5 + d/2 - s/2) :
			//Colour::HSVA(z*2./3, s, 1.0 - d*0.5, 0.5 + d/2 - s/2));
			Colour::HSVA(z*2./3, 0.5 + s/5 + d*0.3, 1.0 - d*0.5, 0.5 + d/2 - s/2));

		cx->set_source_rgba(Z.r(), Z.g(), Z.b(), Z.a());
		cx->fill_preserve();
		
		if (rv) cx->set_source_rgba(0,0,0,0.2);
		else cx->set_source_rgba(1,1,1,0.18);

		cx->set_line_width(0.001);
		cx->stroke();
	};
}

Material Scam::make_mistfila_material(bool rv, double a, double b)
{
	return [rv] (Info info, Context cx)
	{
		auto d_ = info.get<double>("density");
		auto z_ = info.get<double>("z");
		double d;
		if (d_) d = sqrt(*d_);
		d = std::max(10., d); d = std::min(30., d);
		d = (d - 10.) / 20.;
		cx->set_line_cap(Cairo::LINE_CAP_ROUND);

		if (rv)
		{
			cx->set_line_width(0.023 * d);
			cx->set_source_rgba(0,0,0,0.5);
			cx->stroke_preserve();
			cx->set_line_width(0.02 * d);
			cx->set_source_rgba(1,1,1,0.7);
			cx->stroke();
		}
		else 
		{
			cx->set_line_width(0.023 * d);
			cx->set_source_rgba(1,1,1,0.5);
			cx->stroke_preserve();
			cx->set_line_width(0.02 * d);
			cx->set_source_rgba(0,0,0,0.7);
			cx->stroke();
		}
	};
}

