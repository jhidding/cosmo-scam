/*
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

#include <ctime>
#include <sstream>
#include <iomanip>

using namespace System;
using namespace Scam;
using namespace TwoMass;

Vertex make_label(double sg_ra, double sg_dec, std::string const &name, double ox, double oy)
{
	Vertex v(Point(90,90,90) + spherical_to_cartesian(sg_ra, sg_dec, 1.0));
	v.set_info("ox", ox);
	v.set_info("oy", oy);
	v.set_info("name", name);
	return v;
}

struct Cluster
{
	double ra, dec, v;
	std::string name;
	double ox, oy;

	bool ok(double v_min, double v_max) const
	{
		return (v + 5.0 > v_min) and (v_max > v - 5.0);
	}

	Vertex vertex() const
	{
		double l, b; eq_to_sg(radians(ra), radians(dec), l, b);
		return make_label(l, b, name, ox, oy);
	}
};

extern void command_cosmic(int argc_, char **argv_);
*/
