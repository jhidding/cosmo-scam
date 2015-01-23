#pragma once

#include "base/common.hh"
#include "base/argv.hh"
#include "base/format.hh"
#include "geometry/geometry.hh"
#include "render/render.hh"
#include "render/map_projection.hh"
#include "two_mass.hh"
#include "material/colour.hh"
#include "tex/tex.hh"

namespace Scam {
	using TwoMass::Spherical_rotation;
	Array<Segment> make_meridian(Spherical_rotation const &sph, double longitude, double r = 1e6);
	Array<Segment> make_parallel(Spherical_rotation const &sph, double lattitude, double r = 1e6);
	void add_coordinates(Array<ptr<RenderObject>> scene, Spherical_rotation const &sph);
}

