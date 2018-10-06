#pragma once

#include <string>
#include "ply/ply.hh"
#include "geometry/geometry.hh"
#include "render/render.hh"

#include "base/common.hh"
#include "base/header.hh"
#include "base/history.hh"
#include "base/mvector.hh"

namespace Scam
{
	extern Array<Vertex>  read_vertices(std::istream &fi, ptr<PLY::PLY> ply, PLY::Format format);
	extern Array<Polygon> read_polygons(std::istream &fi, ptr<PLY::PLY> ply, Array<Vertex> vertices, PLY::Format format, double lim = 0.0);
	extern Array<Segment> read_segments(std::istream &fi, ptr<PLY::PLY> ply, Array<Vertex> vertices, PLY::Format format, double lim = 0.0);
	extern Array<Vertex>  read_abell(std::string const &fn, bool dist=false, bool rv=true);
	extern Array<Vertex>  read_haloes(std::string const &fn);
	extern Array<Vertex>  read_fornax(std::string const &fn);
	extern std::function<void (Context)> prepare_context(unsigned w, unsigned h, double r);
	extern std::function<void (Context)> prepare_context_rv(unsigned w, unsigned h, double r);
	extern Array<Vertex> read_velocities(std::string const &fn_i);
}

