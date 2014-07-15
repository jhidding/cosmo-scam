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
#include "cosmic.hh"

using namespace System;
using namespace Scam;
using namespace TwoMass;

Scam::Array<Scam::Vertex> Scam::read_vertices(std::istream &fi, Scam::ptr<PLY::PLY> ply, PLY::Format format)
{
	Scam::Array<Scam::Vertex> result;
	auto block = PLY::read_element(fi, (*ply)["vertex"], format);
	
	for (auto const &item : block)
	{
		double x = item.get<double>("x"),
		       y = item.get<double>("y"),
		       z = item.get<double>("z");

		result.push_back(Scam::Vertex(x, y, z));
	}

	return result;
}

Scam::Array<Scam::Polygon> Scam::read_polygons(std::istream &fi, Scam::ptr<PLY::PLY> ply, Scam::Array<Scam::Vertex> vertices, PLY::Format format, double lim)
{
	Scam::Array<Scam::Polygon> result;
	auto block = PLY::read_element(fi, (*ply)["face"], format);
	size_t i = 0;
	for (auto const &item : block)
	{
		std::vector<unsigned> indices = item.get_vector<unsigned>("vertex_index");
		double density = item.get<double>("density");
		if (density < lim) continue;

		auto V = Scam::map([vertices] (int i) { return vertices[i]; }, indices);
		Polygon P(V);
		P.set_info("density", density);
		result.push_back(P);
	}

	return result;
}

Scam::Array<Scam::Segment> Scam::read_segments(std::istream &fi, Scam::ptr<PLY::PLY> ply, Scam::Array<Scam::Vertex> vertices, PLY::Format format, double lim)
{
	Scam::Array<Scam::Segment> result;
	auto block = PLY::read_element(fi, (*ply)["edge"], format);
	size_t i = 0;
	for (auto const &item : block)
	{
		unsigned v1 = item.get<unsigned>("vertex1");
		unsigned v2 = item.get<unsigned>("vertex2");
		double density = item.get<double>("density");
		if (density < lim) continue;
		Segment S(vertices[v1], vertices[v2]);
		S.set_info("density", density);
		result.push_back(S);
	}

	return result;
}

Array<Vertex> Scam::read_abell(std::string const &fn, bool dist, bool rv)
{
	std::cerr << "Reading Abell cluster catalog upto redshift 0.03 ... ";
	Array<Vertex> A;
	std::ifstream fi(fn);
	while (!fi.eof())
	{
		Abell::Cluster C; fi >> C;
		if (C.A_ID == "") continue;

		double sg_ra, sg_dec;
		ga_to_sg(radians(C.l), radians(C.b), sg_ra, sg_dec);
		double d = C.z * 2997.92458;
		Vertex v(Point(90,90,90) + spherical_to_cartesian(sg_ra, sg_dec, (dist ? d : 1.0)));
		v.set_info("r", d);

		std::string name;
		if (C.name != "")
		{
			name = C.name;
			v.set_info("name", C.name);
		}
		else
		{
			name = "A"+C.A_ID;
			v.set_info("name", "A"+C.A_ID);
		}

		v.set_info("ox", C.ox); v.set_info("oy", C.oy);
		TeX label; 
		if (rv) label << "\\color{White}" << name;
		else label << name;

		label.make_svg();
		v.set_info("svg", label.file_name());

		A.push_back(v);
	}
	
	Vertex v(Point(90,90,90));
	v.set_info("r", 0);
	v.set_info("name", "Home");
	v.set_info("ox", 0); v.set_info("oy", -0.1);
	TeX label;	
	if (rv) label << "\\color{White}Home";
	else label << "Home";
	label.make_svg();
	v.set_info("svg", label.file_name());
	A.push_back(v);

	std::cerr << "Ok\n";
	return A;
}

template <unsigned R>
using dVector = Scam::mVector<double, R>;

template <unsigned R>
struct VelocityInfo
{
	dVector<R> x, v;
	double     mass;
	int	   type;
};

Array<Vertex> Scam::read_haloes(std::string const &fn)
{
	Array<Vertex> A;
	std::ifstream fi(fn);
	while (!fi.eof())
	{
		FOF::Halo H; fi >> H;
		
		Vertex v(H.x, H.y, H.z);
		v.set_info("mass", H.mass);
		A.push_back(v);
	}
	return A;
}

template <typename T>
Scam::Array<T> load_from_file(std::istream &fi, std::string const &name)
{
	auto pos = fi.tellg();

	while (fi.good())
	{
		Header S(fi);
		if (S["name"] != name)
		{
			skip_block(fi);
		}
		else
		{
			Array<T> data(fi);
			fi.seekg(pos, std::ios::beg);
			return data;
		}
	}

	throw "Couldn't find record " + name + " in file.";
}

extern Scam::Array<Vertex> read_velocities(std::string const &fn_i);

Array<Vertex> Scam::read_velocities(std::string const &fn_i)
{
	std::ifstream fi(fn_i);
	Header H(fi); History I(fi);

	auto P = load_from_file<VelocityInfo<3>>(fi, "nodes");

	Array<Vertex> V;
	for (VelocityInfo<3> const &p : P)
	{
		Vertex v(p.x[0], p.x[1], p.x[2]);
		v.set_info("vx", p.v[0]);
		v.set_info("vy", p.v[1]);
		v.set_info("vz", p.v[2]);
		v.set_info("mass", p.mass);
		V->push_back(v);
	}

	return V;
}

