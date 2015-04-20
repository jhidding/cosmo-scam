#include "pnm.hh"

using namespace Scam;
using namespace System;
using namespace TwoMass;

Material Scam::scale_material(Material M, double scale)
{
	return [M, scale] (Info info, Context cx)
	{
		cx->save();
		cx->scale(scale, scale);
		M(info, cx);
		cx->restore();
	};
}
/*
class RA_Dec_filter
{
	Plane	A1, A2, B1, B2;

	public:
		RA_Dec_filter(Spherical_rotation const &T, double ra_min, double ra_max, double dec_min, double dec_max)
		{
			auto vec = [T] (double l, double b)
			{
				double ra, dec;
				T(l, b, ra, dec);
				return spherical_to_cartesian(ra, dec, 1.0);
			};

			Vector  z  = vec(0, M_PI/2),
				p1 = vec(ra_min, 0),
				p2 = vec(ra_max, 0),
				q1 = vec((ra_min + ra_max)/2, dec_min),
				q2 = vec((ra_min + ra_max)/2, dec_max),
				x  = Vector::cross(q1, q2);

			Point earth(90,90,90);
			A1 = Plane(earth, Vector::cross(z, p1));
			A2 = Plane(earth, Vector::cross(p2, z));
			B1 = Plane(earth, Vector::cross(q1, x));
			B2 = Plane(earth, Vector::cross(x, q2));
		}

		Array<Vertex> operator()(Array<Vertex> A) const
		{
			Array<Vertex> B;
			for (Vertex const &a : A)
			{
				if (A1.is_below(a) and !A2.is_below(a) and B1.is_below(a) and !B2.is_below(a))
					B.push_back(a);
			}
			return B;
		}

		Array<Polygon> operator()(Array<Polygon> A) const
		{
			Array<Polygon> B;

			for (Polygon const &a : A)
			{
				auto T = A1.split_polygon(a);
					if (not T.first) continue;
				auto U = A2.split_polygon(*T.first);
					if (not U.second) continue;
				auto V = B1.split_polygon(*U.second);
					if (not V.first) continue;
				auto W = B2.split_polygon(*V.first);
					if (not W.second) continue;
				B.push_back(*W.second);
			}

			return B;
		}

		Array<Segment> operator()(Array<Segment> A) const
		{
			Array<Segment> B;

			for (Segment const &a : A)
			{
				auto T = A1.split_segment(a);
					if (not T.first) continue;
				auto U = A2.split_segment(*T.first);
					if (not U.second) continue;
				auto V = B1.split_segment(*U.second);
					if (not V.first) continue;
				auto W = B2.split_segment(*V.first);
					if (not W.second) continue;
				B.push_back(*W.second);
			}

			return B;
		}
};
*/

