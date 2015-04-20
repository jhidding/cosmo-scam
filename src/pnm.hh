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

#include "support.hh"
#include "cosmic.hh"
#include "materials.hh"
#include "coma_filter.hh"
#include "make_coords.hh"

namespace Scam {
	extern Material scale_material(Material M, double scale);

	class Shell_filter
	{
		Sphere 	S1, S2;

		public:
			Shell_filter(double r1, double r2):
				S1(Point(90,90,90), r1),
				S2(Point(90,90,90), r2)
			{}


			Array<Vertex> operator()(Array<Vertex> A) const
			{
				Array<Vertex> B;
				for (Vertex const &a : A)
				{
					if (S2.is_below(a) and !S1.is_below(a))
						B.push_back(a);
				}
				return B;
			}

			Array<Polygon> operator()(Array<Polygon> A) const
			{
				Array<Polygon> B;

				for (Polygon const &a : A)
				{
					auto T = S2.split_polygon(a);
						if (not T.first) continue;
					auto U = S1.split_polygon(*T.first);
						if (not U.second) continue;
					B.push_back(*U.second);
				}

				return B;
			}

			Array<Segment> operator()(Array<Segment> A) const
			{
				Array<Segment> B;

				for (Segment const &a : A)
				{
					auto T = S2.split_segment(a);
						if (not T.first) continue;
					auto U = S1.split_segment(*T.first);
						if (not U.second) continue;
					B.push_back(*U.second);
				}

				return B;
			}
	};

	template <typename T>
	void add_filtered_coordinates(Array<ptr<RenderObject>> scene, Spherical_rotation const &sph, double r, T f)
	{
		for (int b = -75; b <= 75; b += 15)
		{
			if (b==0)
			scene.push_back(ptr<RenderObject>(new SegmentObject(
				f(make_parallel(sph, b, r)), [] (Info info, Context cx)
				{
					cx->set_source_rgb(0.5, 0.5, 0.5);
					cx->set_line_width(0.004);
					cx->stroke();
				})));
			else
			scene.push_back(ptr<RenderObject>(new SegmentObject(
				f(make_parallel(sph, b, r)), [] (Info info, Context cx)
				{
					cx->set_source_rgb(0.5, 0.5, 0.5);
					cx->set_line_width(0.002);
					cx->stroke();
				})));
		}

		for (int l = -180; l < 180; l += 15)
		{
			scene.push_back(ptr<RenderObject>(new SegmentObject(
				f(make_meridian(sph, l, r)), [] (Info info, Context cx)
				{
					cx->set_source_rgb(0.5, 0.5, 0.5);
					cx->set_line_width(0.002);
					cx->stroke();
				})));
		}
	}
}

