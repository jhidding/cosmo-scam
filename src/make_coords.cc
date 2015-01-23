#include "make_coords.hh"
using namespace Scam;
using namespace TwoMass;

Array<Segment> Scam::make_meridian(Spherical_rotation const &sph, double longitude, double r)
{
	Array<Segment> A;
	Point O(0,0,0);
	for (int lattitude = -90; lattitude < 90; lattitude += 2)
	{
		double l1, b1; sph(radians(longitude), radians(lattitude), l1, b1);
		double l2, b2; sph(radians(longitude), radians(lattitude+2), l2, b2);
		Point p = O + spherical_to_cartesian(l1, b1, r),
		      q = O + spherical_to_cartesian(l2, b2, r);
		A.push_back(Segment(p, q));
	}
	return A;
}

Array<Segment> Scam::make_parallel(Spherical_rotation const &sph, double lattitude, double r)
{
	Array<Segment> A;
	Point O(0,0,0);
	for (int longitude = -180; longitude < 180; longitude += 2)
	{
		double l1, b1; sph(radians(longitude), radians(lattitude), l1, b1);
		double l2, b2; sph(radians(longitude+2), radians(lattitude), l2, b2);
		Point p = O + spherical_to_cartesian(l1, b1, r),
		      q = O + spherical_to_cartesian(l2, b2, r);
		A.push_back(Segment(p, q));
	}
	return A;
}

void Scam::add_coordinates(Array<ptr<RenderObject>> scene, Spherical_rotation const &sph)
{
	for (int b = -75; b <= 75; b += 15)
	{
		if (b==0)
		scene.push_back(ptr<RenderObject>(new SegmentObject(
			make_parallel(sph, b), [] (Info info, Context cx)
			{
				cx->set_source_rgb(0.5, 0.5, 0.5);
				cx->set_line_width(0.004);
				cx->stroke();
			})));
		else
		scene.push_back(ptr<RenderObject>(new SegmentObject(
			make_parallel(sph, b), [] (Info info, Context cx)
			{
				cx->set_source_rgb(0.5, 0.5, 0.5);
				cx->set_line_width(0.002);
				cx->stroke();
			})));
	}

	for (int l = -180; l < 180; l += 15)
	{
		scene.push_back(ptr<RenderObject>(new SegmentObject(
			make_meridian(sph, l), [] (Info info, Context cx)
			{
				cx->set_source_rgb(0.5, 0.5, 0.5);
				cx->set_line_width(0.002);
				cx->stroke();
			})));
	}
}

