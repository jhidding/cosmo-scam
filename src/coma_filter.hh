#pragma once
#include "base/common.hh"
#include "geometry/geometry.hh"
#include "two_mass.hh"

namespace Scam
{
	using TwoMass::radians;
	using TwoMass::ga_to_sg;
	using TwoMass::spherical_to_cartesian;

	class Tully_filter
	{
  		// A262	136.59	-25.09	 0.0161	
  		// A426	150.39	-13.38	 0.0183	Perseus
  		// A347	141.17	-17.63	 0.0187	

		Sphere 	S;
		Vector  hub;
		//Array<Segment> fogs;

		public:                 
			Tully_filter(double radius)
			{
				S = Sphere(Point(75, 90, 125), radius);
				hub = Vector(0, 0, 1);
			}

			Point center() const { return S.origin(); }
			Vector shub() const { return hub; } 
			Array<Segment> fog() const { return Array<Segment>(0); }

			Array<Vertex> operator()(Array<Vertex> A) const
			{
				Array<Vertex> B;
				for (Vertex const &a : A)
				{
					if (S.is_below(a)) // and Top.is_below(a) and Bottom.is_below(a))
						B.push_back(a);
				}
				return B;
			}

			Array<Polygon> operator()(Array<Polygon> A) const
			{
				Array<Polygon> B;

				for (Polygon const &a : A)
				{
					auto T = S.split_polygon(a);
						if (not T.first) continue;
					/*auto U = Top.split_polygon(*T.first);
						if (not U.first) continue;
					auto V = Bottom.split_polygon(*U.first);
						if (not V.first) continue;*/
					B.push_back(*T.first);
				}

				return B;
			}

			Array<Segment> operator()(Array<Segment> A) const
			{
				Array<Segment> B;

				for (Segment const &a : A)
				{
					auto T = S.split_segment(a);
						if (not T.first) continue;
					/*auto U = Top.split_segment(*T.first);
						if (not U.first) continue;
					auto V = Bottom.split_segment(*U.first);
						if (not V.first) continue;*/
					B.push_back(*T.first);
				}

				return B;
			}
	};

	class Perseus_filter
	{
  		// A262	136.59	-25.09	 0.0161	
  		// A426	150.39	-13.38	 0.0183	Perseus
  		// A347	141.17	-17.63	 0.0187	

		Sphere 	S;
		Vector  hub;
		Array<Segment> fogs;

		public:                 
			Perseus_filter(double radius)
			{
				Point Earth(90,90,90);

				double ra, dec;
				ga_to_sg(radians(136.59), radians(-25.09), ra, dec);
				Point A262 = Earth + spherical_to_cartesian(ra, dec, 0.0161 * 2997.92458);
				ga_to_sg(radians(150.39), radians(-13.38), ra, dec);
				Point A426 = Earth + spherical_to_cartesian(ra, dec, 0.0183 * 2997.92458);
				ga_to_sg(radians(141.17), radians(-17.63), ra, dec);
				Point A347 = Earth + spherical_to_cartesian(ra, dec, 0.0187 * 2997.92458);

				S = Sphere(A262, radius);
				hub = A262 - A426;

				Vector  v = (A262 - Earth).normalize(),
					w = (A426 - Earth).normalize(),
					x = (A347 - Earth).normalize();

				for (int i = -20; i <= 20; ++i) {
					fogs.push_back(Segment(A262 + v * (i * radius/10.), A262 + v * ((i+1) * radius/10.)));
					fogs.push_back(Segment(A426 + w * (i * radius/10.), A426 + w * ((i+1) * radius/10.)));
					fogs.push_back(Segment(A347 + x * (i * radius/10.), A347 + x * ((i+1) * radius/10.)));
				}
			}

			Point center() const { return S.origin(); }
			Vector shub() const { return hub; } 
			Array<Segment> fog() const { return fogs; }

			Array<Vertex> operator()(Array<Vertex> A) const
			{
				Array<Vertex> B;
				for (Vertex const &a : A)
				{
					if (S.is_below(a)) // and Top.is_below(a) and Bottom.is_below(a))
						B.push_back(a);
				}
				return B;
			}

			Array<Polygon> operator()(Array<Polygon> A) const
			{
				Array<Polygon> B;

				for (Polygon const &a : A)
				{
					auto T = S.split_polygon(a);
						if (not T.first) continue;
					/*auto U = Top.split_polygon(*T.first);
						if (not U.first) continue;
					auto V = Bottom.split_polygon(*U.first);
						if (not V.first) continue;*/
					B.push_back(*T.first);
				}

				return B;
			}

			Array<Segment> operator()(Array<Segment> A) const
			{
				Array<Segment> B;

				for (Segment const &a : A)
				{
					auto T = S.split_segment(a);
						if (not T.first) continue;
					/*auto U = Top.split_segment(*T.first);
						if (not U.first) continue;
					auto V = Bottom.split_segment(*U.first);
						if (not V.first) continue;*/
					B.push_back(*T.first);
				}

				return B;
			}
	};

	class Coma_filter
	{
		// the plane of the Coma great wall aligns with these
		// three Abell clusters, A779, Leo and Coma.
		//  A779	191.07	 44.41	0.0226 
		// A1367	234.81	 73.03	0.0215 Leo
		// A1656	 58.09	 87.96	0.0232 Coma

		Plane 	P, Top, Bottom;
		Sphere 	S;
		Vector  hub;
		Array<Segment> fogs;

		public:                 
			Coma_filter(double width, double radius)
			{
				double ra, dec;
				//ga_to_sg(radians(191.07), radians(44.41), ra, dec);
				//Point A779  = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0266 * 2997.92458);
				Point Q(45, 160, 155);
				ga_to_sg(radians(234.81), radians(73.03), ra, dec);
				Point A1367 = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0215 * 2997.92458);
				ga_to_sg( radians(58.09), radians(87.96), ra, dec);
				Point A1656 = Point(90,90,90) + spherical_to_cartesian(ra, dec, 0.0232 * 2997.92458);

				P = Plane(A1656, Vector::cross(Q - A1656, A1367 - A1656).normalize());
				Top = Plane(P.origin() + P.normal() * (width/2), P.normal());
				Bottom = Plane(P.origin() - P.normal() * (width/2), -P.normal());
				S = Sphere(P.origin(), radius);

				hub = A1367 - A1656;
				Point Earth(90,90,90);

				Vector  v = (A1656 - Earth).normalize(),
					w = (A1367 - Earth).normalize();

				for (int i = -20; i <= 20; ++i) {
					fogs.push_back(Segment(A1656 + v * (i * radius/10.), A1656 + v * ((i+1) * radius/10.)));
					fogs.push_back(Segment(A1367 + w * (i * radius/10.), A1367 + w * ((i+1) * radius/10.)));
				}
			}

			Point center() const { return P.origin(); }
			Vector shub() const { return hub; } 
			Array<Segment> fog() const { return fogs; }

			Array<Vertex> operator()(Array<Vertex> A) const
			{
				Array<Vertex> B;
				for (Vertex const &a : A)
				{
					if (S.is_below(a) and Top.is_below(a) and Bottom.is_below(a))
						B.push_back(a);
				}
				return B;
			}

			Array<Polygon> operator()(Array<Polygon> A) const
			{
				Array<Polygon> B;

				for (Polygon const &a : A)
				{
					auto T = S.split_polygon(a);
						if (not T.first) continue;
					auto U = Top.split_polygon(*T.first);
						if (not U.first) continue;
					auto V = Bottom.split_polygon(*U.first);
						if (not V.first) continue;
					B.push_back(*V.first);
				}

				return B;
			}

			Array<Segment> operator()(Array<Segment> A) const
			{
				Array<Segment> B;

				for (Segment const &a : A)
				{
					auto T = S.split_segment(a);
						if (not T.first) continue;
					auto U = Top.split_segment(*T.first);
						if (not U.first) continue;
					auto V = Bottom.split_segment(*U.first);
						if (not V.first) continue;
					B.push_back(*V.first);
				}

				return B;
			}
	};
}

