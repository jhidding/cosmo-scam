#pragma once

#include "render/render.hh"

namespace Scam
{
	extern Material make_cluster_label_material(bool rv);
	extern Material make_wall_material(bool rv, double a = 5.0, double b = 10.0);
	extern Material make_filament_material(bool rv, double a = 10.0, double b = 30.0);
	extern Material make_galaxy_material(bool rv);
	extern Material make_halo_material(bool rv);
	extern Material make_vel_material(bool rv, double w);
	extern Material make_rainbow_material(bool rv, double a=40., double b=70.);
	extern Material make_mistwall_material(bool rv, double a=40., double b=70.);
	extern Material make_mistfila_material(bool rv, double a=40., double b=70.);
}
