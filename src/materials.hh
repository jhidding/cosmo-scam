#pragma once

#include "render/render.hh"

namespace Scam
{
	extern Material make_cluster_label_material(bool rv);
	extern Material make_wall_material(bool rv);
	extern Material make_filament_material(bool rv);
	extern Material make_galaxy_material(bool rv);
	extern Material make_halo_material(bool rv);
	extern Material make_vel_material(bool rv, double w);
}
