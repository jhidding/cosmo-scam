#pragma once

#include <string>

namespace System
{
	extern std::string date_string();
	extern std::string seconds_since_epoch_string();

	extern std::string 
		timed_filename(std::string const &id, std::string const &stage, float b);
}

