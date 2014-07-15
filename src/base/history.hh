#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <map>

#include "argv.hh"

namespace System
{
	class History: public std::multimap<time_t, std::string>
	{
		public:
			History() {}
			History(std::istream &in);
			void update(std::string const &command);
			void from_file(std::istream &in);
			void to_file(std::ostream &out) const;
			void to_txt_file(std::ostream &out) const;

			History &operator<<(System::Argv const &C);
	};
}

