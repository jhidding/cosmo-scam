/*
 * Reads one block from a file and interprets it as a header.
 * a header contains <name>=<value> pairs in ASCII separated
 * by semicolons ';'. The entire block is saved as a Fortran
 * array, encapsulated by a 32bit unsigned integer containing
 * the size of the block.
 */

#pragma once

#include <iosfwd>
#include <string>
#include <map>

#include "common.hh"
#include "argv.hh"

namespace System
{
	class Header: public std::map<std::string, std::string>
	{
		public:
			Header() {}
			Header(std::istream &f);
			void to_file(std::ostream &f) const;
			void to_txt_file(std::ostream &f) const;
			void from_file(std::istream &f);
		
			std::string const &operator[](std::string const &key) const
			{ return find(key)->second; }

			std::string &operator[](std::string const &key)
			{ return std::map<std::string, std::string>::operator[](key); }

			template <typename T>
			T get(std::string const &name) const
			{
				auto i = find(name);
				if (i == end())
					throw "variable [" + name + "] not found in config.";

				return from_string<T>(i->second);
			}

			Header &operator<<(Argv const &C);
	};
} // namespace Conan
