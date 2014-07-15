#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include "header.hh"
#include "common.hh"

using namespace System;

Header::Header(std::istream &f)
{
	from_file(f);
}

inline std::string strip_newlines(std::string const &str)
{
	size_t pos = str.find_first_of('\n');
	return str.substr(0, pos);
}

Header &Header::operator<<(Argv const &C)
{
	for (auto &c : C)
		(*this)[c.long_name] = c.value;

	return *this;
}

void Header::to_file(std::ostream &f) const
{
	std::ostringstream s;
	for (const_iterator i = begin(); i != end(); ++i)
	{
		s << i->first << "=" << strip_newlines(i->second) << std::endl;
	}

	std::string output(s.str());
	std::vector<char> raw(output.length());
	std::copy(output.begin(), output.end(), raw.begin());
	write_block(f, raw);
}

void Header::to_txt_file(std::ostream &f) const
{
	f << "# header\n";
	for (const_iterator i = begin(); i != end(); ++i)
	{
		f << "# " << i->first << "=" << strip_newlines(i->second) << std::endl;
	}
	f << std::endl;
}

void Header::from_file(std::istream &f)
{
	clear();

	std::vector<char> raw;
	read_block(f, raw);
	std::string input(raw.begin(), raw.end());

	unsigned pos = 0;
	while (input.length() > pos)
	{
		int p1 = input.find('=', pos),
		    p2 = input.find('\n', pos);
		
		std::string name  = input.substr(pos, p1 - pos),
		            value = input.substr(p1 + 1, p2 - p1 - 1);

		insert(value_type(name, value));
		pos = p2 + 1;
	}
}

