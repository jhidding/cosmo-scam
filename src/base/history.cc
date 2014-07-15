#include "history.hh"
#include "common.hh"
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace System;

inline std::string strip_newlines(std::string const &str)
{
	size_t pos = str.find_first_of('\n');
	return str.substr(0, pos);
}

History::History(std::istream &in)
{
	from_file(in);
}

void History::update(std::string const &command)
{
	insert(std::pair<time_t, std::string>(time(NULL), command));
}

void History::to_file(std::ostream &out) const
{
	std::ostringstream s;
	for (const_iterator i = begin(); i != end(); ++i)
	{
		s << i->first << ":" << strip_newlines(i->second) << std::endl;
	}

	std::string output(s.str());
	std::vector<char> raw(output.length());
	std::copy(output.begin(), output.end(), raw.begin());

	write_block(out, raw);
}

void History::to_txt_file(std::ostream &out) const
{
	out << "# history\n";
	for (const_iterator i = begin(); i != end(); ++i)
	{
		out << "# " << i->first << ":" << strip_newlines(i->second) << std::endl;
	}
	out << std::endl;
}

void History::from_file(std::istream &in)
{
	std::vector<char> raw;
	read_block(in, raw);
	std::string input(raw.begin(), raw.end());

	unsigned pos = 0;
	while (input.length() > pos)
	{
		int p1 = input.find(':', pos),
		    p2 = input.find('\n', pos);
		
		std::string t_str = input.substr(pos, p1 - pos),
		            cmd   = input.substr(p1 + 1, p2 - p1 - 1);

		time_t time = from_string<time_t>(t_str);
		insert(value_type(time, cmd));
		pos = p2 + 1;
	}
}

History &History::operator<<(Argv const &C)
{
	update(C.command_line());
	return *this;
}

