/*
 * date() creates a date time stamp, for use in filenames
 */

#include "date.hh"

#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;
using namespace System;

string System::date_string()
{
	time_t tsec = time(NULL);
	struct tm *T = gmtime(&tsec);
	ostringstream s;
	s << setfill('0') 
		<< setw(2) << T->tm_year - 100
		<< setw(2) << T->tm_mon + 1 
		<< setw(2) << T->tm_mday; 
	return s.str();
}

string System::seconds_since_epoch_string()
{
	ostringstream ss;
	ss << time(NULL);
	return ss.str();
}

std::string System::timed_filename(std::string const &id, std::string const &stage, float b)
{
	ostringstream s;

	if (b < 0.0)
	{
		s 	<< id << "." << stage << "." << "init" 
			<< ".conan";
	}
	else
	{
		s 	<< id << "." << stage << "." << setfill('0') << setw(5)
			<< static_cast<int>(round(b * 10000)) 
			<< ".conan";
	}
	return s.str();
}

