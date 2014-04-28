#include "support.hh"
#include <string>
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cmath>

std::string date_string()
{
	time_t tsec = time(NULL);
	struct tm *T = gmtime(&tsec);
	std::ostringstream s;
	s << std::setfill('0') 
		<< std::setw(2) << T->tm_year - 100
		<< std::setw(2) << T->tm_mon + 1 
		<< std::setw(2) << T->tm_mday; 
	return s.str();
}

std::string seconds_since_epoch_string()
{
	std::ostringstream ss;
	ss << time(NULL);
	return ss.str();
}

std::string time_string(double t)
{
	std::ostringstream s;
	s << std::setfill('0') << std::setw(5) << static_cast<int>(round(t * 10000));
	return s.str();
}


std::string timed_filename(std::string const &id, std::string const &stage, float b, std::string const &ext)
{
	std::ostringstream s;

	if (b < 0.0)
	{
		s 	<< id << "." << stage << "." << "init" 
			<< "." << ext;
	}
	else
	{
		s 	<< id << "." << stage << "." << std::setfill('0') << std::setw(5)
			<< static_cast<int>(round(b * 10000)) 
			<< "." << ext;
	}
	return s.str();
}

