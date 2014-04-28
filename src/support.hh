#pragma once
#include <string>

extern std::string date_string();
extern std::string seconds_since_epoch_string();
extern std::string time_string(double t);
extern std::string timed_filename(std::string const &id, std::string const &stage, float b, std::string const &ext);

