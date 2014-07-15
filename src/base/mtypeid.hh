#pragma once
#include <typeinfo>
#include <utility>
#include <vector>
#include <map>
#include <string>

#include "common.hh"

namespace System
{
	class TypeRegister
	{
		typedef std::pair<std::string, std::string> pair;
		typedef std::map<std::string, std::string> map;
		static std::vector<pair> prefab;

		static ptr<map> _id_map;
		static map &id_map() { 
			if (not _id_map) _id_map = ptr<map>(
				new map(prefab.begin(), prefab.end()));
			return *_id_map; }

		public:
			template <typename T>
			static std::string const &name()
			{ return id_map()[typeid(T).name()]; }

			template <typename T>
			static void set_name(std::string const &name)
			{ id_map()[typeid(T).name()] = name; }
	};
}

