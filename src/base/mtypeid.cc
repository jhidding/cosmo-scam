#include "mtypeid.hh"
#include "mvector.hh"
#include <complex>

using namespace System;
using namespace Scam;

ptr<std::map<std::string, std::string>> TypeRegister::_id_map;

template <typename T>
std::pair<std::string, std::string> makeTypeId(std::string const &name)
{
	return std::pair<std::string, std::string>(typeid(T).name(), name);
}

std::vector<std::pair<std::string, std::string>> 
	TypeRegister::prefab = 
{
	makeTypeId<int>("int32"),
	makeTypeId<unsigned>("uint32"),
	makeTypeId<float>("float32"),
	makeTypeId<double>("float64"),
	makeTypeId<std::complex<float>>("complex64"),
	makeTypeId<std::complex<double>>("complex128"),
	makeTypeId<mVector<int,3>>("3i4"),
	makeTypeId<mVector<int,2>>("2i4"),
	makeTypeId<mVector<float,3>>("3f4"),
	makeTypeId<mVector<float,2>>("2f4"),
	makeTypeId<mVector<double,3>>("3f8"),
	makeTypeId<mVector<double,2>>("2f8"),
};

