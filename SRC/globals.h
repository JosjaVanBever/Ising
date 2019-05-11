
#ifndef GLOBALS
#define GLOBALS

#include <iostream>
#include <bitset>
#include <string>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "Htfim.h"
using namespace std;

const string syntax =
	"Syntaxis: ising.exe -L <sites> [-J <J> -g <g>]";
typedef boost::dynamic_bitset<> Bitset;
typedef std::vector<double> State;

#endif