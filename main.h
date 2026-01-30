#ifndef main_H
#define main_H


#include "custom_math.h"
using custom_math::vector_3;
using custom_math::vector_4;

using custom_math::line_segment_3;

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;
using std::to_string;

#include <sstream>
using std::ostringstream;
using std::istringstream;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <set>
using std::set;

#include <map>
using std::map;

#include <utility>
using std::pair;

#include <mutex>
using std::mutex;

#include <thread>
using std::thread;

#include <random>
std::mt19937 generator(0);
std::uniform_real_distribution<real_type> dis(0.0, 1.0);

#include <optional>
#include <utility>
using namespace std;

const real_type pi = 4.0 * atan(1.0);



#endif
