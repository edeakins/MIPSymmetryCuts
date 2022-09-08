#ifndef READ_ORBIT_H_
#define READ_ORBIT_H_

#include <fstream>
#include <ostream>
#include <iostream>
#include <sstream>
#include <string>
#include "orbital_partition.hpp"

class orbit_reader{
public:
    orbit_reader(const char* filename, int size);
    orbital_partition& get_orbital_partition();
    orbital_partition op;
};

#endif