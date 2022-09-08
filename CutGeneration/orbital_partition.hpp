#ifndef ORBITAL_PARTITION_H_
#define ORBITAL_PARTITION_H_

#include <vector>

class orbital_partition{
public:
    std::vector<int> element;
    std::vector<int> orbit;
    std::vector<int> orbit_start;
};

#endif