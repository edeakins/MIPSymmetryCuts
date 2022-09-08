# include "read_orbits.hpp"

orbit_reader::orbit_reader(const char* filename, int size){
    std::ifstream orbit_file(filename);
    std::string file_line;
    int num_orbits = 0;
    op.orbit_start.push_back(0);
    op.orbit.resize(size);
    while (std::getline(orbit_file, file_line)){
        std::istringstream iss(file_line);
        int node;
        while(iss >> node){
            op.element.push_back(node);
            op.orbit.at(node) = num_orbits;
        }
        num_orbits++;
        op.orbit_start.push_back(op.element.size());
    }
}

orbital_partition& orbit_reader::get_orbital_partition(){
    return op;
}