#ifndef WRITEVTK_HPP
#define WRITEVTK_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "mesh2d.hpp"

namespace challenge3 {

    void generateVTKFile(const std::string & filename,const challenge3::Mesh2D & mesh, const std::vector<std::vector<double>> & solution, 
                     size_t nx, size_t ny, double hx, double hy);
}

#endif 