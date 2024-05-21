#ifndef WRITEVTK_HPP
#define WRITEVTK_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "mesh2d.hpp"

namespace challenge3 {

    /// @brief write vtk file format of the solution as scalar field
    /// @param filename name of the file produced
    /// @param x0 x origin
    /// @param y0 y origin
    /// @param solution values assumed on the mesh generated by nx,ny point
    /// @param nx number of points along x
    /// @param ny number of points along y
    /// @param hx spacing along x
    /// @param hy spacing along y
    void generateVTKFile(const std::string & filename,double x0, double y0, const std::vector<std::vector<double>> & solution, 
                     size_t nx, size_t ny, double hx, double hy);
}

#endif 