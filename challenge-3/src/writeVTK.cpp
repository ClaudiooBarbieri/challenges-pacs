#include "writeVTK.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


namespace challenge3 {

    // Generates a STRUCTURES VTK file with a mesh
    void generateVTKFile(const std::string &filename, double x0, double y0, const std::vector<std::vector<double>> &solution,
                         size_t nx, size_t ny, double hx, double hy) {

        // Opens the files
        std::ofstream vtkFile(filename + "_heatmap.vtk"); ///< sol as 2D heatmap
        std::ofstream svtkFile(filename + "_surface.vtk"); ///< sol as 3D surface

        // Check if the files were opened
        if (!vtkFile.is_open() || !svtkFile.is_open()) {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }

        // Write VTK header for heatmap
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Heatmap\n";
        vtkFile << "ASCII\n";                                // File format

        // Write VTK header for surface
        svtkFile << "# vtk DataFile Version 3.0\n";
        svtkFile << "Surface\n";
        svtkFile << "ASCII\n";                                // File format

        // Write grid data for heatmap
        vtkFile << "DATASET STRUCTURED_POINTS\n";                             // Format of the dataset
        vtkFile << "DIMENSIONS " << nx << " " << ny << " 1\n";  // Number of points in each direction
        vtkFile << "ORIGIN " << x0 << " " << y0 << " 0\n";  // Lower-left corner of the structured grid
        vtkFile << "SPACING " << hx << " " << hy << " 1\n";   // Spacing between points in each direction
        vtkFile << "POINT_DATA " << (nx) * (ny) << "\n";                  // Number of points

        // Write scalar field data for heatmap
        vtkFile << "SCALARS solution double\n";               // Description of the scalar field
        vtkFile << "LOOKUP_TABLE default\n";                 // Color table

        // Write grid data for surface
        svtkFile << "DATASET STRUCTURED_GRID\n";
        svtkFile << "DIMENSIONS " << nx << " " << ny << " 1\n";
        svtkFile << "POINTS " << nx * ny << " double\n";

        // Write the points and scalar values
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                vtkFile << solution[i][j] << "\n";
                svtkFile << x0 + i * hx << " " << y0 + j * hy << " " << solution[i][j] << "\n";
            }
        }

        // Include scalar field data for the surface file
        svtkFile << "POINT_DATA " << nx * ny << "\n";
        svtkFile << "SCALARS solution double 1\n";
        svtkFile << "LOOKUP_TABLE default\n";
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << solution[i][j] << "\n";
            }
        }

        std::cout << "VTK heatmap written" << std::endl;
        vtkFile.close();
        std::cout << "VTK surface written" << std::endl;
        svtkFile.close();
    }
}