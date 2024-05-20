#include "writeVTK.hpp"

namespace challenge3{

    // generates a STRUCTURES VTK file with a mesh
    void generateVTKFile(const std::string & filename, const challenge3::Mesh2D & mesh, const std::vector<std::vector<double>> & solution, 
                     size_t nx, size_t ny, double hx, double hy) {

        // opens the file
        std::ofstream vtkFile(filename);

        // check if the file was opened
        if (!vtkFile.is_open()) {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }

        // Write VTK header
        vtkFile <<  "# vtk DataFile Version 3.0\n";
        vtkFile << "Solution on mesh2D\n";
        vtkFile << "ASCII\n";                                // file format


        // Write grid data
        vtkFile << "DATASET STRUCTURED_POINTS\n";                             // format of the dataset
        vtkFile << "DIMENSIONS " << nx << " " << ny << " " << 1 << "\n";  // number of points in each direction
        vtkFile << "ORIGIN " << mesh.getMinX() << " " << mesh.getMinY() << " 0\n";  // lower-left corner of the structured grid
        vtkFile << "SPACING" << " " << hx << " " << hy << " " << 1 << "\n";   // spacing between points in each direction
        vtkFile << "POINT_DATA " << (nx) * (ny) << "\n";                  // number of points

        // Write scalar field data
        vtkFile << "SCALARS solution double\n";               // description of the scalar field
        vtkFile << "LOOKUP_TABLE default\n";                 // color table

        // Write vector field data
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                vtkFile <<  solution[i][j] << "\n";
            }
        }
        std::cout << "VTK written" << std::endl;
    }

}