#include "mesh2d.hpp"

namespace challenge3 {
    /// factory method to create the mesh given the number of points
    Mesh2D Mesh2D::createWithPoints(const point & bl, const point & tr, size_t nx, size_t ny) {
        return Mesh2D(bl,tr, nx, ny);
    };

    /// factory method to create the mesh given the spacing
    Mesh2D Mesh2D::createWithSpacing(const point & bl, const point & tr, double hx, double hy) {
        return Mesh2D(bl,tr, hx, hy);
    }

    /// constructor given the domain and number of points along each direction
    Mesh2D::Mesh2D(const point & bl, const point & tr, size_t nx_, size_t ny_) : Domain(bl,tr),nx{nx_},ny{ny_} { 
        hx = (maxX-minX)/(nx-1); ///< setting x spacing according to number of points wanted
        hy = (maxY-minY)/(ny-1); ///< setting y spacing according to number of points wanted

        /// initialize points in the grid with thei x and y coordinate
        coordinate.reserve(nx*ny);
        for(double y = minY ; y <= maxY ; y+=hy){
            for(double x = minX ; x <= maxX ; x+=hx){
                //std::cout << "(" << x << "," << y << ")" << std::endl;
                coordinate.push_back({x,y});
                //std::cout << "(" << (*(coordinate.end()-1))[0] << "," << (*(coordinate.end()-1))[1] << ")" << std::endl;
            }
        }
    };

    /// constructor given the domain and the spacing along each direction
    Mesh2D::Mesh2D(const point & bl, const point & tr, double hx_, double hy_) : Mesh2D(bl,tr,static_cast<size_t>(std::ceil((maxX - minX) / hx_)+1),
             static_cast<size_t>(std::ceil((maxY - minY) / hy_)+1)) {
        
        if (std::fmod((maxX - minX),hx_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along x with given spacing
            std::cout << "Cannot evenly divide with given x spacing, correction occurred on spacing from " << hx_ << "---> " << hx << std::endl;
        }
        if (std::fmod((maxY - minY),hy_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along y with given spacing
            std::cout << "Cannot evenly divide with given y spacing, corrrection occured on spacing from " << hy_ << "---> " << hy << std::endl;
        }
    };

    /// access point operator
    point Mesh2D::operator()(size_t i, size_t j) const {
        if(i-1<ny && j-1<nx){ ///< check correct indexing value
            return coordinate[(i-1)*nx+j-1]; ///! row index span y and column index span x
        }
            throw std::out_of_range("Invalid index!");
        
    };

    /// print mesh information
    std::ostream& operator<<(std::ostream& os, const Mesh2D & mesh) {
        os << static_cast<const Domain &>(mesh) 
            << "Number of points in x direction: " << mesh.getNx() << " , Spacing: " << mesh.getHx() << std::endl
            << "Number of points in y direction: " << mesh.getNy() << " , Spacing: " << mesh.getHy() << std::endl;
        return os;
    };
} 