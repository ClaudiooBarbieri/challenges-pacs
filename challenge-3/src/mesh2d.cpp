#include "mesh2d.hpp"

namespace challenge3 {

    Mesh2D Mesh2D::createWithPoints(double minX_, double maxX_, double minY_, double maxY_, size_t nx, size_t ny) {
        return Mesh2D(minX_, maxX_, minY_, maxY_, nx, ny);
    };

    Mesh2D Mesh2D::createWithSpacing(double minX_, double maxX_, double minY_, double maxY_, double hx, double hy) {
        return Mesh2D(minX_, maxX_, minY_, maxY_, hx, hy);
    };

    Mesh2D Mesh2D::createWithPoints(const Point2D & tr, const Point2D & bl, size_t nx, size_t ny) {
        return Mesh2D(tr, bl, nx, ny);
    };

    Mesh2D Mesh2D::createWithSpacing(const Point2D & tr, const Point2D & bl, double hx, double hy) {
        return Mesh2D(tr, bl, hx, hy);
    };

    Mesh2D::Mesh2D(double minX_, double maxX_, double minY_, double maxY_, size_t nx_, size_t ny_) : Domain2D(minX_, maxX_, minY_, maxY_),nx{nx_},ny{ny_} { 

        hx = (maxX-minX)/(nx-1); ///< setting x spacing according to number of points 
        hy = (maxY-minY)/(ny-1); ///< setting y spacing according to number of points

        coordinate.reserve(nx*ny); ///< reserve space storing the points

        /// populate points grid
        for(double y = minY ; y <= maxY ; y+=hy){
            for(double x = minX ; x <= maxX ; x+=hx){
                coordinate.emplace_back(x,y);
            }
        }

    };

    Mesh2D::Mesh2D(double minX_, double maxX_, double minY_, double maxY_, double hx_, double hy_) : 
                Mesh2D(minX_, maxX_, minY_, maxY_, static_cast<size_t>(std::ceil((maxX - minX) / hx_)+1), static_cast<size_t>(std::ceil((maxY - minY) / hy_)+1)) {
        
        if (std::fmod((maxX - minX),hx_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along x with given spacing
            std::cout << "Cannot evenly divide with given x spacing, correction occurred on spacing from " << hx_ << "---> " << hx << std::endl;
        }

        if (std::fmod((maxY - minY),hy_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along y with given spacing
            std::cout << "Cannot evenly divide with given y spacing, corrrection occured on spacing from " << hy_ << "---> " << hy << std::endl;
        }

    };

    Mesh2D::Mesh2D(const Point2D & tr, const Point2D & bl, size_t nx_, size_t ny_) : 
                Mesh2D(bl.getX(), tr.getX(), bl.getY(), tr.getY(), nx_, ny) { 
    };

    Mesh2D::Mesh2D(const Point2D & tr, const Point2D & bl, double hx_, double hy_) : 
                Mesh2D(bl.getX(), tr.getX(), bl.getY(), tr.getY(), hx_, hy_) {
    };

    const Point2D & Mesh2D::operator()(size_t i, size_t j) const {

        if(i>=ny && j>=nx){ ///< check correct indexing value
           throw std::out_of_range("Invalid index!");
        }

        return coordinate[i*nx+j]; //! row index span y and column index span x

    };

    std::ostream& operator<<(std::ostream& os, const Mesh2D & mesh) {

        os << std::fixed <<std::setprecision(2) << static_cast<const Domain2D &>(mesh) 
            << "Number of points in x direction: " << mesh.getNx() << " , Spacing: " << mesh.getHx() << std::endl
            << "Number of points in y direction: " << mesh.getNy() << " , Spacing: " << mesh.getHy() << std::endl << std::endl;

        std::cout << std::setw(3);

        /// print coordinates point in the mesh
        for(int i = mesh.ny-1 ; i >= 0 ; --i){
            for(int j = 0 ; j < mesh.nx ; ++j){
                std::cout << "(" << mesh(i,j).getX() << "," <<  mesh(i,j).getY() << ") " << std::setw(3);
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;

        return os;
        
    };
} 