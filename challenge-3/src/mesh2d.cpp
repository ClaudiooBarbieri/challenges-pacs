#include "mesh2d.hpp"

namespace challenge3 {
    /// factory method to create the mesh given the number of points
    Mesh2D Mesh2D::createWithPoints(const Domain & domain, size_t nx, size_t ny) {
        return Mesh2D(domain, nx, ny);
    };

    /// factory method to create the mesh given the spacing
    Mesh2D Mesh2D::createWithSpacing(const Domain & domain, double hx, double hy) {
        return Mesh2D(domain, hx, hy);
    }

    /// constructor given the domain and number of points along each direction
    Mesh2D::Mesh2D(const Domain & domain_, size_t nx_, size_t ny_) : domain{domain_},nx{nx_},ny{ny_} { 
        hx = (domain.getMaxX()-domain.getMinX())/(nx-1); ///< setting x spacing according to number of points wanted
        hy = (domain.getMaxY()-domain.getMinY())/(ny-1); ///< setting y spacing according to number of points wanted

        /// initialize points in the grid with thei x and y coordinate
        coordinate.reserve(nx*ny);
        for(double y = domain.getMinY() ; y <= domain.getMaxY() ; y+=hy){
            for(double x = domain.getMinX() ; x <= domain.getMaxX() ; x+=hx){
                //std::cout << "(" << x << "," << y << ")" << std::endl;
                coordinate.push_back({x,y});
                //std::cout << "(" << (*(coordinate.end()-1))[0] << "," << (*(coordinate.end()-1))[1] << ")" << std::endl;
            }
        }
    };

    /// constructor given the domain and the spacing along each direction
    Mesh2D::Mesh2D(const Domain & domain_, double hx_, double hy_) : Mesh2D(domain_,static_cast<size_t>(std::ceil((domain_.getMaxX() - domain_.getMinX()) / hx_)+1),
             static_cast<size_t>(std::ceil((domain_.getMaxY() - domain_.getMinY()) / hy_)+1)) {
        
        if (std::fmod((domain.getMaxX() - domain.getMinX()),hx_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along x with given spacing
            std::cout << "Cannot evenly divide with given x spacing, correction occurred on spacing from " << hx_ << "---> " << hx << std::endl;
        }
        if (std::fmod((domain.getMaxY() - domain.getMinY()),hy_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along y with given spacing
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
        os << mesh.getDomain() 
            << "Number of points in x direction: " << mesh.getNx() << " , Spacing: " << mesh.getHx() << std::endl
            << "Number of points in y direction: " << mesh.getNy() << " , Spacing: " << mesh.getHy() << std::endl;
        return os;
    };
} 