#include "domain2d.hpp"

namespace challenge3 {

    Domain2D::Domain2D(double minX_, double maxX_, double minY_, double maxY_){

        if(isValid(minX_, maxX_, minY_, maxY_)){ ///< check that is a proper rectangular domain
            minX = minX_;
            maxX = maxX_;
            minY = minY_;
            maxY = maxY_;
        }
        else{
            throw std::invalid_argument("Not a proper rectangular domain!");
        } 

    };

    Domain2D::Domain2D(const Point2D & tr , const Point2D & bl) : Domain2D(bl.getX(), tr.getX(), bl.getY(), tr.getY()) {}

    std::ostream& operator<<(std::ostream& os, const Domain2D & domain) {

        os << std::fixed << std::setprecision(2) << "Domain vertices: " << std::endl
            << "(" << domain.minX << "," << domain.maxY << ") ----- " 
            << "(" << domain.maxX << "," << domain.maxY << ")" << std::endl
            << " | " << std::setw(23) << " | " << std::endl
            << "(" << domain.getMaxX() << "," << domain.getMaxY() << ") ----- " 
            << "(" << domain.getMinX() << "," << domain.getMaxY() << ")"
            << std::endl << std::endl;
        return os;

    };

}