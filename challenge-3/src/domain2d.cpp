#include "domain2d.hpp"

namespace challenge3 {

    /// print domain four vertices
    std::ostream& operator<<(std::ostream& os, const Domain & domain) {
        os << "Vertices of the domain: " 
            << "(" << domain.getMinX()<< "," << domain.getMinY() << ") - " 
            << "(" << domain.getMaxX() << "," << domain.getMinY() << ") - "
            << "(" << domain.getMaxX() << "," << domain.getMaxY() << ") - " 
            << "(" << domain.getMinX() << "," << domain.getMaxY() << ")"
            << std::endl;
        return os;
    };

    /// constructor with bottom-left top-right points
    Domain::Domain(const point &bl, const point &tr) : bottomLeft{bl},topRight{tr} {
        minX = bottomLeft[0];  ///< minimum x value of the domain
        minY = bottomLeft[1]; ///< minimum y value of the domain
        maxX = topRight[0];  ///< maximum x value of the domain
        maxY = topRight[1];  ///< maximum y value of the domain
    };

}