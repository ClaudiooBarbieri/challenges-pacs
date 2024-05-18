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

}