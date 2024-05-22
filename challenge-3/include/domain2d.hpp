#ifndef DOMAIN_2D_HPP
#define DOMAIN_2D_HPP

#include "point2d.hpp"

#include <iostream>
#include <iomanip>

namespace challenge3
{
    
    /**
     * @brief rectangular 2 dimensional domain 
     * @note only proper rectangular can be built, no lines
    */

    class Domain2D {

    protected:
        double minX{0}; ///< minimum x value of the domain
        double maxX{0}; ///< maximum x value of the domain
        double minY{0}; ///< minimum y value of the domain
        double maxY{0}; ///< maximum y value of the domain

    public:

        /**
         * @brief constructor given x range and y range
        */
        Domain2D(double minX_, double maxX_, double minY_, double maxY_);

        /**
         * @brief constructor given top right  bottom left vertices
         * @note if no bottom left corner provided origin is considered
        */
        Domain2D(const Point2D & tr , const Point2D & bl = Point2D(0,0)); 

        // getters
        inline double getMinX() const { return minX; }; 
        inline double getMinY() const { return minY; }; 
        inline double getMaxX() const { return maxX; }; 
        inline double getMaxY() const { return maxY; }; 

        /**
         * @brief print domain four vertices 
        */
        friend std::ostream& operator<<(std::ostream& os, const Domain2D & domain);

    private:

        inline bool isValid(double minX_, double maxX_, double minY_, double maxY_) const { return minX_ < maxX_ && minY_ < maxY_; }; ///< check if given ranges form a valid domain

    };

} 


#endif