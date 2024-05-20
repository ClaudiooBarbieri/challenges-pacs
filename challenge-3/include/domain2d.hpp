#ifndef DOMAIN_2D_HPP
#define DOMAIN_2D_HPP

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

namespace challenge3
{
    using point = std::array<double,2>;
    /**
     * @brief rectangular domain defined by its bottom left and top right points
    */
    class Domain {
        using point = std::array<double,2>;
    private:
        point  bottomLeft{0,0}; ///< bottom left point corner of the rectangular domain
        point  topRight{0,0}; ///< bottom left point corner of the rectangular domain
    protected:
        double minX{0}; ///< minimum x value of the domain
        double maxX{0}; ///< maximum x value of the domain
        double minY{0}; ///< minimum y value of the domain
        double maxY{0}; ///< maximum y value of the domain
    public:
        Domain() = default;
        Domain(const point & bl, const point & tr);

        inline double getMinX() const { return minX; }; ///< minimum x value of the domain
        inline double getMinY() const { return minY; }; ///< minimum y value of the domain
        inline double getMaxX() const { return maxX; }; ///< maximum x value of the domain
        inline double getMaxY() const { return maxY; }; ///< maximum y value of the domain

        /// print domain four vertices counterclockwise
        friend std::ostream& operator<<(std::ostream& os, const Domain & domain);
    };
} 


#endif