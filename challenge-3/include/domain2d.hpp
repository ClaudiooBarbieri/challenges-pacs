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
     * @brief define the rectangular domain by its bottom left and top right points
    */
    class Domain{
        using point = std::array<double,2>;
    private:
        point  bottomLeft{0,0}; ///< bottom left point corner of the rectangular domain
        point  topRight{0,0}; ///< bottom left point corner of the rectangular domain
    public:
        Domain() = default;
        Domain(const point & bl, const point & tr) : bottomLeft{bl},topRight{tr} {};

        inline double getMinX() const { return bottomLeft[0]; }; ///< minimum x value of the domain
        inline double getMinY() const { return bottomLeft[1]; }; ///< minimum y value of the domain
        inline double getMaxX() const { return topRight[0]; }; ///< maximum x value of the domain
        inline double getMaxY() const { return topRight[1]; }; ///< maximum y value of the domain

        /// print domain four vertices counterclockwise
        friend std::ostream& operator<<(std::ostream& os, const Domain & domain);
    };
} 


#endif