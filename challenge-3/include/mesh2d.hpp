#ifndef MESH_2D_HPP
#define MESH__2D_HPP

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include "domain2d.hpp"


namespace challenge3 {

    /**
     * @brief Cartesian 2D grid mesh 
    */
    class Mesh2D {
        using point = std::array<double,2>;
    private:
        std::vector<point> points; ///< matrix containing points of the grid
        Domain domain; //< bottom-left and top-right vertex of the rectangular domain 
        size_t nx; //< number of split along x
        size_t ny; //< number of split along y
        double hx; //< x spacing
        double hy; //< y spacing

        /**
         * @brief constructor given the domain and number of points along each direction
        */
        Mesh2D (const Domain & domain_, size_t nx_, size_t ny_);

        /**
         * @brief constructor given the domain and the spacing along each direction
         * @note if no possible evenly divison reduce the spacing to suitable value
        */
        Mesh2D (const Domain & domain_, double hx_, double hy_);

    public:

        Mesh2D() = default; //< default empty constructor

        /// factory methods to select explicityle the strategy of the creation of the mesh
        static Mesh2D createWithPoints(const Domain & domain, size_t nx, size_t ny);
        static Mesh2D createWithSpacing(const Domain & domain, double hx, double hy);

        /// print mesh information 
        friend std::ostream& operator<<(std::ostream& os, const Mesh2D & mesh);


        /// getters
        inline Domain getDomain() const { return domain; }
        inline size_t getNx() const { return nx; }
        inline size_t getNy() const { return ny; }
        inline double getHx() const { return hx; }
        inline double getHy() const { return hy; }
    };
}

/*
inline double getMinX() const { return domain.getMinX(); }; ///< minimum x value of the domain
        inline double getMinY() const { return domain.getMinY(); }; ///< minimum y value of the domain
        inline double getMaxX() const { return domain.getMaxX(); }; ///< maximum x value of the domain
        inline double getMaxY() const { return domain.getMaxY(); }; ///< maximum y value of the domain
*/

#endif