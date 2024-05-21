#ifndef MESH_2D_HPP
#define MESH_2D_HPP

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include "domain2d.hpp"


namespace challenge3 {

    /**
     * @brief Cartesian 2D grid mesh on the rectangular domain
    */
    class Mesh2D : public Domain{
        using point = std::array<double,2>;
    private:
        std::vector<point> coordinate; ///< coordinates of the points in the grid stored left to right, bottom to top
        //Domain domain; //< bottom-left and top-right vertex of the rectangular domain 
        size_t nx; //< number of split along x
        size_t ny; //< number of split along y
        double hx; //< x spacing
        double hy; //< y spacing

        /**
         * @brief constructor given bottom left and top right points of the domain and number of points along each direction
        */
        Mesh2D (const point & bl, const point & tr, size_t nx_, size_t ny_);

        /**
         * @brief constructor given bottom left and top right points of the domain and the spacing along each direction
         * @note if no possible evenly divison reduce the spacing to suitable value
        */
        Mesh2D (const point & bl, const point & tr, double hx_, double hy_);

    public:

        Mesh2D() = default; //< default empty constructor

        /// factory methods to select explicityle the strategy of the creation of the mesh
        static Mesh2D createWithPoints(const point & bl, const point & tr, size_t nx, size_t ny);
        static Mesh2D createWithSpacing(const point & bl, const point & tr, double hx, double hy);

        /// access point by its indexes
        point operator()(size_t i, size_t j) const;

        /// print mesh information 
        friend std::ostream& operator<<(std::ostream& os, const Mesh2D & mesh);

        /// getters
        inline size_t getNx() const { return nx; }
        inline size_t getNy() const { return ny; }
        inline double getHx() const { return hx; }
        inline double getHy() const { return hy; }
    };
}

#endif