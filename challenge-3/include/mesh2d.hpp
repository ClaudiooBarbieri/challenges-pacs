#ifndef MESH_2D_HPP
#define MESH_2D_HPP

#include "domain2d.hpp"

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>


namespace challenge3 {

    /**
     * @brief Cartesian 2D grid mesh on rectangular domain
     * @note points are stored left to right bottom to top
    */
    class Mesh2D : public Domain2D{

    private:

        std::vector<Point2D> coordinate; ///< coordinates of the Point2Ds in the grid stored left to right, bottom to top 
        size_t nx; //< number of points in x direction
        size_t ny; //< number of points in y direction
        double hx; //< spacing along x
        double hy; //< spacing along y

        /**
         * @brief constructor given ranges of x and y and number of points along each direction
        */
        Mesh2D (double minX_, double maxX_, double minY_, double maxY_, size_t nx_, size_t ny_);

        /**
         * @brief constructor given ranges of x and y and spacing along each direction
         * @note if no possible evenly divison reduce the spacing to suitable value
        */
        Mesh2D (double minX_, double maxX_, double minY_, double maxY_, double hx_, double hy_);

        /**
         * @brief constructor given top right and bottom left points and number of points along each direction
        */
        Mesh2D (const Point2D & tr, const Point2D & bl, size_t nx_, size_t ny_);

        /**
         * @brief constructor given top right and bottom left points and spacing along each direction
         * @note if no possible evenly divison reduce the spacing to suitable value
        */
        Mesh2D (const Point2D & tr, const Point2D & bl, double hx_, double hy_);

    public:

        // ***** factory methods to select explicitly the strategy of the creation of the mesh (ranges) ***** //

        /**
         * @brief factory given ranges and number of points
        */
        static Mesh2D createWithPoints(double minX_, double maxX_, double minY_, double maxY_, size_t nx, size_t ny);

        /**
         * @brief factory given ranges and spacing
        */
        static Mesh2D createWithSpacing(double minX_, double maxX_, double minY_, double maxY_, double hx, double hy);

        // ***** factory methods to select explicitly the strategy of the creation of the mesh (vertices) ***** //

        /**
         * @brief factory given vertices and number of points
        */
        static Mesh2D createWithPoints(const Point2D & tr, const Point2D & bl, size_t nx, size_t ny);

        /**
         * @brief factory given vertices and spacing
        */
        static Mesh2D createWithSpacing(const Point2D & tr, const Point2D & bl, double hx, double hy);

        /**
         * @brief point by its indexes in the mesh
        */
        const Point2D & operator()(size_t i, size_t j) const;

        /**
         * @brief print the mesh
        */
        friend std::ostream& operator<<(std::ostream& os, const Mesh2D & mesh);

        /// getters
        inline size_t getNx() const { return nx; }
        inline size_t getNy() const { return ny; }
        inline double getHx() const { return hx; }
        inline double getHy() const { return hy; }

    };

}

#endif