#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include "mesh2d.hpp"

#include <muParser.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>

namespace challenge3 {

    /**
     * @brief Jacobi Iteration methoed for solving Laplace equation
     * @note accepts only squared and evenly spaced mesh
    */
    class JacobiSolver {

        using solution = std::vector<std::vector<double>>; ///< discrete solution matrix

    private:

        Mesh2D mesh; ///< mesh on which construct the solution
        solution sol; ///< discrete solution 

        mu::Parser f; ///< function for which solving
        mu::Parser fBound; ///< function used for setting boundary condition
        double x; ///< x variable for parser
        double y; ///< y variable for parser

        unsigned int nMax; ///< maximum number of iteration
        double tol; ///< tolerance of the method

        size_t nx; ///< points in x direction
        size_t ny; ///< points in y direction
        double hx; ///< x spacing between points
        double hy; ///< y spacing between points

    public:
        /**
         * @brief solver constructor
         * @param mesh_ mesh on which construct the solution
         * @param function expression of the function
         * @param bCond expression for the function used for boundary condition
         * @param nMax_ max number of iteration
         * @param tol_ tolerance
        */
        JacobiSolver(const Mesh2D & mesh_, const std::string & function, unsigned int nMax_, double tol_ , const std::string & bCond = "0");
        
        /**
         * @brief solver method
         * @note in stencil formula the h^2 in ../doc documentation become hx*hy
         * @note automatically generate the VTK file format of it
        */
        void solve();

        /**
         * @brief solution getter
        */
        inline solution getSolution() { return sol;}

    private:

        /**
         * @brief initialize solution to zero and set Dirichlet boundary condition
         * 
        */
        void initSolution();

        /**
         * @brief set solution size according to mesh dimension
        */
        void initDim(solution & s);

        /**
         * @brief tie x,y variable to the parser and set its expression
        */
        void setParser(mu::Parser & parser, const std::string & expr);

        /**
         * @brief l2 matrix norm of the difference
         * @note using the maximum spacing as h in given formula in ../doc
        */
        double norm(const std::vector<std::vector<double>> & m1 , const std::vector<std::vector<double>> & m2) const;

        /**
         * @brief generate VTK file of the solution
        */
        void generateVTK() const;

    };
}

#endif