#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include <muParser.h>
#include "mesh2d.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "writeVTK.hpp"

namespace challenge3 {
    using solution = std::vector<std::vector<double>>; ///< alias for scalar field of the solution

    /**
     * @brief Jacobi Iteration methoed for solving Laplace equation
     * @note accepts only squared and evenly spaced mesh
    */
    class JacobiSolver {
    private:
        Mesh2D mesh; ///< mesh on which construct the solution of the Laplace equation
        mu::Parser f; ///< function for which solving the Laplace equation
        solution sol; ///< discrete solution of the Laplace equation
        double x; ///< x variable for parser
        double y; ///< y variable for parser
        unsigned int nMax; ///< maximum number of iteration
        double tol; ///< tolerance of the method
        size_t n; ///< points in both x and y direction
        double h; ///< spacing between points

        /// initialize solution with zeros and boundary condition defined in doc 
        void initSolution();

        /// update solution and return the norm of the increment
        double updateSol(); 

        /// matrix L2 norm of the difference 
        double norm(const std::vector<std::vector<double>> & m1 , const std::vector<std::vector<double>> & m2) const;

    public:
        /**
         * @brief solver constructor
         * @param mesh_ mesh on which construct the solution
         * @param expr expression of the function
         * @param nMax_ max number of iteration
         * @param tol_ tolerance
        */
        JacobiSolver(const Mesh2D & mesh_, const std::string & expr, unsigned int nMax_, double tol_);
        
        /// solve method and write VTK file solution
        void solve();

        /// set function given its expression 
        void setFunction(const std::string & expr);

    };
}

#endif