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
#include <mpi.h>
#include <omp.h>

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
         * @return solution in vector of vector form
        */
        std::vector<std::vector<double>> solve();

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
        inline void setDim(solution & s) { s.resize(ny,std::vector<double>(nx,0)); };

        /**
         * @brief tie x,y variable to the parser and set its expression
        */
        void setParser(mu::Parser & parser, const std::string & expr);
        
    };

    /**
    * @brief l2 matrix norm of the difference 
    * @note vector has to be considered as matrix
    * @param h factor accounting the spacing
    * @param nR number of rows of underlying matrix
    * @param nC number of columns of underlying matrix
    */
    double norm(const std::vector<double> & v1 , const std::vector<double> & v2, double h, size_t nR, size_t nC);

    /**
    * @brief l2 matrix norm of the difference
    * @param h factor accounting the spacing
    */
    double norm(const std::vector<std::vector<double>> & m1 , const std::vector<std::vector<double>> & m2, double h);

    /**
     * @brief parallel solver
     * @note deal only with zero as boundary condition, cannot set it
     * @param mesh mesh on whic solve the problem
     * @param function function of the problem
     * @param nMax maximum number of iterations
     * @param tol tolerance on the update of the solution
     * @return solution in vectorial form
    */
    std::vector<double> parallelSolve(const Mesh2D & mesh, std::string & function, unsigned int nMax, double tol);

    /**
    * @brief generate VTK file of the solution
    * @param values values in the form of a vector (storing a matrix) to be translated in vtk file
    * @param x0 x coordinate of the origin of the grid
    * @param y0 y coordinate of the origin of the grid
    * @param nx number of points along x direction
    * @param ny number of points along y direction
    * @param hx spacing in x direction
    * @param hy spacing in y direction
    * @param extra string put as "solution"+extra+".vtk"
    */
    void generateVTK(const std::vector<double> & values, double x0, double y0, size_t nx, size_t ny, double hx, double hy,  const std::string & extra = "");

    /**
    * @brief generate VTK file of the solution
    * @param values values in the form of matrix translated in vtk file
    * @param x0 x coordinate of the origin of the grid
    * @param y0 y coordinate of the origin of the grid
    * @param nx number of points along x direction
    * @param ny number of points along y direction
    * @param hx spacing in x direction
    * @param hy spacing in y direction
    * @param extra string put as "solution"+extra+".vtk"
    */
    void generateVTK(const std::vector<std::vector<double>> & sol, double x0, double y0, size_t nx, size_t ny, double hx, double hy,  const std::string & extra = "");

    /**
     * @brief compare the given solution vectorial form to the exact one wrt l2 norm
     * @note do not take in account bondary condition
     * @param mesh on which evaluate
     * @param computed solution
     * @param exact expression of exact solution
     * @param nProcess number of processor used to compute solution
    */
    void compareSolution(const Mesh2D & mesh, const std::vector<double> & solution, const std::string & exact, int nProcess = 1);

    /**
     * @brief compare the given solution amtricial form to the exact one wrt l2 norm
     * @note do not take in account bondary condition
     * @param mesh on which evaluate
     * @param computed solution
     * @param exact expression of exact solution
     * @param nProcess number of processor used to compute solution
    */
    void compareSolution(const Mesh2D & mesh, const std::vector<std::vector<double>> & solution, const std::string & exact, int nProcess = 1);

}

#endif