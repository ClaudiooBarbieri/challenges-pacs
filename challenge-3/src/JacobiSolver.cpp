#include "JacobiSolver.hpp"

namespace challenge3{

    JacobiSolver::JacobiSolver(const Mesh2D & mesh_, const std::string & expr, size_t nMax_, double tol_) : mesh{mesh_},nMax{nMax_},tol{tol_} {
        /// check mesh condition
        if (mesh.getHx()!=mesh.getHy() || mesh.getNx()!=mesh.getNy()) {
            throw std::invalid_argument("Not squared and/or evenly spaced mesh");
        }
        /// set number of points and spacing
        n = mesh.getNx();
        h = mesh.getHx();
        /// define parser variable
        f.DefineVar("x", &x); 
        f.DefineVar("y", &y);
        /// define parser expression
        setFunction(expr);
        initSolution();
    }

    void JacobiSolver::initSolution(){
        /// solution size according to the mesh
        sol.resize(n);
        for(auto & row : sol)
            row.resize(n);
    }

    void JacobiSolver::solve(){
        double e = tol;
        for(size_t k = 0 ; k < nMax && e >= tol; ++k){
            e = updateSol(); ///< get the norm of the increment of updated solution
        }
        generateVTKFile("../VTK/solution.vtk",mesh,sol,n,n,h,h);
    }

    double JacobiSolver::updateSol(){
        /// init matrix to contain new solution
        solution newSol;
        newSol.resize(n);
        for(auto & row : newSol)
            row.resize(n);
        /// Jacobi method for new solution k+1
        for(size_t i = 1 ; i < n-1 ; ++i){
            for(size_t j = 1 ; j < n-1; ++j){
                x = mesh(i,j)[0]; ///< coordinate x value on the mesh
                y = mesh(i,j)[1]; ///< coordinate y value on the mesh
                newSol[i][j] = .25*(sol[i+1][j]+sol[i-1][j]+sol[i][j+1]+sol[i][j-1]+h*h*f.Eval()); ///< Jacobi iteration formula, four point stencil average
            }
        }
        std::swap(sol,newSol); ///< set current solution to the one calculated
        return JacobiSolver::incrNorm(newSol,sol);
    }

    double JacobiSolver::incrNorm(const solution & sol1 , const solution & sol2) const {
        double ss{0.};
        for(size_t i = 0 ; i < n ; ++i){
            for(size_t j = 0 ; j < n; ++j){
               ss+=(sol1[i][j]-sol2[i][j])*(sol1[i][j]-sol2[i][j]);
            }
        }
        return sqrt(h*ss);
    }

    void JacobiSolver::setFunction(const std::string & expr){
        try {
            f.SetExpr(expr);
        }
        catch (mu::Parser::exception_type &e) {
            std::cerr << "Error: " << e.GetMsg() << std::endl;
        }
    }
}