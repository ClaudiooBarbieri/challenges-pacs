#include "JacobiSolver.hpp"

namespace challenge3{

    JacobiSolver::JacobiSolver(const Mesh2D & mesh_, const std::string & function, unsigned int nMax_, double tol_, const std::string & bCond) : mesh{mesh_},nMax{nMax_},tol{tol_} {
        
        /// set number of points and spacing
        nx = mesh.getNx();
        hx = mesh.getHx();
        ny = mesh.getNy();
        hy = mesh.getHy();

        /// define parsers expression
        setParser(this->f,function);
        setParser(this->fBound,bCond);

        initSolution(); ///< initialize the solution matrix
    };

    void JacobiSolver::setParser(mu::Parser & parser, const std::string & expr){
        try {
            
            /// set parser references variable
            parser.DefineVar("x", &x); 
            parser.DefineVar("y", &y);

            parser.SetExpr(expr); ///< set parser expression

        }
        catch (mu::Parser::exception_type &e) {

            std::cerr << "Error: " << e.GetMsg() << std::endl;

        }
    };

    void JacobiSolver::initSolution(){

        initDim(sol); ///< adjust dimension

        // set boundary condition

        for(size_t j = 0 ; j < nx ; ++j){
            x = mesh(0,j).getX();
            y = mesh(0,j).getY();
            sol[0][j] = fBound.Eval(); ///< first row
            x = mesh(ny-1,j).getX();
            y = mesh(ny-1,j).getY();
            sol[ny-1][j] = fBound.Eval();; ///< last row
        }

        for(size_t i = 0 ; i < ny ; ++i){
            x = mesh(i,0).getX();
            y = mesh(i,0).getY();
            sol[i][0] = fBound.Eval(); ///< first column
            x = mesh(i,nx-1).getX();
            y = mesh(i,nx-1).getY();
            sol[i][nx-1] = fBound.Eval(); ///< last column
        }    

    };

    void JacobiSolver::solve(){
        
        solution newSol;
        initDim(newSol);
        double e = tol;

        // Jacobi iteration
        for(size_t k = 0 ; k < nMax && e >= tol; ++k){
            for(size_t i = 1 ; i < ny-1 ; ++i){
                for(size_t j = 1 ; j < nx-1; ++j){
                    x = mesh(i,j).getX(); ///< coordinate x value on the mesh
                    y = mesh(i,j).getY(); ///< coordinate y value on the mesh
                    newSol[i][j] = .25*(sol[i+1][j]+sol[i-1][j]+sol[i][j+1]+sol[i][j-1]+hx*hy*f.Eval()); ///<  four point stencil
                }
            }
            std::swap(sol,newSol); ///< set current solution to the one calculated
            e = JacobiSolver::norm(newSol,sol); ///< get the norm of the increment of updated solution
        }

        generateVTK(); ///< generate VTK file

    };

    void JacobiSolver::initDim(solution & s){

        s.resize(ny); ///< number of rows
        for(auto & row : s)
            row.resize(nx); ///< number of columns
        
    }

    double JacobiSolver::norm(const std::vector<std::vector<double>> & m1 , const std::vector<std::vector<double>> & m2) const {

        double ss{0.}; ///< init SS

        for(size_t i = 0 ; i < ny ; ++i){
            for(size_t j = 0 ; j < nx ; ++j){
               ss += (m1[i][j]-m2[i][j])*(m1[i][j]-m2[i][j]); ///< sum of SS 
            }
        }

        return sqrt(std::max(hx,hy)*ss); ///< multiply by h and square the SS according to instruction in ../doc

    };

    void JacobiSolver::generateVTK() const {

        std::string filename = "../VTK/solution.vtk"; ///< name of the solution file

        std::ofstream svtkFile(filename); ///< open file

        // check if the files were opened
        if (!svtkFile.is_open()) {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }

        // write VTK header 
        svtkFile << "# vtk DataFile Version 3.0\n";
        svtkFile << "Discrete solution\n";
        svtkFile << "ASCII\n";                                

        // write grid data 
        svtkFile << "DATASET STRUCTURED_GRID\n";
        svtkFile << "DIMENSIONS " << nx << " " << ny << " 1\n";
        svtkFile << "POINTS " << nx * ny << " double\n";

        // grid origin
        double x0 = mesh.getMinX();
        double y0 = mesh.getMinY();

        // Write the points and scalar values
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << x0 + i * hx << " " << y0 + j * hy << " " << sol[i][j] << "\n";
            }
        }

        // include scalar field for coloring solution
        svtkFile << "POINT_DATA " << nx * ny << "\n";
        svtkFile << "SCALARS solution double 1\n";
        svtkFile << "LOOKUP_TABLE default\n";
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << sol[i][j] << "\n";
            }
        }

        std::cout << "VTK produced" << std::endl;
        svtkFile.close();

    };

}