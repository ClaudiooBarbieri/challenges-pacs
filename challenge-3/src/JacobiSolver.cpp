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

        setDim(sol); ///< adjust dimension

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
        setDim(newSol);
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
            e = norm(newSol,sol,std::max(hx,hy)); ///< get the norm of the increment of updated solution, use h as th maximum between the spacings
        }

        generateVTK(sol,mesh.getMinX(),mesh.getMinY(),nx,ny,hx,hy); ///< generate VTK file

    };


    double norm(const std::vector<std::vector<double>> & m1 , const std::vector<std::vector<double>> & m2, double h) {

        size_t nRows = m1.size();

        if(nRows != m2.size()) ///< check row size compatibility
            throw std::invalid_argument("norm size mismatch - row");

        double ss{0.}; ///< init SS

        for(size_t i = 0 ; i < m1.size() ; ++i){

            size_t nCols = m1[i].size(); //

            if(nCols != m2[i].size())  ///< check col size compatibility
                throw std::invalid_argument("norm size mismatch - col");

            for(size_t j = 0 ; j < nCols ; ++j){
               ss += (m1[i][j]-m2[i][j])*(m1[i][j]-m2[i][j]); ///< sum of SS 
            }
        }

        return sqrt(h*ss); ///< multiply by h and square the SS according to instruction in ../doc

    };

    void parallelSolve(const Mesh2D & mesh, std::string & function, unsigned int nMax, double tol){

        // setup MPI variables
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        //std::cout << "rank " << rank << " mesh " << mesh << std::endl;

        // set number of points and spacing
        size_t nx = mesh.getNx();
        double hx = mesh.getHx();
        size_t ny = mesh.getNy();
        double hy = mesh.getHy();

        // if not at least two row per process
        if(ny / size < 2){
            std::cerr << "Too much processes wrt to rows, at least two rows at each process" << std::endl;
            return;
        }

        //std::cout << "rank " << rank << " nx " << nx << std::endl;
        //std::cout << "rank " << rank << " ny " << hx << std::endl;
        //std::cout << "rank " << rank << " hx " << ny << std::endl;
        //std::cout << "rank " << rank << " hy " << hy << std::endl;

        mu::Parser f; ///< parser of the function

        // parser variable
        double x;
        double y;

        // set the parser
        try {
            
            // set parser references variable
            f.DefineVar("x", &x); 
            f.DefineVar("y", &y);

            f.SetExpr(function); ///< set parser expression

        }
        catch (mu::Parser::exception_type &e) {

            std::cerr << "Error: " << e.GetMsg() << std::endl;

        }

        int allConverged = 0; ///< to check convergence of all processes

        int nBaseRows = ny / size; ///< common number of rows for process, even spread part
        //std::cout << "Process: " << rank << "- nabserows :" << nBaseRows << std::endl;
        int nExtraRows = ny % size; ///< extra number of rows to be accounted
        //std::cout << "Process: " << rank << "- nextraows :" << nExtraRows << std::endl;

        int firstRow = rank * nBaseRows + std::min(rank, nExtraRows); ///< set first row, consider the number of extra rows already assigned
        //std::cout << "Process: " << rank << "- firstrow :" << firstRow << std::endl;
        int lastRow = firstRow + ((nExtraRows > rank) ? nBaseRows + 1 : nBaseRows) - 1; ///< last row considring if it has extra rows to be considered
        //std::cout << "Process: " << rank << "- lastrow :" << lastRow << std::endl;

        int localRows = lastRow - firstRow + 1; ///< number of local rows
        //std::cout << "Process: " << rank << "- localrows :" << localRows << std::endl;

        std::vector<std::vector<double>> localSol; ///local solution
        localSol.resize(localRows, std::vector<double>(nx,0)); ///< size solution
        std::vector<std::vector<double>> localNewSol; ///< local new solution for update
        localNewSol.resize(localRows, std::vector<double>(nx,0)); ///< size new solution

        std::vector<double> belowRow(nx,0); ///< for receiving row below evaluated by previous rank process
        //std::cout << "Process: " << rank << "- belowsize :" << belowRow.size() << std::endl;
        std::vector<double> aboveRow(nx,0); ///< for receiving row above evaluated by succesive rank process
        //std::cout << "Process: " << rank << "- above :" << aboveRow.size() << std::endl;

        // Jacobi iteration
        for(size_t k = 0; k < nMax && allConverged != size ; ++k){

            //std::cout << "rank " << rank << " entered loop " << std::endl;

            allConverged = 0; ///< reset convergence flag

            if(rank != size - 1){ ///< if not last process send solution of last of your row and receive the first row of the following process
                MPI_Sendrecv(localSol[localRows-1].data(), nx, MPI_DOUBLE, rank+1, 0, aboveRow.data(), nx, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //std::cout << "rank " << rank << " size sol above " << localSol[0].size() << std::endl;
            }
            if(rank != 0){ ///< if not first process send solution of first of your row and receive the last row of the previous process
                MPI_Sendrecv(localSol[0].data(), nx, MPI_DOUBLE, rank-1, 0, belowRow.data(), nx, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //std::cout << "rank " << rank << " size sol below " << localSol[0].size() << std::endl;
            }

            // last row to be considered for non last rank process (boundary)
            if(rank != size - 1){
                int i = localRows - 1;
                //std::cout << "rank " << rank << " last row index " << i << std::endl;
                for(size_t j = 1 ; j < nx-1; ++j){

                    x = mesh(lastRow,j).getX(); ///< coordinate x value last row mesh correspondence
                    y = mesh(lastRow,j).getY(); ///< coordinate y value last row mesh correspondence

                    //std::cout << "Process: " << rank << "- new sol above" << std::endl;
                    localNewSol[i][j] = .25*(aboveRow[j]+localSol[i-1][j]+localSol[i][j+1]+localSol[i][j-1]+hx*hy*f.Eval()); ///<  four point stencil
                }
            }
            
            // first row to be considered for non first process (boundary)
            if(rank != 0){
                int i = 0;
                //std::cout << "rank " << rank << " first row index " << i << std::endl;
                for(size_t j = 1 ; j < nx-1; ++j){

                    x = mesh(firstRow,j).getX(); ///< coordinate x value first row mesh correspondence
                    y = mesh(firstRow,j).getY(); ///< coordinate y value first row mesh correspondence

                    //std::cout << "Process: " << rank << "- new sol below" << std::endl;
                    localNewSol[i][j] = .25*(localSol[i+1][j]+belowRow[j]+localSol[i][j+1]+localSol[i][j-1]+hx*hy*f.Eval()); ///<  four point stencil
                }
            }
        
            // middle rows
            for(size_t i = 1 ; i < localRows - 1 ; ++i){
                //std::cout << "rank " << rank << " middle" << i << std::endl;

                for(size_t j = 1 ; j < nx-1; ++j){

                    x = mesh(firstRow + i,j).getX(); ///< coordinate x value on the mesh, offset wrt to first row of the current process
                    y = mesh(firstRow + i,j).getY(); ///< coordinate y value on the mesh, offset wrt to first row of the current process

                    //std::cout << "Process: " << rank << "- middle new sol" << std::endl;
                    localNewSol[i][j] = .25*(localSol[i+1][j]+localSol[i-1][j]+localSol[i][j+1]+localSol[i][j-1]+hx*hy*f.Eval()); ///<  four point stencil
                }
            }

            std::swap(localSol,localNewSol); ///< set current solution to the one calculated
            //std::cout << "Process: " << rank << "- sol swapped" << std::endl;

            if(norm(localNewSol,localSol,std::max(hx,hy)) < tol){ ///< check local tolerance convergence
                //std::cout << "Process: " << rank << "- convergence" << std::endl;
                allConverged = 1;
            }

            MPI_Barrier(MPI_COMM_WORLD); ///< wait for tolerance of all processes
            
            MPI_Allreduce(MPI_IN_PLACE, &allConverged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); ///< sum to consider how much processes converged
            //std::cout << "Process: " << rank << "- allreduce" << std::endl;
            
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        // rank 0 assemble the global solution and generate VTK
        std::vector<int> nReceive(size,0); ///< vector containing number of rows to be received from rank [.]
        //std::cout << "Process: " << rank << "- nreceive " << nReceive.size() << std::endl;
        std::vector<int> vDisplacement(size,0); ///< vector containing displacement of where to start writing recived row from rank [.]
        //std::cout << "Process: " << rank << "- displacement " << nReceive.size() << std::endl;
        
        int nElem = localRows * nx; ///< total number of element to send
        //std::cout << "Process: " << rank << "- nel " << nElem << std::endl;
        int displacement = firstRow * nx; ///< evaluate relative displacement
        //std::cout << "Process: " << rank << "- displc " << displacement << std::endl;
        
        MPI_Gather(&nElem, 1, MPI_INT, nReceive.data(), 1, MPI_INT, 0, MPI_COMM_WORLD); ///< gather the number of rows to receive for each rank
        MPI_Gather(&displacement, 1, MPI_INT, vDisplacement.data(), 1, MPI_INT, 0, MPI_COMM_WORLD); ///< gather the number of rows to receive for each rank

        std::vector<double>  globalSol;
        if(rank == 0)
            globalSol.resize(nx*ny);

        // flatten local solution into a vector for gathering
        std::vector<double> localS;
        for (const auto& row : localSol) 
            localS.insert(localS.end(), row.begin(), row.end());
  
        MPI_Gatherv(localS.data(), nElem, MPI_DOUBLE, globalSol.data(), nReceive.data(), vDisplacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if(rank == 0){

            std::vector<std::vector<double>> vecOfVecs;
            size_t subSize = nx;

            for (size_t i = 0; i <globalSol.size(); i += subSize) {
                std::vector<double> subVec(globalSol.begin() + i, globalSol.begin() + std::min(i + subSize, globalSol.size()));
                vecOfVecs.push_back(subVec);
            }
    
        }
            
    };

    void generateVTK(const std::vector<std::vector<double>> & values, double x0, double y0, size_t nx, size_t ny, double hx, double hy,  const std::string & extra){

        std::string filename = "../VTK/solution"+extra+".vtk"; ///< name of the solution file

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


        // Write the points and scalar values
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << x0 + i * hx << " " << y0 + j * hy << " " << values[i][j] << "\n";
            }
        }

        // include scalar field for coloring solution
        svtkFile << "POINT_DATA " << nx * ny << "\n";
        svtkFile << "SCALARS solution double 1\n";
        svtkFile << "LOOKUP_TABLE default\n";
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << values[i][j] << "\n";
            }
        }

        std::cout << "VTK produced matrix" << std::endl;
        svtkFile.close();

    };

    void generateVTK(const std::vector<double> & values, double x0, double y0, size_t nx, size_t ny, double hx, double hy,  const std::string & extra){

        std::string filename = "../VTK/solution"+extra+".vtk"; ///< name of the solution file

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


        // Write the points and scalar values
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << x0 + i * hx << " " << y0 + j * hy << " " << values[i*nx+j] << "\n";
            }
        }

        // include scalar field for coloring solution
        svtkFile << "POINT_DATA " << nx * ny << "\n";
        svtkFile << "SCALARS solution double 1\n";
        svtkFile << "LOOKUP_TABLE default\n";
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << values[i*nx+j] << "\n";
            }
        }

        std::cout << "VTK produced vector" << std::endl;
        svtkFile.close();

    };


}