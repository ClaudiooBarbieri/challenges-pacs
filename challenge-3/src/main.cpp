#include <iostream>
#include "domain2d.hpp"
#include "mesh2d.hpp"
#include "JacobiSolver.hpp"
#include "writeVTK.hpp"
#include "json.hpp"

using challenge3::Mesh2D;
using challenge3::point;
using challenge3::JacobiSolver;
using challenge3::Domain;

int main(int argc, char * argv[]){
    Mesh2D mymesh = Mesh2D::createWithPoints({0,0},{1,1},16,16);
    JacobiSolver js(mymesh,"8*_pi^2*sin(2*_pi*x)*sin(2*_pi*y)",1000,1e-6);
    js.solve();
    return 0;
}

/*Parameters readParameters(const std::string & parFileName){
    Parameters parameters; ///< Parameters object to hold the read parameters
    Parameters defaultParameters; ///< Default parameters used when values are not provided in the JSON file

    // Reading JSON file with data
    std::ifstream f(parFileName);
    nlohmann::json parFile = nlohmann::json::parse(f);      

    // Read options
    parameters.maxIter = parFile["option"].value("n_max_iter", defaultParameters.maxIter);
    parameters.resTol = parFile["option"].value("tol_res", defaultParameters.resTol);
    parameters.stepTol = parFile["option"].value("tol_step", defaultParameters.stepTol);
    parameters.alpha0 = parFile["option"].value("alpha0", defaultParameters.alpha0);
    parameters.mu = parFile["option"].value("mu", defaultParameters.mu);
    parameters.sigma = parFile["option"].value("sigma", defaultParameters.sigma);

    // Read starting point
    parameters.x = {parFile["point"].value("x1", 0.0), parFile["point"].value("x2", 0.0)};

    parameters.s = stringToStrategy(parFile["strategy"]); ///< Convert strategy string to enum

    return parameters; ///< Return the read parameters
}*/

/*
  8  9  10 11  20 21 22 23  
  4  5  6  7   10 11 12 13
  0  1  2  3   00 01 02 03
*/