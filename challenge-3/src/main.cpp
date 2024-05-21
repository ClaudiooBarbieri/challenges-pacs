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

/// parameters for mesh, function and method option
struct Parameters{
  /// mesh  
  double x0; ///< min x
  double y0; ///< min y
  double xn; ///< max x
  double yn; ///< max y
  size_t n; ///< number of points in each direction
  double h; ///< spacing
  /// method option
  unsigned int maxIter; ///< number of maximum iteration
  double tol; ///< tolerance
  /// function
  std::string f;
};

/// function to read parameters
Parameters readParameters(const std::string & parFileName);

int main(int argc, char * argv[]){
  Parameters p = readParameters("../data/param.json");
  Mesh2D mesh = Mesh2D::createWithPoints({p.x0,p.y0},{p.xn,p.yn},p.n,p.n);
  JacobiSolver js(mesh,p.f,p.maxIter,p.tol);
  js.solve();
  return 0;
}

Parameters readParameters(const std::string & parFileName){
  Parameters parameters; ///< Parameters object to hold the read parameters

  // Reading JSON file with data
  std::ifstream f(parFileName);
  nlohmann::json parFile = nlohmann::json::parse(f);      
  /// read and set values
  parameters.x0 = parFile["mesh"].value("minX", 0);
  parameters.y0 = parFile["mesh"].value("minY", 0);
  parameters.xn = parFile["mesh"].value("maxX", 1);
  parameters.yn = parFile["mesh"].value("maxY", 1);
  parameters.n = parFile["mesh"].value("n", 2);
  parameters.h = parFile["mesh"].value("h", 0.5);
  parameters.maxIter = parFile["option"].value("iteration", 1);
  parameters.tol = parFile["option"].value("tolerance", 1e-3);
  parameters.f = parFile["function"];

  return parameters; ///< Return the read parameters
}