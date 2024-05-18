#include <iostream>
#include "domain2d.hpp"
#include "mesh2d.hpp"
#include "JacobiSolver.hpp"

using challenge3::Mesh2D;
using challenge3::point;
using challenge3::JacobiSolver;
using challenge3::Domain;

int main(int argc, char * argv[]){
    Domain mydomain({0,0},{5,10});
    Mesh2D mymesh = Mesh2D::createWithSpacing(mydomain,2,2);
    std::cout << mymesh;
    return 0;
}