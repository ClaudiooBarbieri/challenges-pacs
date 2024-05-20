#include <iostream>
#include "domain2d.hpp"
#include "mesh2d.hpp"
#include "JacobiSolver.hpp"

using challenge3::Mesh2D;
using challenge3::point;
using challenge3::JacobiSolver;
using challenge3::Domain;

int main(int argc, char * argv[]){
    Mesh2D mymesh = Mesh2D::createWithPoints({0,0},{3,2},4,3);
    std::cout << mymesh;
    std::cout << mymesh(2,3)[0] << "," << mymesh(2,3)[1] << std::endl;
    return 0;
}

/*

  8  9  10 11  20 21 22 23  
  4  5  6  7   10 11 12 13
  0  1  2  3   00 01 02 03
*/