#include <vector>
#include "functions.hpp"

/// Type definition for a point represented as a vector of doubles

typedef std::vector<double> Point;

/// Call operator to evaluate the function.

double Function::operator()(const Point & x) const {
    // x1*x2 + 4*(x1^4) + (x2)^2 + 3*x1
    return x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0];
}

/// Call operator to evaluate the gradient of the function.

Point Gradient::operator()(const Point & x) const{
    // (x2 + 16*(x1)^3 + 3 , x1 + 2*x2)
    return {x[1]+16*x[0]*x[0]*x[0]+3,x[0]+2*x[1]};
}