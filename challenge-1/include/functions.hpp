#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>

/// Type definition for a point represented as a vector of doubles

typedef std::vector<double> Point;

/// Struct representing a mathematical function.

struct Function{
    /// Call operator to evaluate the function.
    double operator()(const Point & x) const;
};

/// Struct representing the gradient of a mathematical function.

struct Gradient{
    /// Call operator to evaluate the gradient of the function.
    Point operator()(const Point & x) const;
};

#endif