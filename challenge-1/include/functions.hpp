#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>

typedef std::vector<double> Point;

struct Function{
    double operator()(const Point & x) const;
};

struct Gradient{
    Point operator()(const Point & x) const;
};

#endif