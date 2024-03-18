#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <string>

typedef std::vector<double> Point;

struct Parameters{
    unsigned int maxIter = 100;
    double resTol = 1e-3;
    double stepTol= 1e-3;
    double alpha0 = 0.5;
    double mu = 0.5;
    double sigma = 0.3;
    Point x{0,0};
    void print() const;
};

Parameters readParameters(const std::string & parFileName);

#endif