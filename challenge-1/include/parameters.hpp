#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <string>

typedef std::vector<double> Point;

enum class Strategy {
    EXPONENTIAL,
    INVERSE,
    ARMIJO,
    NONE
};

struct Parameters{
    unsigned int maxIter = 100;
    double resTol = 1e-3;
    double stepTol= 1e-3;
    double alpha0 = 0.5;
    double mu = 0.5;
    double sigma = 0.3;
    Point x{0,0};
    Strategy s = Strategy::EXPONENTIAL;
    void print() const;
};

Parameters readParameters(const std::string & parFileName);

std::string strategyToString(const Strategy & s);

Strategy stringToStrategy(const std::string & str);

#endif