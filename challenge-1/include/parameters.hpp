#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <string>

/// Type definition for a point represented as a vector of doubles
typedef std::vector<double> Point;

/// Enumeration representing optimization strategies
enum class Strategy {
    EXPONENTIAL, ///< Exponential strategy
    INVERSE,     ///< Inverse strategy
    ARMIJO,      ///< Armijo strategy
    NONE         ///< No strategy (invalid)
};

/// Struct representing optimization parameters
struct Parameters{
    unsigned int maxIter = 100; ///< Maximum number of iterations
    double resTol = 1e-3;       ///< Tolerance on the residual
    double stepTol = 1e-3;      ///< Tolerance on the steps
    double alpha0 = 0.5;        ///< Initial step size
    double mu = 0.5;            ///< Decay parameter for step size
    double sigma = 0.3;         ///< Constant parameter for Armijo rule
    Point x{0,0};               ///< Starting point
    Strategy s = Strategy::EXPONENTIAL; ///< Optimization strategy
    
    /// Print parameters to standard output.

    void print() const;
};

/// Read parameters from a JSON file.

Parameters readParameters(const std::string & parFileName);

/// Convert Strategy enum to string.

std::string strategyToString(const Strategy & s);

/// Convert string to Strategy enum.

Strategy stringToStrategy(const std::string & str);

#endif
