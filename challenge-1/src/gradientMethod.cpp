#include <vector>
#include <cmath>
#include <iostream>
#include <utility>
#include "gradientMethod.hpp"
#include "parameters.hpp"
#include "functions.hpp"

/// Type definition for a point represented as a vector of doubles
typedef std::vector<double> Point;

/// Template function to minimize a function using various optimization strategies.
/**
 * This function minimizes a given function using one of the specified optimization strategies.
 * @tparam S Strategy enum specifying the optimization strategy to be used (EXPONENTIAL, INVERSE, ARMIJO)
 * @param parameters Parameters struct containing optimization parameters
 * @param f Function object representing the function to be minimized
 * @param df Gradient object representing the gradient of the function
 * @return A pair containing the minimum point and the number of iterations taken to converge
 */
template<Strategy S>
std::pair<Point,unsigned int> argmin(const Parameters & parameters, const Function & f , const Gradient & df){
    unsigned int iter = 0; ///< Iteration counter
    double alpha = parameters.alpha0; ///< Initial step size
    Point xNew = parameters.x; ///< Initialize the new point as starting point to have the update in the loop correct
    Point x{0,0}; ///< Initialization

    // Main optimization loop
    do {
        x = xNew; ///< Update the current point
        xNew[0] = x[0] - alpha * df(x)[0]; ///< Update the first coordinate of the new point
        xNew[1] = x[1] - alpha * df(x)[1]; ///< Update the second coordinate of the new point
        ++iter; ///< Increment the iteration counter

        // Determine the step size based on the selected strategy
        if constexpr (S == Strategy::EXPONENTIAL) {
            alpha = expDecay(parameters, iter);
        } else if constexpr (S == Strategy::INVERSE) {
            alpha = invDecay(parameters, iter);
        } else if constexpr (S == Strategy::ARMIJO) {
            alpha = lineSearch(parameters, f, df, x);
        } else {
            std::cout << "INVALID STRATEGY!" << std::endl;
            return {x, iter};
        }
    } while(!(iter > parameters.maxIter || norm(xNew, x) < parameters.stepTol || std::fabs(f(xNew) - f(x)) < parameters.resTol));

    return {xNew, iter}; ///< Return the minimum point and the number of iterations
}

/// Template instantiation for ARMIJO strategy
template std::pair<Point, unsigned int> argmin<Strategy::ARMIJO>(const Parameters & parameters, const Function & f , const Gradient & df);
/// Template instantiation for EXPONENTIAL strategy
template std::pair<Point, unsigned int> argmin<Strategy::EXPONENTIAL>(const Parameters & parameters, const Function & f , const Gradient & df);
/// Template instantiation for INVERSE strategy
template std::pair<Point, unsigned int> argmin<Strategy::INVERSE>(const Parameters & parameters, const Function & f , const Gradient & df);

/// Exponential decay function for determining step size
double expDecay(const Parameters & parameters, unsigned int k) {
    return parameters.alpha0 * exp(-parameters.mu * k);
}

/// Inverse decay function for determining step size
double invDecay(const Parameters & parameters, unsigned int k) {
    return parameters.alpha0 / (1 + parameters.mu * k);
}

/// Line search function using Armijo rule
double lineSearch(const Parameters & parameters, const Function & f , const Gradient & df, const Point & xk) {
    double alphak = parameters.alpha0; ///< Initial step size
    while (!ArmijoRule(f, df, xk, alphak, parameters.sigma)) {
        alphak /= 2; ///< Reduce step size
    }
    return alphak; ///< Return the determined step size
}

/// Armijo rule for line search
bool ArmijoRule(const Function & f , const Gradient & df, const Point & xk, double alpha, double sigma) {
    return f(xk) - f({xk[0] - alpha * df(xk)[0], xk[1] - alpha * df(xk)[1]}) >= sigma * alpha * std::pow(norm(df(xk)), 2);
}

/// Euclidean distance (L2 norm) between two points
double norm(const Point & x, const Point & y) {
    return std::sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]));
}
