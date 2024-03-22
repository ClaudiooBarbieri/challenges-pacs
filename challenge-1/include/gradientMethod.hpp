#ifndef GRADIENT_METHOD_HPP
#define GRADIENT_METHOD_HPP

#include <vector>
#include <utility>
#include "parameters.hpp"
#include "functions.hpp"

/// Type definition for a point represented as a vector of doubles
typedef std::vector<double> Point;

/// Template function to minimize a function using various optimization strategies.

template<Strategy S>
std::pair<Point,unsigned int> argmin(const Parameters & parameters, const Function & f , const Gradient & df);

/// Extern template instantiation for ARMIJO strategy
extern template std::pair<Point, unsigned int> argmin<Strategy::ARMIJO>(const Parameters & parameters, const Function & f , const Gradient & df);
/// Extern template instantiation for EXPONENTIAL strategy
extern template std::pair<Point, unsigned int> argmin<Strategy::EXPONENTIAL>(const Parameters & parameters, const Function & f , const Gradient & df);
/// Extern template instantiation for INVERSE strategy
extern template std::pair<Point, unsigned int> argmin<Strategy::INVERSE>(const Parameters & parameters, const Function & f , const Gradient & df);

/// Exponential decay function for determining step size

double expDecay(const Parameters & parameters, unsigned int k);

/// Inverse decay function for determining step size

double invDecay(const Parameters & parameters, unsigned int k);

/// Line search function using Armijo rule

double lineSearch(const Parameters & parameters, const Function & f , const Gradient & df, const Point & xk);

/// Armijo rule for line search

bool ArmijoRule(const Function & f , const Gradient & df, const Point & xk, double alpha, double sigma);

/// Euclidean distance (L2 norm) between two points (second origin by default)

double norm(const Point & x, const Point & y = {0,0});

#endif
