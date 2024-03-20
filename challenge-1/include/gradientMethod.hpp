#ifndef GRADIENT_METHOD_HPP
#define GRADIENT_METHOD_HPP

#include <vector>
#include <utility>
#include "parameters.hpp"
#include "functions.hpp"

typedef std::vector<double> Point;

template<Strategy S>
std::pair<Point,unsigned int> argmin(const Parameters & parameters, const Function & f , const Gradient & df);

double expDecay(const Parameters & parameters, unsigned int k);

double invDecay(const Parameters & parameters, unsigned int k);

double lineSearch(const Parameters & parameters, const Function & f , const Gradient & df, const Point & xk);

bool ArmijoRule(const Function & f , const Gradient & df, const Point & xk, double alpha, double sigma);

double norm(const Point & x, const Point & y = {0,0});

#endif