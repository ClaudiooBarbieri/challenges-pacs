#include <vector>
#include <cmath>
#include <iostream>
#include <utility>
#include "gradientMethod.hpp"
#include "parameters.hpp"
#include "functions.hpp"

typedef std::vector<double> Point;

template<Strategy S>
std::pair<Point,unsigned int> argmin(const Parameters & parameters, const Function & f , const Gradient & df){
    unsigned int iter = 0;
    double alpha = parameters.alpha0;
    Point xNew = parameters.x;
    Point x{0,0};
    do{
        x = xNew;
        xNew[0] = x[0]- alpha*df(x)[0];
        xNew[1] = x[1]- alpha*df(x)[1];
        ++iter;
        if constexpr (S == Strategy::EXPONENTIAL){
            alpha = expDecay(parameters,iter);
        }
        else if constexpr (S == Strategy::INVERSE){
            alpha = invDecay(parameters,iter);
        }
        else if constexpr (S == Strategy::ARMIJO){
            alpha = lineSearch(parameters,f,df,x);
        }
        else{
            std::cout << "INVALID STRATEGY!" << std::endl;
            return {x,iter};
        }
    }
    while(!(iter > parameters.maxIter || norm(xNew,x) < parameters.stepTol || std::fabs(f(xNew)-f(x))< parameters.resTol));
    return {xNew,iter};
}

template std::pair<Point,unsigned int> argmin<Strategy::ARMIJO>(const Parameters & parameters, const Function & f , const Gradient & df);
template std::pair<Point,unsigned int> argmin<Strategy::EXPONENTIAL>(const Parameters & parameters, const Function & f , const Gradient & df);
template std::pair<Point,unsigned int> argmin<Strategy::INVERSE>(const Parameters & parameters, const Function & f , const Gradient & df);

double expDecay(const Parameters & parameters, unsigned int k){
    return parameters.alpha0*exp(-parameters.mu*k);
}

double invDecay(const Parameters & parameters, unsigned int k){
    return parameters.alpha0/(1+parameters.mu*k);
}

double lineSearch(const Parameters & parameters, const Function & f , const Gradient & df, const Point & xk){
    double alphak = parameters.alpha0;
    while(!ArmijoRule(f,df,xk,alphak,parameters.sigma)){
        alphak/=2;
    }
    return alphak;
}

bool ArmijoRule(const Function & f , const Gradient & df, const Point & xk, double alpha, double sigma){
    std::cout << "n2 "<<  std::pow(norm(df(xk)),2) << std::endl;
    return f(xk)-f({xk[0]-alpha*df(xk)[0],xk[1]-alpha*df(xk)[1]})>=sigma*alpha*std::pow(norm(df(xk)),2);
}

double norm(const Point & x, const Point & y){
    return std::sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]));
}
