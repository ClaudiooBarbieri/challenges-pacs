#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include <vector>
#include <iostream>

typedef std::vector<double> point;

struct parameters{
    //options
    const unsigned int n_max_iter = 100;
    const double tol_res = 1e-3;
    const double tol_step = 1e-3;
    const double alpha0 = 0.5;
    const double mu = 0.5;
    const double sigma = 0.5;
    //functions

    //point
    point x;

    parameters(unsigned int it, double tol_r, double tol_s, double a, double m, double s) :
                n_max_iter(it),tol_res(tol_r),tol_step(tol_s),alpha0(a),mu(m),sigma(s){ std::cout << "strcut constructed" << std::endl;}
};


parameters read_parameters(const std::string & file_name);

#endif