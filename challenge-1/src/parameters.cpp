#include <string>
#include <iostream>
#include <fstream> 
#include "parameters.hpp"
#include "json.hpp"


Parameters readParameters(const std::string & parFileName){
    //
    Parameters parameters;
    Parameters defaultParameters;

    //reading json file with data
    std::ifstream f(parFileName);
    nlohmann::json parFile = nlohmann::json::parse(f);      

    //option
    parameters.maxIter = parFile["option"].value("n_max_iter", defaultParameters.maxIter);
    parameters.resTol = parFile["option"].value("tol_res", defaultParameters.resTol);
    parameters.stepTol = parFile["option"].value("tol_step", defaultParameters.stepTol);
    parameters.alpha0 = parFile["option"].value("alpha0", defaultParameters.alpha0);
    parameters.mu = parFile["option"].value("mu", defaultParameters.mu);
    parameters.sigma = parFile["option"].value("sigma", defaultParameters.sigma);

    //point
    parameters.x = {parFile["point"].value("x1", 0.0),parFile["point"].value("x2", 0.0)};

    return parameters;
}

void Parameters::print() const{
    std::cout << std::endl;
    std::cout << " \t ### OPTION ###" << std::endl;
    std::cout << " Maximum number of iterations : " << maxIter << std::endl;
    std::cout << " Tolerance on the residual : " << resTol << std::endl;
    std::cout << " Tolerance on the steps : " << resTol << std::endl;
    std::cout << " alpha 0 : " << alpha0 << std::endl;
    std::cout << " mu = " << mu << std::endl;
    std::cout << " sigma = " << sigma << std::endl;
    std::cout << std::endl;
    std::cout << " \t ### STARTING POINT ###" << std::endl;
    std::cout << "x = (" << x[0] << "," << x[1]<< ")" << std::endl;
    std::cout << std::endl;
}