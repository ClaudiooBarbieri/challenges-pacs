#include <string>
#include <iostream>
#include <fstream> 
#include "parameters.hpp"
#include "json.hpp"

/// Read parameters from a JSON file.

Parameters readParameters(const std::string & parFileName){
    Parameters parameters; ///< Parameters object to hold the read parameters
    Parameters defaultParameters; ///< Default parameters used when values are not provided in the JSON file

    // Reading JSON file with data
    std::ifstream f(parFileName);
    nlohmann::json parFile = nlohmann::json::parse(f);      

    // Read options
    parameters.maxIter = parFile["option"].value("n_max_iter", defaultParameters.maxIter);
    parameters.resTol = parFile["option"].value("tol_res", defaultParameters.resTol);
    parameters.stepTol = parFile["option"].value("tol_step", defaultParameters.stepTol);
    parameters.alpha0 = parFile["option"].value("alpha0", defaultParameters.alpha0);
    parameters.mu = parFile["option"].value("mu", defaultParameters.mu);
    parameters.sigma = parFile["option"].value("sigma", defaultParameters.sigma);

    // Read starting point
    parameters.x = {parFile["point"].value("x1", 0.0), parFile["point"].value("x2", 0.0)};

    parameters.s = stringToStrategy(parFile["strategy"]); ///< Convert strategy string to enum

    return parameters; ///< Return the read parameters
}

/// Print parameters to standard output.

void Parameters::print() const{
    std::cout << std::endl;
    std::cout << " \t ### OPTION ###" << std::endl;
    std::cout << " Maximum number of iterations : " << maxIter << std::endl;
    std::cout << " Tolerance on the residual : " << resTol << std::endl;
    std::cout << " Tolerance on the steps : " << resTol << std::endl;
    std::cout << " alpha 0 : " << alpha0 << std::endl;
    std::cout << " mu : " << mu << std::endl;
    std::cout << " sigma : " << sigma << std::endl;
    std::cout << " \t ### STRATEGY ###" << std::endl;
    std::cout << strategyToString(s) << std::endl;
    std::cout << " \t ### STARTING POINT ###" << std::endl;
    std::cout << "x = (" << x[0] << "," << x[1]<< ")" << std::endl;
    std::cout << std::endl;
}

/// Convert Strategy enum to string.

std::string strategyToString(const Strategy & s){
    switch (s) {
        case Strategy::EXPONENTIAL:
            return {"exponential"};
        case Strategy::INVERSE:
            return {"inverse"};
        case Strategy::ARMIJO:
            return {"Armijo"};
    }
    return "not valid strategy!";
}

/// Convert string to Strategy enum.

Strategy stringToStrategy(const std::string & str){
    if (str == "exponential") {
        return Strategy::EXPONENTIAL;
    } else if (str == "inverse") {
        return Strategy::INVERSE;
    } else if (str == "Armijo") {
        return Strategy::ARMIJO;
    } else {
        return Strategy::NONE;
    }
}