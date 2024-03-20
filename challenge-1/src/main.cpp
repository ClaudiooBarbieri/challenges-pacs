#include <iostream>
#include <utility>
#include "parameters.hpp"
#include "functions.hpp"
#include "gradientMethod.hpp"

/// Prints the result of optimization.

void printResult(const std::pair<Point, unsigned int> & result){
    std::cout <<  "argmin = (" << result.first[0] << " , " << result.first[1] << ")" <<std::endl;
    std::cout << "Reached in " << result.second << " iteration" << std::endl;
}

int main(){
    // Read parameters from JSON file
    Parameters parameters = readParameters("../data/parameters.json");
    // Print parameters
    parameters.print();
    // Create function and gradient objects
    Function f;
    Gradient df;
    // Declare result
    std::pair<Point,unsigned int> result;
    // Perform optimization based on selected strategy
    switch (parameters.s) {
        case Strategy::EXPONENTIAL:
            result = argmin<Strategy::EXPONENTIAL>(parameters, f, df);
            break;
        case Strategy::INVERSE:
            result = argmin<Strategy::INVERSE>(parameters, f, df);
            break;
        case Strategy::ARMIJO:
            result = argmin<Strategy::ARMIJO>(parameters, f, df);
            break;
    }
    // Print optimization result
    printResult(result);
}
