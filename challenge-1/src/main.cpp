#include <iostream>
#include <utility>
#include "parameters.hpp"
#include "functions.hpp"
#include "gradientMethod.hpp"

void printResult(const std::pair<Point, unsigned int> & result){
    std::cout <<  "argmin = (" << result.first[0] << " , " << result.first[1] << ")" <<std::endl;
    std::cout << "Reached in " << result.second << " iteration" << std::endl;
}

int main(){
    Parameters parameters = readParameters("../data/parameters.json");
    parameters.print();
    Function f;
    Gradient df;
    auto result = argmin<Strategy::INVERSE>(parameters,f,df);
    printResult(result);
    return 0;
}
