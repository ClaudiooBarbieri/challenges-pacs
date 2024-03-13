#include <iostream>
#include <fstream>
#include "parameters.hpp"
#include "json.hpp"


int main(){
    parameters p = read_parameters("../data/data.json");
    return 0;
}