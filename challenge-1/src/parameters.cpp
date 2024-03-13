#include "parameters.hpp"
#include <string>
#include <fstream> 
#include "json.hpp"

using json = nlohmann::json;

parameters read_parameters(const std::string & file_name){
    //reading json file with data
    std::ifstream f(file_name);
    json data = json::parse(f);
    const unsigned int n_max_iter = data["option"].value("n_max_iter", 0.0);
    const double tol_res = data["option"].value("tol_res", 0.0);
    const double tol_step = data["option"].value("tol_step", 0.0);
    const double alpha0 = data["option"].value("alpha0", 0.0);
    const double mu = data["option"].value("mu", 0.0);
    const double sigma = data["option"].value("sigma", 0.0);

    return {n_max_iter,tol_res,tol_step,alpha0,mu,sigma};
}