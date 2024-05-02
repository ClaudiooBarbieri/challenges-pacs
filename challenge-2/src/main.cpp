#include <iostream>
#include "matrix.hpp"

int main(){
    std::map<std::array<size_t,2>,double> trial;
    trial[{1,2}]=5.1;
    algebra::Matrix<double,algebra::StorageOrder::RowWise> m(3,3,trial);
    std::cout << m.isCompressed() << std::endl;
    return 0;
}