#include <iostream>
#include "matrix.hpp"

int main(){
    algebra::nonZeroElem<int> m_data = {
        {{0, 0}, 1},
        {{0, 2}, 5},
        {{3, 1}, 7},
        {{2, 1}, 6},
        {{3, 2}, 2},
        {{3, 3}, 3}
    };

    // Instantiate Matrix with RowWise storage order
    algebra::Matrix<int, algebra::StorageOrder::RowWise> m1(4, 4, m_data);
    

    // Instantiate Matrix with ColWise storage order
    algebra::Matrix<int, algebra::StorageOrder::ColWise> m2(4, 4, m_data);

    /*algebra::nonZeroElem<int> row(m_data.lower_bound({0,0}),m_data.lower_bound({0,3}));
    for(size_t k=0;k<3;++k){
        algebra::nonZeroElem<int> row(m_data.lower_bound({k,0}),m_data.lower_bound({k,3}));
        for( auto v : row){
        std::cout << v.second <<  "-";
        }
        std::cout << row.size() << std::endl;
    }*/
    

    // Print matrices
    std::cout << "Matrix with RowWise storage order:" << std::endl;
    std::cout << m1 << std::endl;
    std::cout << "Compress" << std::endl;
    m1.compress();
    std::cout << m1 << std::endl;
    std::cout << "UNcompress" << std::endl;
    m1.uncompress();
    std::cout << m1 << std::endl;

    std::cout << "Matrix with ColWise storage order:" << std::endl;
    std::cout << m2 << std::endl;
    std::cout << "Compress" << std::endl;
    m2.compress();
    std::cout << m2 << std::endl;
    std::cout << "UNcompress" << std::endl;
    m2.uncompress();
    std::cout << m2 << std::endl;

    //algebra::Matrix<double, algebra::StorageOrder::RowWise> m;

    /*
    std::string file = "../data/lns__131.mtx";

    m.readMatrix(file);
    std::cout << "Matrix read:" << std::endl;
    std::cout << m << std::endl;
    */
    return 0;
}