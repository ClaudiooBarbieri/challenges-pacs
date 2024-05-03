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

    std::cout << m1 << std::endl;
    //m2.compress();
    m1(0,1) = 3;
    std::cout << m1 << std::endl;

    
    /*
    std::cout << "-------" << std::endl;
    std::cout << m1(0,0) << std::endl;
    std::cout << m1(0,3) << std::endl;
    std::cout << m1(2,3) << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << m2(0,0) << std::endl;
    std::cout << m2(0,3) << std::endl;
    std::cout << m2(3,1) << std::endl;
    std::cout << "-------" << std::endl;
*/
    algebra::Matrix<double, algebra::StorageOrder::RowWise> lns__131;

    std::string file = "../data/lns__131.mtx";

    lns__131.readMatrix(file);
    std::cout << "Matrix read:" << std::endl;
    //std::cout << lns__131 << std::endl;
    return 0;
}