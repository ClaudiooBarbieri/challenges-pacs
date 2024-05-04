#include <iostream>
#include "matrix.hpp"

using ElementType = double;
using Order = algebra::StorageOrder;
template<typename T, Order O>
using Matrix = algebra::Matrix<T,O>;

int main(){

    Matrix<ElementType, Order::RowWise> rw_lns__131;
    Matrix<ElementType, Order::ColWise> cw_lns__131;
    std::vector<ElementType> res;
    std::string file = "../data/lns__131.mtx";

    rw_lns__131.readMatrix(file);
    cw_lns__131.readMatrix(file);

    std::vector<ElementType> v(rw_lns__131.getNCols(),1);

    res = rw_lns__131 * v;

    for( size_t i = 0 ; i < res.size(); ++i)
        std::cout << res[i] <<std::endl;
    


    /* TEST
    algebra::nonZeroElem<int> m_data = {
        {{0, 0}, 1},
        {{0, 2}, 5},
        {{3, 1}, 7},
        {{2, 1}, 6},
        {{3, 2}, 2},
        {{3, 3}, 3},
        {{4, 0}, 8}
    };

    algebra::Matrix<int, algebra::StorageOrder::RowWise> m1(5, 4, m_data);
    
    algebra::Matrix<int, algebra::StorageOrder::ColWise> m2(5, 4, m_data);  

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
    m2.compress();
    m1(0,1) = 3;
    std::cout << m1 << std::endl;

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

    

    std::cout << "***********************************" << std::endl;
    std::vector<int> v(4,1);
    std::vector<int> res1 = m1*v;
    for( const auto & x : res1)
        std::cout << x <<std::endl;
    std::cout << "*************************************" << std::endl;
    std::vector<int> res2 = m2*v;
    for( const auto & x : res2)
        std::cout << x <<std::endl;
    m1.compress();
    m2.compress();
    std::cout << "CCC***********************************" << std::endl;
    std::vector<int> resc1 = m1*v;
    for( const auto & x : resc1)
        std::cout << x <<std::endl;
    std::cout << "CCC***********************************" << std::endl;
    std::vector<int> resc2 = m2*v;
    for( const auto & x : resc2)
        std::cout << x <<std::endl;
    */
    return 0;
}