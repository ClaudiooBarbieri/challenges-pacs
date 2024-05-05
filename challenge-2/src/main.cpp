#include <iostream>
#include "chrono.hpp"
#include "matrix.hpp"

using ElementType = double; ///< type of element stored in the matrix
using Order = algebra::StorageOrder; ///< order of storage in memory of the matrix
template<typename T, Order O>
using Matrix = algebra::Matrix<T,O>; 

/**
 * @brief prints the results of matrix-vector multiplication tests for each combination of StorageOrder and state of the matrix, compute also the %gain
 * @param urwTime UNcompressed RowWise time 
 * @param crwTime compressed RowWise time
 * @param ucwTime UNcompressed ColWise time
 * @param ccwTime compressed ColWise time
 * @param nIt number of itereation made in the time calculation
 */
void printResults(double urwTime, double crwTime, double ucwTime, double ccwTime, size_t nIt);

int main(){

    constexpr size_t TEST_ITERATIONS = 30;  ///< number of testing iteration

    Matrix<ElementType, Order::RowWise> rwMatrix; ///< matrix stored RowWise 
    Matrix<ElementType, Order::ColWise> cwMatrix; ///< matrix stored ColWise

    std::string file = "../data/lns__131.mtx"; ///< file of matrix stored in Matrix Market (MM) format

    //read matrices 
    rwMatrix.readMatrix(file);
    cwMatrix.readMatrix(file);

    std::vector<ElementType> testVector(rwMatrix.getNCols(),1); ///< 1s vector of the right size to test matrix-vector multiplication

    Timings::Chrono clock; ///< chrono clock to evaluate times
    double urwTime{0}; ///< UNcompressed RowWise total time
    double crwTime{0}; ///< compressed RowWise total time
    double ucwTime{0}; ///< UNcompressed ColWise total time
    double ccwTime{0}; ///< compressed ColWise total time
    
    // loop to evaluate elapsed time for the matrix-vector multiplication each combination of state and StorageOrder
    for(size_t t = 0 ; t < TEST_ITERATIONS ; ++t){
        clock.start(); rwMatrix * testVector; clock.stop();
        urwTime+=clock.wallTime();
        rwMatrix.compress();
        clock.start(); rwMatrix * testVector; clock.stop();
        crwTime+=clock.wallTime();
        rwMatrix.uncompress();
        clock.start(); cwMatrix * testVector; clock.stop();
        ucwTime+=clock.wallTime();
        cwMatrix.compress();
        clock.start(); cwMatrix * testVector; clock.stop();
        ccwTime+=clock.wallTime();
        cwMatrix.uncompress();
    }

    // print the result
    printResults(urwTime,crwTime,ucwTime,ccwTime,TEST_ITERATIONS);
    return 0;
}

void printResults(double urwTime, double crwTime, double ucwTime, double ccwTime, size_t nIt) {
    std::cout << "Mean elapsed time of matrix-vector product:" << std::endl;
    std::cout << "UNcompressed RowWise: " << urwTime/nIt << " microseconds" << std::endl;
    std::cout << "Compressed RowWise: " << crwTime/nIt << " microseconds" << std::endl;
    double rwGain = ((urwTime - crwTime) / urwTime) * 100;
    std::cout << "RowWise Compression Gain: " << rwGain << "%" << std::endl;
    std::cout << "UNcompressed ColWise: " << ucwTime/nIt << " microseconds" << std::endl;
    std::cout << "Compressed ColWise: " << ccwTime/nIt << " microseconds" << std::endl;
    double cwGain = ((ucwTime - ccwTime) / ucwTime) * 100;
    std::cout << "ColWise Compression Gain: " << cwGain << "%" << std::endl; 
}