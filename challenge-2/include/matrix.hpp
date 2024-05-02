#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <map>
#include <iostream>

namespace algebra{

    /** 
     * @brief Enum Class to specify the order storage of the matrix
    **/
    enum class StorageOrder {RowWise, ColWise};

    /**
     * @brief dynamic matrix class template
     * 
     * @tparam T type of the elements stored in the matrix
     * @tparam Order storage order of the matrix
    */
    template<typename T, StorageOrder Order>
    class Matrix {
    public:
        /**
         * @brief constructor of the matrix in uncompressed state
         * 
         * @param rows_ number of rows of the matrix
         * @param cols_ number of columns of the matrix
         * @param data_ map with indexes and values of non zero elements of the matrix
        */
        Matrix(size_t rows_, size_t cols_, std::map<std::array<size_t,2>,T> data_) : rows{rows_},cols{cols_},data{data_},compressed{false} {};

        // check the state of the matrix
        const bool isCompressed() const {return compressed;};

        // resize the dimension of the matrix leaveing it in uncompressed state
        void resize(size_t nrows, size_t ncols);

        // add or change values in uncompressed state, only change existing in compressed state
        T& operator()(size_t i, size_t j);

        // access in bound elements of the matrix
        const T& operator()(size_t i, size_t j) const;

        // method that converts the storage from uncompressed state to compressed state
        void compress();

        // method that brings back the matrix to the uncompressed state emptying the storage used for compressed format
        void uncompress();

        // matrix-vector product
        friend std::vector<T> operator*(const Matrix<T, Order>& m, const std::vector<T>& v);

        // reader of matrices written in matrix market format
        void read();

    private:
        size_t rows{0},cols{0}; ///< matrix dimensions
        bool compressed{false}; ///< boolean to describe the state of the matrix representation

        // uncompressed state variable
        std::map<std::array<size_t,2>,T> data; ///< non zero elements storage for compressed state, indexes and values
        
        // compressed state variable
        std::vector<size_t> innerIdx; ///< inner indexes vector for compressed state
        std::vector<size_t> outerIdx; ///< outer indexes vector for compressed state
        std::vector<T> values; ///< vector of values for compressed state
        
        // inline method to check the validity of the indexing
        inline bool validIndex(size_t i, size_t j) const {return i<rows && j<cols;};

        // empty variable unused after the change of state
        void clear();

    };

}

#endif