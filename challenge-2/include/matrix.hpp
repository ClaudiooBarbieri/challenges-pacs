#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <limits>
#include <stdexcept>

namespace algebra{

    template<typename T>
    using nonZeroElem = std::map<std::array<size_t,2>,T>;

    /** 
     * @brief enumerator for the order used to store matrix values
    **/
    enum class StorageOrder {RowWise, ColWise};

    /// comparator
    template <StorageOrder Order>
    struct Comparator {};

    /// array comparison for storage by row , lexicographic wrt row,column 
    template <>
    struct Comparator<StorageOrder::RowWise> {
        bool operator()(const std::array<size_t, 2>& lhs, const std::array<size_t, 2>& rhs) const {
            return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]); //< first compare row index then column
        }
    };

    /// array comparison for storage by column, lexicographic wrt column,row 
    template <>
    struct Comparator<StorageOrder::ColWise> {
        bool operator()(const std::array<size_t, 2>& lhs, const std::array<size_t, 2>& rhs) const {
            return lhs[1] < rhs[1] || (lhs[1] == rhs[1] && lhs[0] < rhs[0]); //< first compare column index then row
        }
    };

    /**
     * @brief struct containing indexes and values of matrix in compressed state
    **/
    template<typename T>
    struct CompressedData{
        std::vector<size_t> innerIdx;
        std::vector<size_t> outerIdx;
        std::vector<T> nzElem;
    };

    /// type of storage for uncompressed matrix
    template<typename T, StorageOrder Order>
    using UncompressedData = std::map<std::array<size_t,2>,T,Comparator<Order>>;

    /**
     * @brief dynamic class template for sparse matrix 
     * @tparam T type of the elements stored in the matrix
     * @tparam Order storage order of the matrix 
    **/
    template<typename T, StorageOrder Order>
    class Matrix {
    private:
        /// matrix dimensions
        size_t rows{0};
        size_t cols{0};

        /// state of the matrix
        bool compressed{false};

        /// storage of non zero elements in uncompressed state
        UncompressedData<T,Order> uData;

        /// storage of non zero elements in compressed state
        CompressedData<T> cData;

        /// check if index is coherent with matrix sizes
        bool indexInBound(size_t i,size_t j) const {
            return i<rows && j<cols;
        }

        /**
         * @brief corresponding index in non zero value vector of compressed storage
         * @return index or -1 if correspond to zero value element
        **/
        int getIndex(size_t i, size_t j) const {
            if constexpr(Order == StorageOrder::RowWise){
                for(size_t r = cData.innerIdx[i] ; r < cData.innerIdx[i+1]; r++){ ///< loop over element of row i
                    if(cData.outerIdx[r] == j) ///< check if the column index correspond to the given one j
                        return r; ///< return index of non zero value in non zero element vector
                    }
            }
            else if constexpr(Order == StorageOrder::ColWise){
                for(size_t r = cData.innerIdx[j] ; r < cData.innerIdx[j+1]; r++){ ///< loop over element of col j
                    if(cData.outerIdx[r] == i) ///< check if the column index correspond to the given one i
                        return r; ///< return index of non zero value in non zero element vector
                    }
                }
            return -1;
        }

        /// print vectors of the compressed state format

        void printCompressedVector() const{
            std::cout << "NZ :    "; ///< value of non zero element
            for(size_t k = 0 ; k < cData.nzElem.size() ; ++k){
                std::cout << cData.nzElem[k] << "   "; 
            }
            std::cout << std::endl;
            
            std::cout << "Outer : "; ///< corresponding row/column index of the value
            for(size_t k = 0 ; k < cData.outerIdx.size() ; ++k){
                std::cout << cData.outerIdx[k] << "   ";
            }
            std::cout << std::endl;
            std::cout << "Inner : "; ///< starting corresponding index
            for(size_t k = 0 ; k < cData.innerIdx.size() ; ++k){
                std::cout << cData.innerIdx[k] << "   ";
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    public:

        /// deafualt constructor
        Matrix() = default;

        /**
         * @brief constructor of the matrix in uncompressed state
         * @note assume coherent values for the number of rows and columns
         * @param nrows number of rows of the matrix
         * @param ncols number of columns of the matrix
         * @param data map containing indexes and values of the non zero elemenents
        **/
        Matrix(size_t nRows, size_t nCols, const nonZeroElem<T> & data) : rows{nRows},cols{nCols},uData(data.begin(),data.end()),compressed{false} { };

        /// read matrix in Matrix Market (MM) format
        void readMatrix(const std::string & fileName){
            std::ifstream ifs(fileName);
            if(!ifs.is_open()){
                std::cerr << "Cannot open the file! " << fileName << std::endl;
                return;
            }
            std::string line;
            while (std::getline(ifs, line)) {
                if (line[0] == '%') ///< skip comments
                    continue; 
                std::istringstream iss(line);
                // set dimensions
                size_t nRows, nCols, nNonZeros;
                iss >> nRows >> nCols >> nNonZeros; ///< read number of rows, columns and non zero elements
                rows = nRows;
                cols = nCols;
                break;
            }
            while (std::getline(ifs, line)) {
                std::istringstream iss(line);
                // set values
                size_t i,j;
                T val;
                iss >> i >> j >> val; ///< read row,column index and corresponding value
                uData.insert({{i-1,j-1},val}); ///< indexing in MM are 1 based
            }
        }

        /**
         * @brief compress the matrix in the compressed form of its StorageOrder
         * @note empty the map storing the values
        */
        void compress(){
            if(this->isCompressed()){
                std::cout << "Already in compressed state!" << std::endl;
                return;
            }
            if constexpr(Order == StorageOrder::RowWise){
                // reserve the space necessary for the vector of compressed state representation
                cData.nzElem.reserve(uData.size());
                cData.outerIdx.reserve(uData.size());
                cData.innerIdx.reserve(rows + 1);
                size_t inserted{0}; ///< keep track of the number of non zero elements inserted in the compressed vector
                cData.innerIdx.emplace_back(inserted); ///< first index is always zero
                for(size_t k = 0 ; k < rows; ++k){
                    nonZeroElem<T> row(uData.lower_bound({k,0}),uData.lower_bound({k,cols})); ///< extract the non zero elements of row k
                    for(const auto & el : row){
                        cData.outerIdx.emplace_back(el.first[1]); ///< add column index of the value to vector of column index
                        cData.nzElem.emplace_back(el.second); ///< add the value to vector of non zero elements
                    }
                    inserted+= row.size(); ///< update the number of inserted values
                    cData.innerIdx.emplace_back(inserted); ///< add index to vector of rows index
                }
            }
            else if constexpr (Order == StorageOrder::ColWise){
                // reserve the space necessary for the vector of compressed state representation
                cData.nzElem.reserve(uData.size());
                cData.outerIdx.reserve(uData.size());
                cData.innerIdx.reserve(cols + 1);
                size_t inserted{0}; ///< keep track of the number of non zero elements inserted in the compressed vector
                cData.innerIdx.emplace_back(inserted); ///< first index is always zero
                for(size_t k = 0 ; k < cols; ++k){
                    nonZeroElem<T> col(uData.lower_bound({0,k}),uData.lower_bound({rows,k})); ///< extract the non zero elements of col k
                    for(const auto & el : col){
                        cData.outerIdx.emplace_back(el.first[0]); ///< add row index of the value to vector of row index
                        cData.nzElem.emplace_back(el.second); ///< add the value to vector of non zero elements
                    }
                    inserted+= col.size(); ///< update the number of inserted values
                    cData.innerIdx.emplace_back(inserted); ///< add index to vector of column index
                }
            }
            uData.clear();
            compressed=true;
        }

        /**
         * @brief uncompress the matrix in the uncompressed form of its StorageOrder
         * @note empty the vectors used by compressed form
        */
        void uncompress(){
            if(!this->isCompressed()){
                std::cout << "Already in uncompressed state!" << std::endl;
                return;
            }
            if constexpr(Order == StorageOrder::RowWise){
                size_t k = 0; ///< index of the non zero element value to insert
                for(size_t i = 0 ; i < rows ; ++i){ ///< loop over rows
                    for(size_t j = 0 ; j < cData.innerIdx[i+1]-cData.innerIdx[i];++j){ ///<loop over number of elements of row k has to be inserted
                        uData.insert({{i,cData.outerIdx[k]},cData.nzElem[k]}); ///< row index from outer loop,column index and value at index k of the corrisponding vector 
                        ++k;
                    }
                }
            }
            else if constexpr (Order == StorageOrder::ColWise){
                size_t k = 0; ///< index of the non zero element value to insert
                for(size_t j = 0 ; j < cols ; ++j){ ///< loop over columns
                    for(size_t i = 0 ; i < cData.innerIdx[j+1]-cData.innerIdx[j];++i){ ///<loop over number of elements of column k has to be inserted
                        uData.insert({{cData.outerIdx[k],j},cData.nzElem[k]}); ///< column index from outer loop,row index and value at index k of the corrisponding vector 
                        ++k;
                    }
                }
            }
            cData.nzElem.clear();
            cData.outerIdx.clear();
            cData.innerIdx.clear();
            compressed=false;
        }

        /**
         * @brief non const call operator add or change value of element at index i,j
         * @note if uncompressed can change and add, if compressed just change existing
        **/
        T& operator()(size_t i, size_t j){
            if(!this->indexInBound(i,j)){
                throw std::out_of_range("Indices are out of bounds"); ///< throw exception, not valid index
            }
            if(this->isCompressed()){
                int idx = this->getIndex(i,j); ///< get index of possible non zero value element
                if(idx==-1){
                    throw std::invalid_argument("Cannot modify non zero value in compressed state"); ///< cannot change a zero valued element
                }
                else{
                    return cData.nzElem[idx];
                }
            }   
            else{
                return uData[{i,j}];
            }
        }


        /**
         * @brief const call operator, give the value at index i,j
         * @note if out of bound throw an exception
        **/
        const T operator()(size_t i, size_t j) const{
            if(this->indexInBound(i,j)){
                if(this->isCompressed()){
                    int idx = this->getIndex(i,j);
                    return idx==-1 ? 0 : cData.nzElem[idx]; ///< return value at index i,j, -1 means zero valued element
                }
                else{
                    const auto it = uData.find({i,j}); ///< iterator to element of index (i,j) if exists
                    return it!=uData.end() ? it->second : 0; ///< if the iteretor is not the end return the corrisponding value of the non zero element otherwise zero
                }
            }
            else{
                throw std::out_of_range("Indices are out of bounds"); ///< throw exception, not valid index
            }
        }


        /// check the state of the matrix
        bool isCompressed() const { return compressed; }

        /// resize the dimension of the matrix (only possible to enlarge it)
        void resize(size_t nRows, size_t nCols){
            if(this->isCompressed()){
                std::cerr << "Cannot resize while in compressed state!" << std::endl;
                return;
            }
            if(nRows>= rows && nCols>=cols){
                rows = nRows;
                cols = nCols;
            }
            else{
                std::cerr << "Cannot shrink the matrix!" << std::endl;
            }
        }

        /// overload stream operator for the matrix
        friend std::ostream& operator<<(std::ostream& os, const Matrix & matrix) {
            if(!matrix.isCompressed()){
                for (const auto& d : matrix.uData) {
                    os << "( " << d.first[0] << " : " << d.first[1] << " )" << " = " << d.second << std::endl; ///< print (row,col) index and value in the stored order
                }
            }
            else{
                std::cout << "Print Compressed" << std::endl;
                matrix.printCompressedVector();
            }
            return os;
        }
    };
}

#endif