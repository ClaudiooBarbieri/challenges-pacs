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

    ///  enumerator for the order used to store the matrix 
    enum class StorageOrder {RowWise, ColWise};

    /// comparator to manage the storage in the correct order
    template <StorageOrder Order>
    struct Comparator {};

    /**
     * @brief comparator for row wise storage
     * @note lexicographic wrt row-column
    */
    template <>
    struct Comparator<StorageOrder::RowWise> {
        bool operator()(const std::array<size_t, 2>& lhs, const std::array<size_t, 2>& rhs) const {
            return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]); //< first compare row index then column
        }
    };

    /**
     * @brief comparator for column wise storage
     * @note lexicographic wrt column-row
    */ 
    template <>
    struct Comparator<StorageOrder::ColWise> {
        bool operator()(const std::array<size_t, 2>& lhs, const std::array<size_t, 2>& rhs) const {
            return lhs[1] < rhs[1] || (lhs[1] == rhs[1] && lhs[0] < rhs[0]); //< first compare column index then row
        }
    };

    /// store indexes and values of non zero element of sparse matrix 
    template<typename T>
    using nonZeroElem = std::map<std::array<size_t,2>,T>;

    /**
     * @brief struct containing the vectors used for compressed storage
    */
    template<typename T>
    struct CompressedData{
        std::vector<T> nzElem; ///< values of the non zero elements
        std::vector<size_t> elemIdx; ///< column(CSR) row(CSC) indexes of the non zero elements
        std::vector<size_t> beginIdx; ///< indexes of start of the row(CSR) column(CSC) in previous vectors
    };

    /// type of storage for uncompressed matrix
    template<typename T, StorageOrder Order>
    using UncompressedData = std::map<std::array<size_t,2>,T,Comparator<Order>>;

    /**
     * @brief dynamic class template for sparse matrix 
     * @tparam T type of the elements stored in the matrix
     * @tparam Order storage order of the matrix 
    */
    template<typename T, StorageOrder Order>
    class Matrix {
    private:
        size_t rows{0}; ///< number of rows
        size_t cols{0}; ///< number of columns
        bool compressed{false}; ///< flag of the compressed state
        UncompressedData<T,Order> uData; ///< storage for the uncompressed state
        CompressedData<T> cData; ///< storage for compressed state

        /// check if the index is coherent with matrix sizes
        bool indexInBound(size_t i,size_t j) const {
            return i<rows && j<cols;
        }

        /// get the index in the vector of non zero element in position (i,j)
        int getIndex(size_t i, size_t j) const {
            if constexpr(Order == StorageOrder::RowWise){
                for(size_t r = cData.beginIdx[i] ; r < cData.beginIdx[i+1]; r++){ ///< loop over element of row i
                    if(cData.elemIdx[r] == j) ///< check if the column index correspond to the given one j
                        return r; ///< return index of non zero value in non zero element vector
                    }
            }
            else if constexpr(Order == StorageOrder::ColWise){
                for(size_t r = cData.beginIdx[j] ; r < cData.beginIdx[j+1]; r++){ ///< loop over element of column j
                    if(cData.elemIdx[r] == i) ///< check if the column index correspond to the given one i
                        return r; ///< return index of non zero value in non zero element vector
                    }
                }
            std::cerr << "Zero element in given position" << std::endl;
            return -1;
        }

        /// get row k, indexes and values of non zero elements
        nonZeroElem<T> getRow(size_t k) const{
            if constexpr (Order == StorageOrder::RowWise){
                return {uData.lower_bound({k,0}),uData.lower_bound({k,cols})};
            }
            else{
                throw std::invalid_argument("Cannot extract row from ColWise storage");
            }
        }

        /// get column k, indexes and values of non zero elements
        nonZeroElem<T> getCol(size_t k) const{
            if constexpr (Order == StorageOrder::ColWise){
                return {uData.lower_bound({0,k}),uData.lower_bound({rows,k})};
            }
            else{
                throw std::invalid_argument("Cannot extract column from RowWise storage");
            }
        }

        /// print vectors used for compressed state format storage
        void printCompressedVector() const{
            std::cout << "NZ :    "; ///< value of non zero element
            for(size_t k = 0 ; k < cData.nzElem.size() ; ++k){
                std::cout << cData.nzElem[k] << "   "; 
            }
            std::cout << std::endl;
            std::cout << "Outer : "; ///< column(CSR) row(CSC) index of the element
            for(size_t k = 0 ; k < cData.elemIdx.size() ; ++k){
                std::cout << cData.elemIdx[k] << "   ";
            }
            std::cout << std::endl;
            std::cout << "Inner : "; ///< starting row(CSR) column(CSC) index
            for(size_t k = 0 ; k < cData.beginIdx.size() ; ++k){
                std::cout << cData.beginIdx[k] << "   ";
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    public:

        /// default constructor
        Matrix() = default;

        /**
         * @brief constructor in uncompressed state
         * @note assume coherent values for the number of rows and columns
         * @param nrows number of rows 
         * @param ncols number of columns
         * @param data indexes and values of the non zero elemenents stored in a nonZeroElem type variable
        */
        Matrix(size_t nRows, size_t nCols, const nonZeroElem<T> & data) : rows{nRows},cols{nCols},uData(data.begin(),data.end()),compressed{false} { };

        /// get number of rows
        size_t getNRows() const { return rows; }

        /// get number of columns
        size_t getNCols() const { return cols; }

        /// check the state of the matrix
        bool isCompressed() const { return compressed; }

        /// read matrix in Matrix Market (MM) format
        void readMatrix(const std::string & fileName){
            std::ifstream ifs(fileName);
            if(!ifs.is_open()){
                std::cerr << "Cannot open the file! " << fileName << std::endl;
                return;
            }
            std::string line;
            // skip comment and read dimension information
            while (std::getline(ifs, line)) {
                if (line[0] == '%') ///< skip comments
                    continue; 
                std::istringstream iss(line);
                // set dimensions
                size_t nRows, nCols, nNonZeros;
                iss >> nRows >> nCols >> nNonZeros; ///< read number of rows, columns and number of non zero elements
                rows = nRows;
                cols = nCols;
                break;
            }
            // read values and indexes of non zero element
            while (std::getline(ifs, line)) {
                std::istringstream iss(line);
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
                // reserve the space necessary for the vector of compressed state
                cData.nzElem.reserve(uData.size());
                cData.elemIdx.reserve(uData.size());
                cData.beginIdx.reserve(rows + 1);
                size_t inserted{0}; ///< keep track of the number of non zero elements inserted in the compressed vector
                cData.beginIdx.emplace_back(inserted); ///< first index is always zero
                for(size_t k = 0 ; k < rows; ++k){
                    nonZeroElem<T> row = this->getRow(k); ///< extract the non zero elements in row k
                    for(const auto & el : row){
                        cData.elemIdx.emplace_back(el.first[1]); ///< add column index of the value to vector of column index
                        cData.nzElem.emplace_back(el.second); ///< add the value to vector of non zero elements
                    }
                    inserted+= row.size(); ///< update the number of inserted values
                    cData.beginIdx.emplace_back(inserted); ///< add starting index of the row k+1
                }
            }
            else if constexpr (Order == StorageOrder::ColWise){
                // reserve the space necessary for the vector of compressed state
                cData.nzElem.reserve(uData.size());
                cData.elemIdx.reserve(uData.size());
                cData.beginIdx.reserve(cols + 1);
                size_t inserted{0}; ///< keep track of the number of non zero elements inserted in the compressed vector
                cData.beginIdx.emplace_back(inserted); ///< first index is always zero
                for(size_t k = 0 ; k < cols; ++k){
                    nonZeroElem<T> col = this->getCol(k); ///< extract the non zero elements in col k
                    for(const auto & el : col){
                        cData.elemIdx.emplace_back(el.first[0]); ///< add row index of the value to vector of row index
                        cData.nzElem.emplace_back(el.second); ///< add the value to vector of non zero elements
                    }
                    inserted+= col.size(); ///< update the number of inserted values
                    cData.beginIdx.emplace_back(inserted); ///< add starting index of the column k+1
                }
            }
            uData.clear(); ///< clear memory used in uncompressed state
            compressed=true; ///< set flag of compressed state
        }

        /**
         * @brief uncompress the matrix in the corresponding form of its StorageOrder
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
                    for(size_t j = 0 ; j < cData.beginIdx[i+1]-cData.beginIdx[i];++j){ ///<loop over number of elements of row k has to be inserted
                        uData.insert({{i,cData.elemIdx[k]},cData.nzElem[k]}); ///< row index from outer loop,column index and value at index k
                        ++k; ///< got to next element to add
                    }
                }
            }
            else if constexpr (Order == StorageOrder::ColWise){
                size_t k = 0; ///< index of the non zero element value to insert
                for(size_t j = 0 ; j < cols ; ++j){ ///< loop over columns
                    for(size_t i = 0 ; i < cData.beginIdx[j+1]-cData.beginIdx[j];++i){ ///<loop over number of elements of column k has to be inserted
                        uData.insert({{cData.elemIdx[k],j},cData.nzElem[k]}); ///< column index from outer loop,row index and value at index k
                        ++k; ///< got to next element to add
                    }
                }
            }
            // clear memory used in compressed state 
            cData.nzElem.clear();
            cData.elemIdx.clear();
            cData.beginIdx.clear();
            compressed=false; ///< unset flag of compressed state 
        }

        /**
         * @brief add or change value of element at index i,j
         * @note if compressed can only change an existing non zero element
        */
        T& operator()(size_t i, size_t j){
            if(!this->indexInBound(i,j)){
                throw std::out_of_range("Indices are out of bounds"); ///< not valid index
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
         * @brief examine the value at index i,j
         * @note if out of bound throw an exception
        */
        const T operator()(size_t i, size_t j) const{
            if(this->indexInBound(i,j)){
                if(this->isCompressed()){
                    int idx = this->getIndex(i,j);
                    return idx==-1 ? 0 : cData.nzElem[idx]; ///< return value at index i,j, -1 means zero valued element
                }
                else{
                    const auto it = uData.find({i,j}); ///< iterator to element of index (i,j) if exists
                    return it!=uData.end() ? it->second : 0; ///< end iterator indicates a zero value (not stored)
                }
            }
            else{
                throw std::out_of_range("Indices are out of bounds"); ///< not valid index
            }
        }

        /**
         * @brief resize the dimension of the matrix 
         * @note only possible to enlarge it
        */
        void resize(size_t nRows, size_t nCols){
            if(this->isCompressed()){
                std::cerr << "Cannot resize while in compressed state!" << std::endl;
                return;
            }
            if(nRows>= rows && nCols>=cols){
                rows = nRows; ///< update row size
                cols = nCols; ///< update column size
            }
            else{
                std::cerr << "Cannot shrink the matrix!" << std::endl;
            }
        }

        /// matrix vector multiplication specification for rowwise
        template<typename U>
        friend std::vector<U> operator*(const Matrix<U,StorageOrder::RowWise> & matrix, const std::vector<U> & v) ;

        /// matrix vector multiplication specification for colwise
        template<typename U>
        friend std::vector<U> operator*(const Matrix<U,StorageOrder::ColWise> & matrix, const std::vector<U> & v) ;

        /// print the matrix according to its StorageOrder
        friend std::ostream& operator<<(std::ostream& os, const Matrix<T,Order> & matrix) {
            if(!matrix.isCompressed()){
                for (const auto& d : matrix.uData) {
                    os << "( " << d.first[0] << " : " << d.first[1] << " )" << " = " << d.second << std::endl; ///< print (row,col) index and value in the stored order
                }
            }
            else{
                std::cout << "Print Compressed" << std::endl;
                matrix.printCompressedVector(); ///< print the vector used for compressed state
            }
            return os;
        }
    };

    /**
     * @brief matrix-vector multiplication 
     * @note RowWise partial specialization
    */
    template<typename T>
    std::vector<T> operator*(const Matrix<T,StorageOrder::RowWise> & matrix, const std::vector<T> & v) {
        size_t nCols = matrix.getNCols();
        if(nCols != v.size())
            throw std::invalid_argument("No matching dimensions between matrix and vector"); ///< if number of matrix column different from the vector(column) size
        size_t nRows = matrix.getNRows(); ///< get resulting vector size
        std::vector<T> res(nRows,0); ///< set zero vector of the right dimension
        if(matrix.isCompressed()){
            for(size_t i = 0 ; i < matrix.getNRows(); ++i){
                for(size_t k = matrix.cData.beginIdx[i] ; k < matrix.cData.beginIdx[i+1]; ++k){ ///< loop over non zero element of row i
                    res[i] += matrix.cData.nzElem[k]*v[matrix.cData.elemIdx[k]]; ///< kth non zero element multiplied with the element of its same column index in the vector and add to row in the resulting vector
                }
            }
        }
        else{
            for(size_t i = 0 ; i < nRows ; ++i){
                nonZeroElem<T> row = matrix.getRow(i); ///< extract row i
                for(const auto & el : row){
                    res[i] += el.second*v[el.first[1]]; ///< product of non zero element with its corresponding value in the vector and sum to the resulting row
                }
            }
        }
        return res;
    }

    /**
     * @brief matrix-vector multiplication 
     * @note ColWise partial specialization
    */
    template<typename T>
    std::vector<T> operator*(const Matrix<T,StorageOrder::ColWise> & matrix, const std::vector<T> & v) {
        size_t nCols = matrix.getNCols();
        if(nCols != v.size())
            throw std::invalid_argument("No matching dimensions");///< if number of matrix column different from the vector(column) size
        size_t nRows = matrix.getNRows(); ///< get resulting vector size 
        std::vector<T> res(nRows,0); ///< set zero vector of the right dimension
        if(matrix.isCompressed()){
            for(size_t j = 0 ; j < matrix.getNCols(); ++j){
                for(size_t k = matrix.cData.beginIdx[j] ; k < matrix.cData.beginIdx[j+1]; ++k){ ///< loop over non zero element of column j
                    res[matrix.cData.elemIdx[k]] += matrix.cData.nzElem[k]*v[j]; ///< each non zero element of column j multiplied by jth row of the vector and added to its corresponding row index in the resulting vector
                }
            }
        }
        else{
            for(size_t j = 0 ; j < nCols ; ++j){
                nonZeroElem<T> col = matrix.getCol(j); ///< extract column j
                for(const auto & el : col){
                    res[el.first[0]] += el.second*v[el.first[1]]; ///< product of non zero element with its corresponding value in the vector and sum to the resulting row
                }
            }
        }
        return res;
    }

}

#endif