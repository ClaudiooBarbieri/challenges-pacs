# Sparse Matrix Class Template

This is a C++ class template for handling sparse matrices. It provides functionalities for matrix-vector multiplication, matrix compression, decompression, reading matrices from Matrix Market format files, and more.

## Features

- **Matrix creation**: Create sparse matrices either by providing non-zero elements directly or by reading from Matrix Market format files.
- **Compression**: Compress matrices in CSR CSC format to optimize memory usage
- **Decompression**: Decompress matrices to access individual elements efficiently.
- **Matrix-vector multiplication**: Multiply matrices with vectors efficiently.
- **Storage order**: Choose between row-wise or column-wise storage order for matrices.

## Example of Usage

To use this library, include the `matrix.hpp` header file in your C++ project and instantiate the `Matrix` class with the desired element type (`ElementType`) and storage order (`Order`). Here's an example:

```cpp
#include "matrix.hpp"

using ElementType = double; 
using Order = algebra::StorageOrder;
template<typename T, Order O>
using Matrix = algebra::Matrix<T,O>; 

int main() {
    // Create a row-wise matrix
    Matrix<ElementType, Order::RowWise matrix(3, 3, {
        {{0, 1}, 3.0},
        {{1, 0}, 2.0},
        {{2, 2}, 1.0}
    });

    // Compress the matrix
    matrix.compress();

    // Perform matrix-vector multiplication
    std::vector<ElementType> vector = {1.0, 2.0, 3.0};
    std::vector<ElementType> result = matrix * vector;

    // Output the result
    for (auto& val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Perform other operation

    return 0;
}
```

## Testing

The main function is made for testing routine to evaluate the performance of matrix-vector multiplication for different combinations of storage order and matrix state (compressed or uncompressed). Here's a brief description of the testing process:

1. **Initialization**: Initialize sparse matrices by reading them from Matrix Market format files. Two matrices are created, one stored in row-wise order and the other in column-wise order.

2. **Matrix-Vector Multiplication**: Perform matrix-vector multiplication using both compressed and uncompressed matrices. Timing information is collected for each combination.

3. **Results**: Calculate and print the mean elapsed time of matrix-vector product for each combination of storage order and matrix state. Additionally, the compression gain percentage is calculated for both row-wise and column-wise storage orders.

This testing routine helps assess the efficiency and effectiveness of matrix compression in terms of computational performance.


