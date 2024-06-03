# A matrix–free parallel solver for the Laplace equation

All the details about the problem can be found in : [Challenge 3](../challenge-3/doc/Challenge23-24_3.pdf)

## Setup Instructions

#### Setup `PACS_ROOT` Variable in the [Makefile](../challenge-3/src/Makefile)

To configure the build environment, you need to specify the path directory to your `Examples` of [PACS](https://github.com/pacs-course/pacs-examples.git), you need Need `muParser` and `MPI` libraries. Also need of `json.hpp`.

```makefile
PACS_ROOT = your path
```
## General Possible Usage

#### Edit [JSON](../challenge-3/data/param.json) parameters accordingly to your need

- `minX`, `minY`: Minimum x and y coordinates of the mesh grid.
- `maxX`, `maxY`: Maximum x and y coordinates of the mesh grid.
- `nx`, `ny`: Number of divisions along the x and y axes, respectively.
- `hx`, `hy`: Step size along the x and y axes, respectively.
- `iteration`: Number of maximum iterations.
- `tolerance`: Tolerance level on solution changes for convergence.
- `f`: The main function to be used in the simulation.
- `fbound`: Function used to setup boundary condition.
- `exact`: Exact solution function for testing.

Can use the function defined in the `main` to read them and put in a struct

#### Create `Mesh2D` on which you want to use the solver
```cpp
// Creating a mesh with specified ranges and number of points
    challenge3::Mesh2D mesh = challenge3::Mesh2D::createWithPoints(minX, maxX, minY, maxY, nx, ny);

// Creating a mesh with specified ranges and spacing
    challenge3::Mesh2D mesh = challenge3::Mesh2D::createWithSpacing(minX, maxX, minY, maxY, hx, hy);

// Creating a mesh with specified vertices and number of points
    challenge3::Point2D topRight(1.0, 1.0), bottomLeft(0.0, 0.0);
    challenge3::Mesh2D mesh = challenge3::Mesh2D::createWithPoints(topRight, bottomLeft, nx, ny);

// Creating a mesh with specified vertices and spacing
    challenge3::Mesh2D mesh = challenge3::Mesh2D::createWithSpacing(topRight, bottomLeft, hx, hy);
```

#### Use the `JacobiSolver` class

```cpp
// Example
    JacobiSolver solver(mesh, "function", nMax, tolerance, "boundary condition function");
```

- `mesh_` (const Mesh2D&): The mesh on which to construct the solution.
- `function` (const std::string&): Expression of the function.
- `bCond` (const std::string&): Expression for the function used for boundary condition. Default value is "0".
- `nMax_` (unsigned int): Maximum number of iterations.
- `tol_` (double): Tolerance.

#### Generate VTK file for visualization of the solution

```cpp
// Example
    generateVTK(values, x0, y0, nx, ny, hx, hy, "extra");
```

- `values` (const &type ): Values in the form of a vector or vector of vector representing a matrix to be translated in VTK file.
- `x0` (double): X coordinate of the origin of the grid.
- `y0` (double): Y coordinate of the origin of the grid.
- `nx` (size_t): Number of points along the x direction.
- `ny` (size_t): Number of points along the y direction.
- `hx` (double): Spacing in the x direction.
- `hy` (double): Spacing in the y direction.
- `extra` (const std::string&): String appended to the output file name. Default value is an empty string.



### Focus on implemented `main`

1. Initializes MPI for parallel computing.

2. Reads parameters for the PDE solver from a JSON file.

3. Initializes a 2D mesh based on the parameters.

4. Measures the start time for performance evaluation.

5. Solves the PDE in parallel and produces VTK files for visualization.

6. Measures the end time for performance evaluation.

7. If the process rank is 0 (master process), compares the computed solution with an 
exact solution and prints error and elapsed time.

8. Synchronizes all processes.

9. Finalizes MPI.

## Testing 

In the folder [test](../challenge-3/test/) you can find a [script](../challenge-3/test/test.sh) that runs a small scalability test with 1,2 and 4 processes according to the parameters contained in the JSON parameter [file](../challenge-3/data/param.json). Ypu can try also with different values.

There [my result](../challenge-3/test/my_test.md) and [hardware info](../challenge-3/HWinfo/)

To run the test, go to `../challenge-3/test/` and run

```bash
./test.sh
```
##### note
not able to implement efficient OMP
Test with `JacobiSover` are made using `chrono.hpp` of PACS library, in order to use it need to link `-lpacs`. [code](../challenge-3/test/mainTest.cpp)



