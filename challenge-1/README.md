# Gradient Descent Implementation

This repository contains an implementation of the gradient descent optimization algorithm.

## Overview

Gradient descent is an iterative optimization algorithm used to minimize some function by iteratively moving in the direction of steepest descent as defined by the negative of the gradient. It's commonly used in machine learning for optimizing parameters of a model.

## Implementation Details

The gradient descent algorithm implemented in this repository is designed to minimize a given objective function. The objective function and its gradient are provided in a C++ file named `functions.cpp`, and the initial parameters are specified in a JSON file named `parameters.json`.

## How to Use

To use the gradient descent implementation in your project, follow these steps:

1. Clone this repository to your local machine:

    ```bash
    git clone https://github.com/your_username/gradient-descent.git
    ```

2. Modify the `parameters.json` file to specify the optimization options and the starting point:

    ```json
    {
        "option": {
            "n_max_iter": 100,
            "tol_res": 1e-6,
            "tol_step": 1e-6,
            "alpha0": 0.25,
            "mu": 0.2,
            "sigma": 0.3
        },
        "point": {
            "x1": 0,
            "x2": 0
        },
        "strategy": "inverse"
    }
    ```

    - `option`: Specifies the optimization options including the maximum number of iterations (`n_max_iter`), tolerance for the residual (`tol_res`), tolerance for the step (`tol_step`), initial step size (`alpha0`), damping parameter for Armijo strategy (`mu`), and scale parameter for Armijo strategy (`sigma`).
    - `point`: Specifies the starting point for optimization.
    - `strategy`: Specifies the optimization strategy which could be "exponential", "inverse", or "Armijo".

3. If you want to change the objective function, edit the `functions.cpp` file and modify the `objective_function` and `gradient_function` functions.

4. Compile the C++ code:

    ```bash
    g++ -o gradient_descent functions.cpp -std=c++11
    ```

5. Run the compiled executable:

    ```bash
    ./gradient_descent
    ```

6. The optimized parameters will be printed to the console.

## Contributing

Contributions are welcome! If you'd like to contribute to this project, please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
