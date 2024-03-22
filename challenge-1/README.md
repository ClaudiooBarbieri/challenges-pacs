# Gradient Descent Implementation

This repository contains an implementation of the gradient descent optimization algorithm given as challenge for the PACS course 2023/24.

## Overview

Gradient descent is an iterative optimization algorithm used to minimize some function by iteratively moving in the direction of steepest descent as defined by the negative of the gradient. More detail about the requests in the doc folder

## Implementation Details

The gradient descent algorithm implemented in this repository is designed to minimize a given objective function. The objective function and its gradient are provided in a C++ file named `functions.cpp`, and the initial parameters are specified in a JSON file named `parameters.json` and stored in `Struct Parameters` defined in `parameters.hpp` . Additionally, the optimization algorithm according to the selected strategy is implemented in `gradientMethod.cpp`.

## How to Use

To use the gradient descent implementation in your project, follow these steps:

1. Clone this repository to your local machine:

    ```bash
    git clone git@github.com:ClaudiooBarbieri/challenges-pacs.git
    ```

2. Modify the `parameters.json` file to specify the optimization options and the starting point (otherwise defaults are provided):

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

    - `option`: Specifies the optimization options including the maximum number of iterations (`n_max_iter`), tolerance for the residual (`tol_res`), tolerance for the step (`tol_step`), initial step size (`alpha0`), damping parameter (`mu`), and scale parameter for Armijo strategy (`sigma`).
    - `point`: Specifies the starting point for optimization.
    - `strategy`: Specifies the optimization strategy which could be "exponential", "inverse", or "Armijo".

3. If you want to change the objective function, edit the `functions.cpp` file and modify `Function` and `Gradient` functors.

4. Set the `EXAMPLES_INCLUDE` and `EXAMPLES_LIB` variables in the Makefile to point to the locations where you have the examples and libraries provided by the `pacs-examples` repository. (it is intended for people who are enrolled in the course, if not the only thing you need is the `json.h` header so set them in order to reach it)

6. Run `make` command in the `src` directory and execute the `main`:

    ```bash
    cd src
    make
    ./main
    ```
