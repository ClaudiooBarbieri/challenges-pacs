#ifndef POINT_2D_HPP
#define POINT_2D_HPP

namespace challenge3{

    /**
    * @brief 2 dimensional point 
    * @note coordinates of double type
    */

    class Point2D{

    private:

        double x; ///< x coordinate
        double y; ///< y coordinate

    public:

        // constructor
        explicit Point2D(double x_, double y_) : x{x_},y{y_} {};

        // getters
        inline double getX() const { return x; }
        inline double getY() const { return y; }

    };
    
}

#endif