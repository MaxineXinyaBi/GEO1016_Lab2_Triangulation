/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
    double fx, double fy,     /// input: the focal lengths (same for both cameras)
    double cx, double cy,     /// input: the principal point (same for both cameras)
    double s,                 /// input: the skew factor (same for both cameras)
    const std::vector<Vector2D>& points_0,  /// input: 2D image points in the 1st image.
    const std::vector<Vector2D>& points_1,  /// input: 2D image points in the 2nd image.
    std::vector<Vector3D>& points_3d,       /// output: reconstructed 3D points
    Matrix33& R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
    Vector3D& t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       'triangulation()'.

    std::cout << "\nTODO: implement the 'triangulation()' function in the file 'Triangulation/triangulation_method.cpp'\n\n";

    std::cout << "[Liangliang]:\n"
        "\tSimilar to the first assignment, basic linear algebra data structures and functions are provided in\n"
        "\tthe following files:\n"
        "\t    - Triangulation/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
        "\t    - Triangulation/vector.h: manages vectors of arbitrary sizes and related functions.\n"
        "\t    - Triangulation/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
        "\tFor more details about these data structures and a complete list of related functions, please\n"
        "\trefer to the header files mentioned above.\n\n"
        "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
        "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
        "\tFor your final submission, adhere to the following guidelines:\n"
        "\t    - submit ONLY the 'Triangulation/triangulation_method.cpp' file.\n"
        "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
        "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
        0, 2.2, 3.3,
        0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
        0, 2.2, 3.3, 1,
        0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, { 0, 1, 2, 3, 4, 5, 6, 7, 8 }); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    if (points_0.size() < 8 || points_0.size() != points_1.size()) {
        std::cerr << "Invalid input: not enough points or point numbers don't match.\n";
        return false;
    }
    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;

    size_t num_points = points_0.size();
    // 1.Calculate the mean
    double mean_x0 = 0.0, mean_y0 = 0.0, mean_x1 = 0.0, mean_y1 = 0.0;
    for (size_t i = 0; i < num_points; ++i) {
        mean_x0 += points_0[i][0];
        mean_y0 += points_0[i][1];
        mean_x1 += points_1[i][0];
        mean_y1 += points_1[i][1];
    }
    mean_x0 /= num_points;
    mean_y0 /= num_points;
    mean_x1 /= num_points;
    mean_y1 /= num_points;

    // 2.Calculate the scale
    double scale0 = 0.0, scale1 = 0.0;
    for (size_t i = 0; i < num_points; ++i) {
        double dx0 = points_0[i][0] - mean_x0;
        double dy0 = points_0[i][1] - mean_y0;
        scale0 += sqrt(dx0 * dx0 + dy0 * dy0);

        double dx1 = points_1[i][0] - mean_x1;
        double dy1 = points_1[i][1] - mean_y1;
        scale1 += sqrt(dx1 * dx1 + dy1 * dy1);
    }
    scale0 = sqrt(2.0) / (scale0 / num_points);
    scale1 = sqrt(2.0) / (scale1 / num_points);

    // 3.Construct a normalized transformation matrix
    Matrix33 T0(
        scale0, 0, -scale0 * mean_x0,
        0, scale0, -scale0 * mean_y0,
        0, 0, 1
    );
    Matrix33 T1(
        scale1, 0, -scale1 * mean_x1,
        0, scale1, -scale1 * mean_y1,
        0, 0, 1
    );

    // 4.Construct the A matrix
    Matrix A(num_points, 9);
    for (size_t i = 0; i < num_points; ++i) {
        double x0 = scale0 * (points_0[i][0] - mean_x0);
        double y0 = scale0 * (points_0[i][1] - mean_y0);
        double x1 = scale1 * (points_1[i][0] - mean_x1);
        double y1 = scale1 * (points_1[i][1] - mean_y1);
        A.set_row(i, { x0 * x1, x0 * y1, x0, y0 * x1, y0 * y1, y0, x1, y1, 1 });
    }

    // 5.Solve F by SVD
    Matrix U, S, V;
    svd_decompose(A, U, S, V);

    Matrix33 F;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            F(i, j) = V(i, 8);

    Matrix Uf, Sf, Vf;
    svd_decompose(F, Uf, Sf, Vf);

    Sf(2, 2) = 0; // enforce rank 2
    F = Uf * Sf * Vf.transpose();

    F = T1.transpose() * F * T0; // denormalize

    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}