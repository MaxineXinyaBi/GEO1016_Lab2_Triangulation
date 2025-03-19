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
#include <iostream>
#include <cmath>
#include <vector>

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

//    ...（原有注释略）...

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
        scale0 += sqrt(pow(points_0[i][0] - mean_x0, 2) + pow(points_0[i][1] - mean_y0, 2));
        scale1 += sqrt(pow(points_1[i][0] - mean_x1, 2) + pow(points_1[i][1] - mean_y1, 2));
    }
    scale0 = (scale0 > 0) ? sqrt(2.0) / (scale0 / num_points) : 1.0;
    scale1 = (scale1 > 0) ? sqrt(2.0) / (scale1 / num_points) : 1.0;

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

    // Check the dimensions of the V matrix
    if (V.cols() < 9) {
        std::cerr << "Error: V matrix does not have 9 columns! V.cols() = " << V.cols() << std::endl;
        return false;
    }

    // Get F, take the last column of V
    Matrix33 F;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            F(i, j) = V(i, V.cols() - 1);  // Ensure that no out-of-bounds access is made

    Matrix Uf, Sf, Vf;
    svd_decompose(F, Uf, Sf, Vf);

    Sf(2, 2) = 0; // enforce rank 2
    F = Uf * Sf * Vf.transpose();

    F = T1.transpose() * F * T0; // denormalize

    //      - compute the essential matrix E;
    //      - recover rotation R and t.
    // 1。 construct intrinsic matrix K
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);

    // 2. calculate E
    Matrix33 E = K.transpose() * F * K;
    // 3. SVD decompose of E
    Matrix33 EU, Esigma, EV;
    svd_decompose(E, EU, Esigma, EV);
    // 4. make sure the E has two identical sigular values and the 3rd is 0
    double new_sigma = (Esigma(0, 0) + Esigma(1, 1)) / 2.0;
    Matrix33 Sigma(new_sigma, 0, 0,
                   0, new_sigma, 0,
                   0, 0, 0);
    E = EU * Sigma * EV.transpose();
    // 4. svd decompose again
    svd_decompose(E, EU, Esigma, EV);
    // 5. rotation matrix
    Matrix33 W(0, -1, 0,
               1, 0, 0,
               0, 0, 1);
    Matrix33 r1 = EU * W * EV.transpose();
    Matrix33 r2 = EU * W.transpose() * EV.transpose();
    double determinant_r1 = determinant(r1);
    double determinant_r2 = determinant(r2);
    Matrix33 R1 = determinant_r1 * r1;
    Matrix33 R2 = determinant_r2 * r2;
    // 6.t is the last column of U
    Vector3D t1 = EU.get_column(EU.cols() - 1);
    Vector3D t2 = -EU.get_column(EU.cols() - 1);
    // there are 4 possible pairs (R1 t1), (R2 t2), (R1 t2), (R2 t1)

    // ============== 新增：获取投影矩阵并进行三角化 + 确定正确(R, t) ============= //

    // 辅助函数：根据K, R, t 构造投影矩阵 P = K * [R|t] (3×4)
    auto getProjectionMatrix = [&](const Matrix33& Kmat, const Matrix33& Rmat, const Vector3D& tvec) -> Matrix34
    {
        // 注意：这里必须预先分配一个 3x4 的矩阵
        Matrix34 RT(3, 4, 0.0);
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                RT(row, col) = Rmat(row, col);
            }
            RT(row, 3) = tvec[row];
        }
        Matrix34 P = Kmat * RT; // Kmat: 3x3, RT: 3x4, 结果为3x4
        std::cout << "Projection matrix computed. Example element P(0,0): " << P(0,0) << std::endl;
        return P;
    };

    // 辅助函数：使用DLT线性方法三角化一对2D点
    auto triangulateDLT = [&](const Matrix34& P1, const Matrix34& P2,
                              const Vector2D& pt1, const Vector2D& pt2) -> Vector3D
    {
        Matrix A_(4, 4, 0.0);
        // 第1行: pt1.y*P1(2,:) - P1(1,:)
        for (int c = 0; c < 4; c++){
            A_(0, c) = pt1[1] * P1(2, c) - P1(1, c);
        }
        // 第2行: P1(0,:) - pt1.x*P1(2,:)
        for (int c = 0; c < 4; c++){
            A_(1, c) = P1(0, c) - pt1[0] * P1(2, c);
        }
        // 第3行: pt2.y*P2(2,:) - P2(1,:)
        for (int c = 0; c < 4; c++){
            A_(2, c) = pt2[1] * P2(2, c) - P2(1, c);
        }
        // 第4行: P2(0,:) - pt2.x*P2(2,:)
        for (int c = 0; c < 4; c++){
            A_(3, c) = P2(0, c) - pt2[0] * P2(2, c);
        }
        // 输出A_的第一行以调试
        std::cout << "Triangulation A_(0,:) = ";
        for (int c = 0; c < 4; c++){
            std::cout << A_(0, c) << " ";
        }
        std::cout << std::endl;

        Matrix UU, SS, VV;
        svd_decompose(A_, UU, SS, VV);
        double x = VV(0, 3), y = VV(1, 3), z = VV(2, 3), w = VV(3, 3);
        if (fabs(w) < 1e-12) w = 1e-12;
        return Vector3D(x/w, y/w, z/w);
    };

    // 手动构造一个 3×3 单位矩阵
    Matrix33 I;
    I(0,0) = 1.0; I(0,1) = 0.0; I(0,2) = 0.0;
    I(1,0) = 0.0; I(1,1) = 1.0; I(1,2) = 0.0;
    I(2,0) = 0.0; I(2,1) = 0.0; I(2,2) = 1.0;

    // 定义第1个摄像机的投影矩阵 (假设C1与世界坐标系重合, R=I, t=0)
    Matrix34 P1_proj = getProjectionMatrix(K, I, Vector3D(0.0, 0.0, 0.0));

    // 收集4个候选 (R, t)
    std::vector<std::pair<Matrix33, Vector3D>> candidates = {
        { R1, t1 }, { R1, t2 }, { R2, t1 }, { R2, t2 }
    };

    int best_in_front_count = -1;
    Matrix33 best_R;
    Vector3D best_t;
    std::vector<Vector3D> best_3d_points;

    // 对每个候选解进行测试，统计三角化后在两摄像机中均有正深度的点数
    for (auto& cnd : candidates) {
        Matrix33 Rc = cnd.first;
        Vector3D tc = cnd.second;

        std::cout << "Testing candidate with R first element: " << Rc(0,0)
                  << ", t: (" << tc[0] << ", " << tc[1] << ", " << tc[2] << ")" << std::endl;

        // 计算第2个摄像机的投影矩阵
        Matrix34 P2_proj = getProjectionMatrix(K, Rc, tc);

        int count_in_front = 0;
        std::vector<Vector3D> tmp_points;
        tmp_points.reserve(num_points);

        for (size_t i = 0; i < num_points; i++) {
            Vector3D X = triangulateDLT(P1_proj, P2_proj, points_0[i], points_1[i]);
            tmp_points.push_back(X);

            // 对于第1个摄像机，世界坐标系即为其坐标系，检查Z > 0
            bool front_c1 = (X[2] > 0.0);
            // 对于第2个摄像机，将 X 转换到 C2 坐标系: X2 = Rc * X + tc
            Vector3D X2 = Rc * X + tc;
            bool front_c2 = (X2[2] > 0.0);

            if (front_c1 && front_c2)
                count_in_front++;
        }

        std::cout << "Candidate valid point count: " << count_in_front << " / " << num_points << std::endl;

        if (count_in_front > best_in_front_count) {
            best_in_front_count = count_in_front;
            best_R = Rc;
            best_t = tc;
            best_3d_points = tmp_points;
        }
    }

    if (best_in_front_count <= 0) {
        std::cerr << "Failed to find a valid relative pose with positive depths.\n";
        return false;
    }

    std::cout << "Best candidate has " << best_in_front_count << " valid points." << std::endl;

    // 更新输出 R, t 和 points_3d
    R = best_R;
    t = best_t;
    points_3d = best_3d_points;

    return !points_3d.empty();
}
