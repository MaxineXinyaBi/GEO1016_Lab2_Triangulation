#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <iostream>
#include <cmath>
#include <vector>

using namespace easy3d;

// 辅助函数：打印矩阵
void printMatrix(const Matrix& M, const std::string& name) {
    std::cout << name << " (" << M.rows() << "x" << M.cols() << "):\n";
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            std::cout << M(i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// 辅助函数：打印 Vector3D
void printVector3D(const Vector3D& v, const std::string& name) {
    std::cout << name << " = (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
}

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
    std::cout << "\nStarting triangulation process...\n\n";

    if (points_0.size() < 8 || points_0.size() != points_1.size()) {
        std::cerr << "Invalid input: not enough points or point numbers don't match.\n";
        return false;
    }
    size_t num_points = points_0.size();

    // 1. 计算均值
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
    std::cout << "Mean for image 0: (" << mean_x0 << ", " << mean_y0 << ")\n";
    std::cout << "Mean for image 1: (" << mean_x1 << ", " << mean_y1 << ")\n\n";

    // 2. 计算尺度因子
    double scale0 = 0.0, scale1 = 0.0;
    for (size_t i = 0; i < num_points; ++i) {
        scale0 += sqrt(pow(points_0[i][0] - mean_x0, 2) + pow(points_0[i][1] - mean_y0, 2));
        scale1 += sqrt(pow(points_1[i][0] - mean_x1, 2) + pow(points_1[i][1] - mean_y1, 2));
    }
    scale0 = (scale0 > 0) ? sqrt(2.0) / (scale0 / num_points) : 1.0;
    scale1 = (scale1 > 0) ? sqrt(2.0) / (scale1 / num_points) : 1.0;
    std::cout << "Scale for image 0: " << scale0 << "\n";
    std::cout << "Scale for image 1: " << scale1 << "\n\n";

    // 3. 构造归一化变换矩阵 T0 和 T1
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
    printMatrix(T0, "T0");
    printMatrix(T1, "T1");

    // 4. 构造 A 矩阵
    Matrix A(num_points, 9);
    for (size_t i = 0; i < num_points; ++i) {
        double x0 = scale0 * (points_0[i][0] - mean_x0);
        double y0 = scale0 * (points_0[i][1] - mean_y0);
        double x1 = scale1 * (points_1[i][0] - mean_x1);
        double y1 = scale1 * (points_1[i][1] - mean_y1);
        A.set_row(i, { x0 * x1, x0 * y1, x0, y0 * x1, y0 * y1, y0, x1, y1, 1 });
    }
    printMatrix(A, "A");

    // 5. 对 A 进行 SVD 分解求解 F
    Matrix U, S, V;
    svd_decompose(A, U, S, V);
    printMatrix(U, "U from A");
    printMatrix(S, "S from A");
    printMatrix(V, "V from A");

    // 从 V 的最后一列构造 F
    Matrix33 F;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            F(i, j) = V(i, V.cols() - 1);
    printMatrix(F, "F (before rank enforcement)");

    // 对 F 进行 SVD 分解以强制其秩为2
    Matrix Uf, Sf, Vf;
    svd_decompose(F, Uf, Sf, Vf);
    Sf(2, 2) = 0; // 强制最低奇异值为0
    F = Uf * Sf * Vf.transpose();
    printMatrix(F, "F (after rank enforcement)");

    // 对 F 进行去归一化
    F = T1.transpose() * F * T0;
    printMatrix(F, "F (after denormalization)");

    // 6. 计算本质矩阵 E
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);
    std::cout << "Intrinsic matrix K:\n";
    printMatrix(K, "K");

    Matrix33 E = K.transpose() * F * K;
    printMatrix(E, "E (initial)");

    // 7. 对 E 进行 SVD 分解
    Matrix33 EU, Esigma, EV;
    svd_decompose(E, EU, Esigma, EV);
    printMatrix(EU, "EU from E");
    printMatrix(Esigma, "Esigma from E");
    printMatrix(EV, "EV from E");

    // 8. 调整 Esigma 使其符合本质矩阵的性质（前两个相等，第三个为零）
    double new_sigma = (Esigma(0, 0) + Esigma(1, 1)) / 2.0;
    Matrix33 Sigma(new_sigma, 0, 0,
                   0, new_sigma, 0,
                   0, 0, 0);
    std::cout << "Sigma for enforcing essential matrix properties:\n";
    printMatrix(Sigma, "Sigma");
    E = EU * Sigma * EV.transpose();
    printMatrix(E, "E (after enforcing properties)");

    // 再次分解 E 以求出 R 和 t
    svd_decompose(E, EU, Esigma, EV);
    printMatrix(EU, "EU (recomputed)");
    printMatrix(Esigma, "Esigma (recomputed)");
    printMatrix(EV, "EV (recomputed)");

    Matrix33 W(0, -1, 0,
               1, 0, 0,
               0, 0, 1);
    Matrix33 r1 = EU * W * EV.transpose();
    Matrix33 r2 = EU * W.transpose() * EV.transpose();
    double determinant_r1 = determinant(r1);
    double determinant_r2 = determinant(r2);
    Matrix33 R1 = determinant_r1 * r1;
    Matrix33 R2 = determinant_r2 * r2;
    printMatrix(R1, "Candidate Rotation R1");
    printMatrix(R2, "Candidate Rotation R2");

    // t 为 EU 的最后一列
    Vector3D t1 = EU.get_column(EU.cols() - 1);
    Vector3D t2 = -EU.get_column(EU.cols() - 1);
    printVector3D(t1, "Candidate translation t1");
    printVector3D(t2, "Candidate translation t2");

    // 9. 构造投影矩阵和进行三角化
    auto getProjectionMatrix = [&](const Matrix33& Kmat, const Matrix33& Rmat, const Vector3D& tvec) -> Matrix34
    {
        Matrix34 RT(3, 4, 0.0);
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                RT(row, col) = Rmat(row, col);
            }
            RT(row, 3) = tvec[row];
        }
        Matrix34 P = Kmat * RT;
        printMatrix(P, "Projection Matrix P");
        return P;
    };

    auto triangulateDLT = [&](const Matrix34& P1, const Matrix34& P2,
                                const Vector2D& pt1, const Vector2D& pt2) -> Vector3D
    {
        Matrix A_(4, 4, 0.0);
        for (int c = 0; c < 4; c++){
            A_(0, c) = pt1[1] * P1(2, c) - P1(1, c);
        }
        for (int c = 0; c < 4; c++){
            A_(1, c) = P1(0, c) - pt1[0] * P1(2, c);
        }
        for (int c = 0; c < 4; c++){
            A_(2, c) = pt2[1] * P2(2, c) - P2(1, c);
        }
        for (int c = 0; c < 4; c++){
            A_(3, c) = P2(0, c) - pt2[0] * P2(2, c);
        }
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

    // 设第1个摄像机的投影矩阵 P1 (C1与世界坐标系重合)
    Matrix33 I;
    I(0,0) = 1.0; I(0,1) = 0.0; I(0,2) = 0.0;
    I(1,0) = 0.0; I(1,1) = 1.0; I(1,2) = 0.0;
    I(2,0) = 0.0; I(2,1) = 0.0; I(2,2) = 1.0;
    Matrix34 P1_proj = getProjectionMatrix(K, I, Vector3D(0.0, 0.0, 0.0));

    // 候选 (R, t) 组合
    std::vector<std::pair<Matrix33, Vector3D>> candidates = {
        { R1, t1 }, { R1, t2 }, { R2, t1 }, { R2, t2 }
    };

    int best_in_front_count = -1;
    Matrix33 best_R;
    Vector3D best_t;
    std::vector<Vector3D> best_3d_points;

    // 对每个候选进行测试
    for (auto& cnd : candidates) {
        Matrix33 Rc = cnd.first;
        Vector3D tc = cnd.second;
        std::cout << "Testing candidate with R(0,0): " << Rc(0,0)
                  << " and t: (" << tc[0] << ", " << tc[1] << ", " << tc[2] << ")\n";

        Matrix34 P2_proj = getProjectionMatrix(K, Rc, tc);

        int count_in_front = 0;
        std::vector<Vector3D> tmp_points;
        tmp_points.reserve(num_points);

        for (size_t i = 0; i < num_points; i++) {
            Vector3D X = triangulateDLT(P1_proj, P2_proj, points_0[i], points_1[i]);
            tmp_points.push_back(X);

            bool front_c1 = (X[2] > 0.0);
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
