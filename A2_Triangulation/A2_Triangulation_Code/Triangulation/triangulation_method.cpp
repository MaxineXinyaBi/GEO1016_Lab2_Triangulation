#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace easy3d;

//------------------------------------------------------------------------------
// 辅助函数：打印矩阵
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// 辅助函数：打印 Vector3D
//------------------------------------------------------------------------------
void printVector3D(const Vector3D& v, const std::string& name) {
    std::cout << name << " = (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
}

//------------------------------------------------------------------------------
// 新增：基于 DLT 的三角测量（求解齐次 3D 点）
//------------------------------------------------------------------------------
static Vector4D triangulatePoint(const Matrix34& P1, const Matrix34& P2,
                                 const Vector2D& pt1, const Vector2D& pt2)
{
    Matrix A(4, 4, 0.0);
    // 对第一幅图像
    for (int c = 0; c < 4; c++){
        A(0, c) = pt1[1] * P1.get_row(2)[c] - P1.get_row(1)[c];
    }
    for (int c = 0; c < 4; c++){
        A(1, c) = P1.get_row(0)[c] - pt1[0] * P1.get_row(2)[c];
    }
    // 对第二幅图像
    for (int c = 0; c < 4; c++){
        A(2, c) = pt2[1] * P2.get_row(2)[c] - P2.get_row(1)[c];
    }
    for (int c = 0; c < 4; c++){
        A(3, c) = P2.get_row(0)[c] - pt2[0] * P2.get_row(2)[c];
    }
    std::cout << "triangulatePoint: Matrix A:" << std::endl;
    printMatrix(A, "A");
    Matrix U, S, V;
    svd_decompose(A, U, S, V);
    Matrix Vt = V.transpose();
    Vector X = Vt.get_row(Vt.rows()-1);
    std::cout << "triangulatePoint: Homogeneous solution X = " << X << std::endl;
    return X;
}

//------------------------------------------------------------------------------
// Step #2: 确定正确的 R 和 t
//------------------------------------------------------------------------------
static void determine_correct_pose(const Matrix33& K,
                                   const Matrix33& R1, const Matrix33& R2,
                                   const Vector3D& t,
                                   const std::vector<Vector2D>& pts1,
                                   const std::vector<Vector2D>& pts2,
                                   Matrix33& best_R, Vector3D& best_t)
{
    std::vector<Matrix33> Rs = {R1, R1, R2, R2};
    std::vector<Vector3D> ts = {t, -t, t, -t};

    int best_count = -1;
    for (int i = 0; i < 4; ++i) {
        // 构造投影矩阵 P1 = K [I|0]
        Matrix34 P1(K(0,0), K(0,1), K(0,2), 0,
                    K(1,0), K(1,1), K(1,2), 0,
                    K(2,0), K(2,1), K(2,2), 0);

        // 构造投影矩阵 P2 = K [R|t]
        Matrix34 P2(
            K(0,0) * Rs[i](0,0) + K(0,1) * Rs[i](1,0) + K(0,2) * Rs[i](2,0),
            K(0,0) * Rs[i](0,1) + K(0,1) * Rs[i](1,1) + K(0,2) * Rs[i](2,1),
            K(0,0) * Rs[i](0,2) + K(0,1) * Rs[i](1,2) + K(0,2) * Rs[i](2,2),
            K(0,0) * ts[i][0]  + K(0,1) * ts[i][1]  + K(0,2) * ts[i][2],
            K(1,0) * Rs[i](0,0) + K(1,1) * Rs[i](1,0) + K(1,2) * Rs[i](2,0),
            K(1,0) * Rs[i](0,1) + K(1,1) * Rs[i](1,1) + K(1,2) * Rs[i](2,1),
            K(1,0) * Rs[i](0,2) + K(1,1) * Rs[i](1,2) + K(1,2) * Rs[i](2,2),
            K(1,0) * ts[i][0]  + K(1,1) * ts[i][1]  + K(1,2) * ts[i][2],
            K(2,0) * Rs[i](0,0) + K(2,1) * Rs[i](1,0) + K(2,2) * Rs[i](2,0),
            K(2,0) * Rs[i](0,1) + K(2,1) * Rs[i](1,1) + K(2,2) * Rs[i](2,1),
            K(2,0) * Rs[i](0,2) + K(2,1) * Rs[i](1,2) + K(2,2) * Rs[i](2,2),
            K(2,0) * ts[i][0]  + K(2,1) * ts[i][1]  + K(2,2) * ts[i][2]
        );

        int count = 0;
        for (size_t j = 0; j < pts1.size(); ++j) {
            Vector4D Xh = triangulatePoint(P1, P2, pts1[j], pts2[j]);
            Vector3D X(Xh[0] / Xh[3], Xh[1] / Xh[3], Xh[2] / Xh[3]);
            if (X[2] > 0)
                count++;
        }
        if (count > best_count) {
            best_count = count;
            best_R = Rs[i];
            best_t = ts[i];
        }
    }
}

//------------------------------------------------------------------------------
// Step #3.1: 计算投影矩阵：P = K [R|t]
//------------------------------------------------------------------------------
static Matrix34 compute_projection_matrix(const Matrix33& K, const Matrix33& R, const Vector3D& t)
{
    Matrix34 P;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            P(i, j) = K(i, 0) * R(0, j) + K(i, 1) * R(1, j) + K(i, 2) * R(2, j);
        }
        P(i, 3) = K(i, 0) * t[0] + K(i, 1) * t[1] + K(i, 2) * t[2];
    }
    return P;
}

//------------------------------------------------------------------------------
// Step #3.2: 线性三角测量（计算所有对应点的 3D 坐标）
//------------------------------------------------------------------------------
static void triangulate_all_points(const Matrix34& P1, const Matrix34& P2,
                                   const std::vector<Vector2D>& pts1,
                                   const std::vector<Vector2D>& pts2,
                                   std::vector<Vector3D>& points_3d)
{
    points_3d.clear();
    for (size_t i = 0; i < pts1.size(); ++i) {
        Vector4D Xh = triangulatePoint(P1, P2, pts1[i], pts2[i]);
        points_3d.emplace_back(Xh[0] / Xh[3], Xh[1] / Xh[3], Xh[2] / Xh[3]);
    }
}

//------------------------------------------------------------------------------
// Step #3.3: [Optional] 非线性最小二乘优化对 3D 点进行精细化
//------------------------------------------------------------------------------
class MyObjective : public Objective_LM {
public:
    // 显式调用 Objective_LM 构造函数，传入残差数量和变量数量
    MyObjective(const Matrix34& P1, const Matrix34& P2,
                const std::vector<Vector2D>& pts1,
                const std::vector<Vector2D>& pts2,
                size_t num_points)
        : Objective_LM(static_cast<int>(num_points * 4), static_cast<int>(num_points * 3)),
          P1_(P1), P2_(P2), pts1_(pts1), pts2_(pts2), num_points_(num_points) {}

    virtual bool evaluate(const std::vector<double>& params, std::vector<double>& residuals) {
        residuals.resize(num_points_ * 4);
        for (size_t i = 0; i < num_points_; ++i) {
            Vector3D X(params[i*3], params[i*3+1], params[i*3+2]);
            Vector4D Xh1 = P1_ * X.homogeneous();
            Vector4D Xh2 = P2_ * X.homogeneous();
            double proj1_x = Xh1[0] / Xh1[2];
            double proj1_y = Xh1[1] / Xh1[2];
            double proj2_x = Xh2[0] / Xh2[2];
            double proj2_y = Xh2[1] / Xh2[2];
            residuals[i*4]     = proj1_x - pts1_[i][0];
            residuals[i*4 + 1] = proj1_y - pts1_[i][1];
            residuals[i*4 + 2] = proj2_x - pts2_[i][0];
            residuals[i*4 + 3] = proj2_y - pts2_[i][1];
        }
        return true;
    }
private:
    Matrix34 P1_, P2_;
    std::vector<Vector2D> pts1_, pts2_;
    size_t num_points_;
};

static void refine_3D_points(const Matrix34& P1, const Matrix34& P2,
                             const std::vector<Vector2D>& pts1,
                             const std::vector<Vector2D>& pts2,
                             std::vector<Vector3D>& points_3d)
{
    size_t N = points_3d.size();
    MyObjective objective(P1, P2, pts1, pts2, N);
    std::vector<double> params(N * 3);
    for (size_t i = 0; i < N; ++i) {
        params[i * 3]     = points_3d[i][0];
        params[i * 3 + 1] = points_3d[i][1];
        params[i * 3 + 2] = points_3d[i][2];
    }
    easy3d::Optimizer_LM optimizer;
    optimizer.optimize(&objective, params, nullptr);
    for (size_t i = 0; i < N; ++i) {
        points_3d[i] = Vector3D(params[i * 3], params[i * 3 + 1], params[i * 3 + 2]);
    }
}

//------------------------------------------------------------------------------
// 主函数：Triangulation::triangulation()
//------------------------------------------------------------------------------
bool Triangulation::triangulation(
    double fx, double fy,     /// input: the focal lengths (same for both cameras)
    double cx, double cy,     /// input: the principal point (same for both cameras)
    double s,                 /// input: the skew factor (same for both cameras)
    const std::vector<Vector2D>& points_0,  /// input: 2D image points in the 1st image.
    const std::vector<Vector2D>& points_1,  /// input: 2D image points in the 2nd image.
    std::vector<Vector3D>& points_3d,       /// output: reconstructed 3D points
    Matrix33& R,   /// output: recovered rotation of 2nd camera
    Vector3D& t    /// output: recovered translation of 2nd camera
) const
{
    std::cout << "\nStarting triangulation process...\n\n";

    if (points_0.size() < 8 || points_0.size() != points_1.size()) {
        std::cerr << "Invalid input: not enough points or point numbers don't match.\n";
        return false;
    }
    size_t num_points = points_0.size();

    // 这里假设已经计算出 F、E 及候选 R1, R2, t_candidate（此处用占位数据演示）
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);
    printMatrix(K, "K");

    // 占位候选解：使用单位矩阵作为候选旋转，t_candidate 为 (0.1, 0, 0)
    Matrix33 R1 = Matrix33::identity(1.0);
    Matrix33 R2 = Matrix33::identity(1.0);
    Vector3D t_candidate(0.1, 0.0, 0.0);

    // --- Step #2: 确定正确的 R 和 t ---
    Matrix33 best_R;
    Vector3D best_t;
    determine_correct_pose(K, R1, R2, t_candidate, points_0, points_1, best_R, best_t);
    R = best_R;
    t = best_t;
    std::cout << "Selected relative pose:" << std::endl;
    printMatrix(R, "R");
    printVector3D(t, "t");

    // --- Step #3.1: 计算投影矩阵 ---
    Matrix34 P1_proj = compute_projection_matrix(K, Matrix33::identity(1.0), Vector3D(0,0,0)); // P1 = K[I|0]
    Matrix34 P2_proj = compute_projection_matrix(K, R, t);  // P2 = K[R|t]

    // --- Step #3.2: 线性三角测量 ---
    triangulate_all_points(P1_proj, P2_proj, points_0, points_1, points_3d);
    std::cout << "After linear triangulation, number of points: " << points_3d.size() << std::endl;

    // --- Step #3.3: [Optional] 非线性优化 3D 点 ---
    refine_3D_points(P1_proj, P2_proj, points_0, points_1, points_3d);
    std::cout << "After non-linear refinement, number of points: " << points_3d.size() << std::endl;

    return points_3d.size() > 0;
}
