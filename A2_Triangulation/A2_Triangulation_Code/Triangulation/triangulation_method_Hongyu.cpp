#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;

//--------------------------------------------
// **Step #2: 确定正确的 R 和 t**
//--------------------------------------------
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
        // 构造投影矩阵
        Matrix34 P1(K(0,0), K(0,1), K(0,2), 0,
                    K(1,0), K(1,1), K(1,2), 0,
                    K(2,0), K(2,1), K(2,2), 0);

        Matrix34 P2(K(0,0) * Rs[i](0,0) + K(0,1) * Rs[i](1,0) + K(0,2) * Rs[i](2,0),
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
                    K(2,0) * ts[i][0]  + K(2,1) * ts[i][1]  + K(2,2) * ts[i][2]);

        // 计算三角测量点的正深度个数
        int count = 0;
        for (size_t j = 0; j < pts1.size(); ++j) {
            Vector4D Xh = triangulatePoint(P1, P2, pts1[j], pts2[j]);
            Vector3D X(Xh[0] / Xh[3], Xh[1] / Xh[3], Xh[2] / Xh[3]);

            if (X[2] > 0) // 确保点位于相机前方
                count++;
        }

        // 选择使最多点位于相机前方的姿态
        if (count > best_count) {
            best_count = count;
            best_R = Rs[i];
            best_t = ts[i];
        }
    }
}

//--------------------------------------------
// **Step #3.1: 计算投影矩阵**
//--------------------------------------------
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

//--------------------------------------------
// **Step #3.2: 线性三角测量**
//--------------------------------------------
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

//--------------------------------------------
// **Step #3.3: 非线性优化 3D 点**
//--------------------------------------------
static void refine_3D_points(const Matrix34& P1, const Matrix34& P2,
                             const std::vector<Vector2D>& pts1,
                             const std::vector<Vector2D>& pts2,
                             std::vector<Vector3D>& points_3d)
{
    easy3d::OptimizerLM optimizer;
    optimizer.set_num_parameters(points_3d.size() * 3);
    optimizer.set_num_residuals(points_3d.size() * 4);

    optimizer.set_function([&](const std::vector<double>& params, std::vector<double>& residuals) {
        for (size_t i = 0; i < points_3d.size(); ++i) {
            Vector3D X(params[i * 3], params[i * 3 + 1], params[i * 3 + 2]);

            Vector4D Xh1 = P1 * X.homogeneous();
            Vector4D Xh2 = P2 * X.homogeneous();

            residuals[i * 4]     = Xh1[0] / Xh1[2] - pts1[i][0];
            residuals[i * 4 + 1] = Xh1[1] / Xh1[2] - pts1[i][1];
            residuals[i * 4 + 2] = Xh2[0] / Xh2[2] - pts2[i][0];
            residuals[i * 4 + 3] = Xh2[1] / Xh2[2] - pts2[i][1];
        }
    });

    std::vector<double> params(points_3d.size() * 3);
    for (size_t i = 0; i < points_3d.size(); ++i) {
        params[i * 3] = points_3d[i][0];
        params[i * 3 + 1] = points_3d[i][1];
        params[i * 3 + 2] = points_3d[i][2];
    }

    optimizer.optimize(params);

    for (size_t i = 0; i < points_3d.size(); ++i) {
        points_3d[i] = Vector3D(params[i * 3], params[i * 3 + 1], params[i * 3 + 2]);
    }
}