/**
 * Copyright (C) 2015 by Liangliang Nan
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
 * Compute the Frobenius norm of a 3x3 matrix.
 * Frobenius norm = sqrt( sum of squares of all entries )
 */
static double frobenius_norm_33(const Matrix33 &M)
{
    double sum = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            sum += M(i,j) * M(i,j);
    return std::sqrt(sum);
}

/**
 * Compute the determinant of a 3x3 matrix.
 */
static double determinant_33(const Matrix33 &M)
{
    double a = M(0,0), b = M(0,1), c = M(0,2);
    double d = M(1,0), e = M(1,1), f = M(1,2);
    double g = M(2,0), h = M(2,1), i = M(2,2);
    return a * (e*i - f*h) - b * (d*i - f*g) + c * (d*h - e*g);
}

/**
 * Normalize 2D image points to improve numerical stability in the 8-point algorithm.
 * 1) Shift the mean of the points to the origin.
 * 2) Scale them so that the average distance from the origin is sqrt(2).
 * 3) Construct the normalization matrix T accordingly.
 */
static void normalizePoints(const std::vector<Vector2D>& pts,
                            std::vector<Vector3D>& pts_h, Matrix33& T)
{
    int N = pts.size();
    double mean_x = 0, mean_y = 0;
    for (int i = 0; i < N; ++i) {
        mean_x += pts[i][0];
        mean_y += pts[i][1];
    }
    mean_x /= N;
    mean_y /= N;

    double avg_dist = 0;
    for (int i = 0; i < N; ++i) {
        double dx = pts[i][0] - mean_x;
        double dy = pts[i][1] - mean_y;
        avg_dist += std::sqrt(dx * dx + dy * dy);
    }
    avg_dist /= N;

    double scale = std::sqrt(2.0) / avg_dist;
    T = Matrix33(scale,      0, -scale * mean_x,
                 0,      scale, -scale * mean_y,
                 0,          0,              1);

    pts_h.resize(N);
    for (int i = 0; i < N; ++i) {
        // Convert to homogeneous coordinates and apply normalization T.
        pts_h[i] = pts[i].homogeneous();
        pts_h[i] = T * pts_h[i];
    }
}

/**
 * Step #1: Estimate the fundamental matrix using the normalized 8-point algorithm.
 *         - Normalize the 2D points
 *         - Formulate and solve the linear system (SVD)
 *         - Enforce the rank-2 constraint on F
 *         - Denormalize F
 */
static Matrix33 estimateFundamentalMatrix(const std::vector<Vector2D>& pts1,
                                          const std::vector<Vector2D>& pts2)
{
    int N = pts1.size();
    if (N < 8) {
        std::cerr << "Error: Not enough points to compute the fundamental matrix." << std::endl;
        return Matrix33();
    }

    // 1) Normalize the input point sets.
    std::vector<Vector3D> pts1_norm, pts2_norm;
    Matrix33 T1, T2;
    normalizePoints(pts1, pts1_norm, T1);
    normalizePoints(pts2, pts2_norm, T2);

    // 2) Construct the linear system A * f = 0, with f being the 9D vector of F.
    Matrix A(N, 9, 0.0);
    for (int i = 0; i < N; ++i) {
        double x1 = pts1_norm[i][0], y1 = pts1_norm[i][1];
        double x2 = pts2_norm[i][0], y2 = pts2_norm[i][1];
        A(i,0) = x2 * x1;
        A(i,1) = x2 * y1;
        A(i,2) = x2;
        A(i,3) = y2 * x1;
        A(i,4) = y2 * y1;
        A(i,5) = y2;
        A(i,6) = x1;
        A(i,7) = y1;
        A(i,8) = 1;
    }

    // 3) Solve A * f = 0 using SVD.
    Matrix U, S, V;
    svd_decompose(A, U, S, V);
    Matrix Vt = V.transpose();
    // The solution is the last row of V^T (or last column of V).
    Vector f = Vt.get_row(Vt.rows() - 1);

    // Reshape the 9D vector f into a 3x3 matrix F_norm.
    Matrix33 F_norm;
    for (int i = 0; i < 9; ++i)
        F_norm(i/3, i%3) = f[i];

    // 4) Enforce rank(F) = 2 by zeroing the smallest singular value.
    Matrix33 Uf, Sf, Vtf;
    svd_decompose(F_norm, Uf, Sf, Vtf);
    Matrix33 Vt_f = Vtf.transpose();
    Sf(2,2) = 0.0;  // Force rank to 2
    Matrix33 F_rank2 = Uf * Sf * Vt_f;

    // 5) Denormalize: F = T2^T * F_rank2 * T1
    Matrix33 F = T2.transpose() * F_rank2 * T1;

    // 6) Optionally, normalize the final F to have a unit Frobenius norm for consistency.
    double normF = frobenius_norm_33(F);
    if (normF != 0)
        F = F / normF;

    return F;
}

/**
 * Construct the 3x4 projection matrix from camera intrinsic K, rotation R, and translation t.
 */
static Matrix34 constructProjectionMatrix(const Matrix33& K, const Matrix33& R, const Vector3D& t)
{
    Matrix34 P(3,4, 0.0);

    // P = K * [R | t], with dimension 3x4.
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double sum = 0;
            for (int k = 0; k < 3; ++k)
                sum += K(i,k) * R(k,j);
            P(i,j) = sum;
        }

    for (int i = 0; i < 3; ++i) {
        double sum = 0;
        for (int k = 0; k < 3; ++k)
            sum += K(i,k) * t[k];
        P(i,3) = sum;
    }
    return P;
}

/**
 * Linear triangulation for a single point:
 *   Solve for X by minimizing ||A * X|| where A is derived from the cross product constraints.
 */
static Vector4D triangulatePoint(const Matrix34& P1, const Matrix34& P2,
                                 const Vector2D& x1, const Vector2D& x2)
{
    // Each point gives two equations per camera, so we have a 4x4 system.
    Matrix A(4, 4, 0.0);

    Vector row0 = x1[0] * P1.get_row(2) - P1.get_row(0);
    Vector row1 = x1[1] * P1.get_row(2) - P1.get_row(1);
    Vector row2 = x2[0] * P2.get_row(2) - P2.get_row(0);
    Vector row3 = x2[1] * P2.get_row(2) - P2.get_row(1);

    A.set_row(0, row0);
    A.set_row(1, row1);
    A.set_row(2, row2);
    A.set_row(3, row3);

    // Solve using SVD.
    Matrix U, S, V;
    svd_decompose(A, U, S, V);
    Matrix Vt = V.transpose();
    Vector X = Vt.get_row(Vt.rows() - 1);
    return X; // Homogeneous 4D representation.
}

/**
 * Compute the cheirality (number of points in front of both cameras).
 * This function also returns the 3D coordinates of the points.
 */
static int computeCheirality(const Matrix34& P1, const Matrix34& P2,
                             const std::vector<Vector2D>& pts1,
                             const std::vector<Vector2D>& pts2,
                             std::vector<Vector3D>& points3d)
{
    int count = 0;
    points3d.clear();
    int N = pts1.size();

    for (int i = 0; i < N; ++i) {
        Vector4D Xh = triangulatePoint(P1, P2, pts1[i], pts2[i]);
        Vector3D X(Xh[0], Xh[1], Xh[2]);
        if (std::fabs(Xh[3]) > 1e-12)  // Avoid division by zero
            X = X / Xh[3];

        points3d.push_back(X);

        // Check depth w.r.t. the first camera
        if (X[2] <= 0)
            continue;

        // Check depth w.r.t. the second camera
        Vector4D Xh2 = P2 * X.homogeneous();
        if (Xh2[2] > 0)
            count++;
    }
    return count;
}

/**
 * Multiply a matrix M by a vector v (for the LM routine).
 */
static std::vector<double> multiplyMatrixVector(const Matrix &M, const std::vector<double> &v)
{
    int m = M.rows(), n = M.cols();
    std::vector<double> res(m, 0.0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            res[i] += M(i,j) * v[j];
    return res;
}

/**
 * Solve a linear system (A * x = b) using Gaussian elimination.
 * Used in the LM refinement.
 */
static bool solveLinearSystem(const Matrix &A_input, const std::vector<double> &b_input, std::vector<double> &x)
{
    int n = A_input.rows();
    Matrix A = A_input; // Make a copy.
    std::vector<double> b = b_input;
    x.resize(n, 0.0);
    const double eps = 1e-12;

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        // Find pivot row
        int pivot = i;
        for (int j = i+1; j < n; ++j) {
            if (std::fabs(A(j,i)) > std::fabs(A(pivot,i)))
                pivot = j;
        }
        if (std::fabs(A(pivot,i)) < eps)
            return false;

        // Swap pivot row into place
        if (pivot != i) {
            for (int k = 0; k < n; ++k)
                std::swap(A(i,k), A(pivot,k));
            std::swap(b[i], b[pivot]);
        }

        // Eliminate below pivot
        for (int j = i+1; j < n; ++j) {
            double factor = A(j,i) / A(i,i);
            for (int k = i; k < n; ++k)
                A(j,k) -= factor * A(i,k);
            b[j] -= factor * b[i];
        }
    }

    // Back-substitution
    for (int i = n-1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i+1; j < n; ++j)
            sum += A(i,j) * x[j];
        x[i] = (b[i] - sum) / A(i,i);
    }
    return true;
}

/**
 * Refine all 3D points simultaneously using nonlinear optimization (Levenberg–Marquardt).
 * Each 2D point is compared with the reprojected point from the 3D guess, and the reprojection errors are minimized.
 */
static void refinePointsLM(const std::vector<Vector2D>& pts0,
                           const std::vector<Vector2D>& pts1,
                           const Matrix34& P1, const Matrix34& P2,
                           std::vector<Vector3D>& points3d)
{
    int N = points3d.size();
    int num_params = 3 * N;
    int num_residuals = 4 * N; // two residuals (x,y) per camera, for two cameras => 4 per point

    std::vector<double> p(num_params, 0.0);
    for (int i = 0; i < N; ++i) {
        p[3*i + 0] = points3d[i][0];
        p[3*i + 1] = points3d[i][1];
        p[3*i + 2] = points3d[i][2];
    }

    // Cost function: r_vec = (u1 - observed_u1, v1 - observed_v1, u2 - observed_u2, v2 - observed_v2, ...)
    auto costFunction = [&](const std::vector<double>& p_vec, std::vector<double>& r_vec) {
        r_vec.resize(num_residuals, 0.0);
        for (int i = 0; i < N; ++i) {
            double x = p_vec[3*i + 0];
            double y = p_vec[3*i + 1];
            double z = p_vec[3*i + 2];
            Vector4D Xh(x, y, z, 1.0);

            // Reprojection on camera 1
            Vector4D proj1 = P1 * Xh;
            double d1 = proj1[2];
            double u1 = proj1[0] / d1;
            double v1 = proj1[1] / d1;
            r_vec[4*i + 0] = u1 - pts0[i][0];
            r_vec[4*i + 1] = v1 - pts0[i][1];

            // Reprojection on camera 2
            Vector4D proj2 = P2 * Xh;
            double d2 = proj2[2];
            double u2 = proj2[0] / d2;
            double v2 = proj2[1] / d2;
            r_vec[4*i + 2] = u2 - pts1[i][0];
            r_vec[4*i + 3] = v2 - pts1[i][1];
        }
    };

    // Compute the Jacobian matrix for LM: partial derivatives of each residual w.r.t. each 3D coordinate
    auto computeJacobian = [&](const std::vector<double>& p_vec, Matrix &J) {
        J = Matrix(num_residuals, num_params, 0.0);
        for (int i = 0; i < N; ++i) {
            double x = p_vec[3*i + 0];
            double y = p_vec[3*i + 1];
            double z = p_vec[3*i + 2];
            Vector4D Xh(x, y, z, 1.0);

            auto computeCamJacobian = [&](const Matrix34 &P, int offset) {
                Vector4D proj = P * Xh;
                double d = proj[2];
                double u = proj[0] / d;
                double v = proj[1] / d;

                Vector P0 = P.get_row(0);
                Vector P1 = P.get_row(1);
                Vector P2 = P.get_row(2);

                for (int j = 0; j < 3; ++j) {
                    double du_dX = (P0[j] * d - proj[0] * P2[j]) / (d*d);
                    double dv_dX = (P1[j] * d - proj[1] * P2[j]) / (d*d);
                    J(offset + 0, 3*i + j) = du_dX;
                    J(offset + 1, 3*i + j) = dv_dX;
                }
            };
            // For each point, fill four rows in the Jacobian (two for camera 1, two for camera 2).
            computeCamJacobian(P1, 4*i);
            computeCamJacobian(P2, 4*i + 2);
        }
    };

    double lambda = 1e-3;
    int maxIter = 50;
    double tol = 1e-6;

    // Compute initial residual and cost
    std::vector<double> r, r_new;
    costFunction(p, r);
    double cost = 0.0;
    for (double ri : r)
        cost += ri * ri;

    // Levenberg–Marquardt iterations
    for (int iter = 0; iter < maxIter; ++iter) {
        Matrix J;
        computeJacobian(p, J);
        Matrix JT = J.transpose();
        Matrix A = JT * J;

        // Damping diagonal
        for (int i = 0; i < A.rows(); ++i)
            A(i,i) += lambda;

        // -JT * r
        std::vector<double> JT_r = multiplyMatrixVector(JT, r);
        std::vector<double> b(A.rows(), 0.0);
        for (int i = 0; i < A.rows(); ++i)
            b[i] = -JT_r[i];

        // Solve the linear system A * delta = b
        std::vector<double> delta;
        bool solved = solveLinearSystem(A, b, delta);
        if (!solved) {
            std::cerr << "LM: Failed to solve the linear system." << std::endl;
            break;
        }

        // Update parameters
        std::vector<double> p_new = p;
        for (int i = 0; i < (int)p_new.size(); ++i)
            p_new[i] += delta[i];

        // Compute new cost
        costFunction(p_new, r_new);
        double cost_new = 0.0;
        for (double ri : r_new)
            cost_new += ri * ri;

        // LM acceptance/rejection
        if (cost_new < cost) {
            p = p_new;
            cost = cost_new;
            lambda *= 0.8;
            if (std::sqrt(cost_new) < tol)
                break;
            r = r_new;
        } else {
            lambda *= 2.0;
        }
    }

    // Update the 3D points
    for (int i = 0; i < N; ++i) {
        points3d[i][0] = p[3*i + 0];
        points3d[i][1] = p[3*i + 1];
        points3d[i][2] = p[3*i + 2];
    }
}

/**
 * Step #2: Recover the relative pose (R, t) from the fundamental matrix.
 *   1) Compute the essential matrix E = K^T * F * K
 *   2) Decompose E using SVD to obtain the candidate rotations (R1, R2) and translation (t)
 *   3) Among the four (R, t) combinations, select the one maximizing the number of points with positive depth.
 */
static void recoverRelativePose(const Matrix33& K, const Matrix33& F,
                                const std::vector<Vector2D>& pts0,
                                const std::vector<Vector2D>& pts1,
                                Matrix33& bestR, Vector3D& bestt)
{
    // Compute E
    Matrix33 E = K.transpose() * F * K;

    // SVD of E
    Matrix U, S_mat, V;
    svd_decompose(E, U, S_mat, V);
    Matrix Vt = V.transpose();

    // Force singular values: rank(E) = 2
    double sigma = (S_mat(0,0) + S_mat(1,1)) / 2.0;
    S_mat(0,0) = sigma;
    S_mat(1,1) = sigma;
    S_mat(2,2) = 0.0;
    E = U * S_mat * Vt;

    // Possible rotation from E using the matrix W
    Matrix33 W(0, -1, 0,
               1,  0, 0,
               0,  0, 1);

    Matrix33 R1 = U * W * Vt;
    Matrix33 R2 = U * W.transpose() * Vt;

    // Ensure the determinant of the rotations is positive
    if (determinant_33(R1) < 0)
        R1 = R1 * (-1.0);
    if (determinant_33(R2) < 0)
        R2 = R2 * (-1.0);

    // The translation is the 3rd column of U
    Vector3D t_candidate = U.get_column(2);

    // Construct four possible pairs
    std::vector<Matrix33> Rs = {R1, R1, R2, R2};
    std::vector<Vector3D> ts = {
        t_candidate,
        t_candidate * (-1.0),
        t_candidate,
        t_candidate * (-1.0)
    };

    // Check which one yields the most points in front of both cameras
    int best_count = -1;

    Matrix33 I_mat(1,0,0, 0,1,0, 0,0,1);
    Vector3D zero(0,0,0);
    Matrix34 P1 = constructProjectionMatrix(K, I_mat, zero);

    for (int i = 0; i < 4; ++i) {
        Matrix34 P2 = constructProjectionMatrix(K, Rs[i], ts[i]);
        std::vector<Vector3D> current_points3d;
        int count = computeCheirality(P1, P2, pts0, pts1, current_points3d);
        if (count > best_count) {
            best_count = count;
            bestR = Rs[i];
            bestt = ts[i];
        }
    }
}

/**
 * Step #3: Linear triangulation of all corresponding image points.
 */
static std::vector<Vector3D> triangulateAllPoints(const Matrix34 &P1, const Matrix34 &P2,
                                                  const std::vector<Vector2D>& pts0,
                                                  const std::vector<Vector2D>& pts1)
{
    std::vector<Vector3D> points3d;
    points3d.reserve(pts0.size());

    for (size_t i = 0; i < pts0.size(); ++i) {
        Vector4D Xh = triangulatePoint(P1, P2, pts0[i], pts1[i]);
        Vector3D X(Xh[0], Xh[1], Xh[2]);
        if (std::fabs(Xh[3]) > 1e-12)  // Avoid division by zero
            X = X / Xh[3];
        points3d.push_back(X);
    }
    return points3d;
}

/**
 * Step #4 (optional): Evaluate the reconstructed points.
 * (You can implement additional functions to measure and assess the quality of the 3D reconstruction.)
 */

// --------------------------------------------------------------
// Main triangulation function that calls the above sub-tasks.
// --------------------------------------------------------------
bool Triangulation::triangulation(
        double fx, double fy,       // focal lengths (same for both cameras)
        double cx, double cy,       // principal point (same for both cameras)
        double s,                   // skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  // 2D points in the first image
        const std::vector<Vector2D> &points_1,  // 2D points in the second image
        std::vector<Vector3D> &points_3d,       // output: reconstructed 3D points
        Matrix33 &R,                // output: recovered rotation of the second camera
        Vector3D &t                 // output: recovered translation of the second camera
) const
{
    // Basic checks
    if (points_0.size() < 8 || points_0.size() != points_1.size()) {
        std::cerr << "Invalid input: need at least 8 corresponding points and the same number in both views." << std::endl;
        return false;
    }

    // Construct intrinsic matrix K
    Matrix33 K(fx, s,  cx,
               0,  fy, cy,
               0,   0,  1);

    // Step #1: Estimate the fundamental matrix F using the normalized 8-point algorithm
    Matrix33 F = estimateFundamentalMatrix(points_0, points_1);

    // Step #2: Recover (R, t) from F
    recoverRelativePose(K, F, points_0, points_1, R, t);

    // Build the two projection matrices
    Matrix33 I_mat(1,0,0, 0,1,0, 0,0,1);
    Vector3D zero(0,0,0);
    Matrix34 P1 = constructProjectionMatrix(K, I_mat, zero);
    Matrix34 P2 = constructProjectionMatrix(K, R, t);

    // Step #3: Triangulate all corresponding points
    points_3d = triangulateAllPoints(P1, P2, points_0, points_1);

    // (Optional) Nonlinear refinement of 3D points
    refinePointsLM(points_0, points_1, P1, P2, points_3d);

    return (!points_3d.empty());
}
