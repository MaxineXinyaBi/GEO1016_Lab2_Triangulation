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

//------------------------------------------------------------------------------
// Function: Calculate the Frobenius norm of a 3x3 matrix
//------------------------------------------------------------------------------
static double frobenius_norm_33(const Matrix33 &M)
{
    double sum = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            sum += M(i,j) * M(i,j);
    return std::sqrt(sum);
}


//------------------------------------------------------------------------------
// Function: Normalize 2D points (convert to homogeneous coordinates)
//------------------------------------------------------------------------------
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
    
    // Construct the normalized matrix T
    T = Matrix33(scale,      0, -scale * mean_x,
                 0,      scale, -scale * mean_y,
                 0,          0,              1);
    pts_h.resize(N);
    for (int i = 0; i < N; ++i) {
        pts_h[i] = pts[i].homogeneous();  // Convert to homogeneous coordinates
        pts_h[i] = T * pts_h[i];          // Normalization
    }
}

//------------------------------------------------------------------------------
// Function: Calculate the basic matrix F (normalized 8-point algorithm)
//------------------------------------------------------------------------------
static Matrix33 computeFundamentalMatrix(const std::vector<Vector2D>& pts1,
                                           const std::vector<Vector2D>& pts2)
{
    int N = pts1.size();
    if (N < 8) {
        std::cerr << "Not enough points to compute F" << std::endl;
        return Matrix33();
    }
    std::vector<Vector3D> pts1_norm, pts2_norm;
    Matrix33 T1, T2;
    normalizePoints(pts1, pts1_norm, T1);
    normalizePoints(pts2, pts2_norm, T2);

    // Construct the design matrix A (N x 9)
    Matrix A(N, 9, 0.0);
    for (int i = 0; i < N; ++i) {
        double x1 = pts1_norm[i][0];
        double y1 = pts1_norm[i][1];
        double x2 = pts2_norm[i][0];
        double y2 = pts2_norm[i][1];
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

    // SVD decomposition    
    Matrix U, S, V;
    U.resize(N, N);
    S.resize(N, 9);
    V.resize(9, 9);
    svd_decompose(A, U, S, V);

    // Compute the transpose of V
    Matrix Vt = V.transpose();
    Vector f = Vt.get_row(Vt.rows()-1);
    Matrix33 F_norm;
    for (int i = 0; i < 9; ++i) {
        F_norm(i/3, i%3) = f[i];
    }

    // Force F to be rank 2: Do SVD on F_norm and set the smallest singular value to zero
    Matrix33 Uf, Sf, Vtf;
    svd_decompose(F_norm, Uf, Sf, Vtf);

    // Vtf is V instead of Vᵀ, transposed again
    Matrix33 Vt_f = Vtf.transpose();
    Sf(2,2) = 0;
    Matrix33 F_norm2 = Uf * Sf * Vt_f;

    // Denormalization
    Matrix33 F = T2.transpose() * F_norm2 * T1;
    double normF = frobenius_norm_33(F);
    if (normF != 0) F = F / normF;
    return F;
}


//------------------------------------------------------------------------------
// Function: Construct projection matrix P = K * [R | t]
//------------------------------------------------------------------------------
static Matrix34 constructProjectionMatrix(const Matrix33& K, const Matrix33& R, const Vector3D& t)
{
    Matrix34 P(3,4, 0.0);
    
    // P[0:2, 0:3] = K * R
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double sum = 0;
            for (int k = 0; k < 3; ++k)
                sum += K(i,k) * R(k,j);
            P(i,j) = sum;
        }
    
    // P[0:2, 3] = K * t
    for (int i = 0; i < 3; ++i) {
        double sum = 0;
        for (int k = 0; k < 3; ++k)
            sum += K(i,k) * t[k];
        P(i,3) = sum;
    }
    return P;
}

//------------------------------------------------------------------------------
// Function: Triangulation (DLT) of a single point
//------------------------------------------------------------------------------
static Vector4D triangulatePoint(const Matrix34& P1, const Matrix34& P2,
                                   const Vector2D& x1, const Vector2D& x2)
{
    Matrix A(4, 4, 0.0);
    Vector row0 = x1[0] * P1.get_row(2) - P1.get_row(0);
    Vector row1 = x1[1] * P1.get_row(2) - P1.get_row(1);
    Vector row2 = x2[0] * P2.get_row(2) - P2.get_row(0);
    Vector row3 = x2[1] * P2.get_row(2) - P2.get_row(1);
    A.set_row(0, row0);
    A.set_row(1, row1);
    A.set_row(2, row2);
    A.set_row(3, row3);

    Matrix44 U, S, V;
    svd_decompose(A, U, S, V);

    Matrix Vt = V.transpose();
    Vector X = Vt.get_row(Vt.rows()-1);
    return X;
}

//------------------------------------------------------------------------------
// Function: calculate the positive position of the candidate solution, return the number of positive depth points, and output all 3D points at the same time
//------------------------------------------------------------------------------
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
        
        // Non-homogeneous coordinates: X = (Xh[0], Xh[1], Xh[2]) / Xh[3]
        Vector3D X(Xh[0], Xh[1], Xh[2]);
        X = X / Xh[3];
        points3d.push_back(X);

        // For the first camera (P1 = K[I|0]), the depth is X[2]
        if (X[2] <= 0)
            continue;

        // For the second camera, calculate P2 * X.homogeneous()
        Vector4D Xh2 = P2 * X.homogeneous();
        if (Xh2[2] > 0)
            count++;
    }
    return count;
}

//------------------------------------------------------------------------------
// Main function：Triangulation::triangulation()
//------------------------------------------------------------------------------
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        double s,                 /// input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    std::cout << "\nTODO: implement the 'triangulation()' function in the file 'Triangulation/triangulation_method.cpp'\n\n";

    if (points_0.size() < 8 || points_0.size() != points_1.size()) {
        std::cerr << "Invalid input: need at least 8 corresponding points and equal number in both views." << std::endl;
        return false;
    }

    // Construct the camera intrinsic parameter matrix K
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0,  0,  1);

    // ----------------------------------------------------------------------
    // Step 1: Estimate the fundamental matrix F (normalized 8-point algorithm)
    // ----------------------------------------------------------------------
    Matrix33 F = computeFundamentalMatrix(points_0, points_1);
    std::cout << "Computed Fundamental Matrix F:" << std::endl;
    std::cout << F << std::endl;

    // ----------------------------------------------------------------------
    // Step 2: Calculate the essential matrix E = K^T * F * K
    // ----------------------------------------------------------------------
    Matrix33 E = K.transpose() * F * K;
    std::cout << "Computed Essential Matrix E:" << std::endl;
    std::cout << E << std::endl;

    // ----------------------------------------------------------------------
    // Step 3: Perform SVD decomposition on E and force its singular values ​​to be [sigma, sigma, 0]
    // ----------------------------------------------------------------------
    Matrix33 U, S_mat, V;
    svd_decompose(E, U, S_mat, V);
    Matrix Vt = V.transpose();
    double sigma = (S_mat(0,0) + S_mat(1,1)) / 2.0;
    S_mat(0,0) = sigma; S_mat(1,1) = sigma; S_mat(2,2) = 0;
    E = U * S_mat * Vt;

    // ----------------------------------------------------------------------
    // Step 4: Recovering Candidate Camera Pose
    // ----------------------------------------------------------------------
    Matrix33 W(0, -1, 0,
               1,  0, 0,
               0,  0, 1);
    Matrix33 R1 = U * W * Vt;
    Matrix33 R2 = U * W.transpose() * Vt;
    if (determinant(R1) < 0)
        R1 = R1 * (-1.0);
    if (determinant(R2) < 0)
        R2 = R2 * (-1.0);
    Vector3D t_candidate = U.get_column(2); // Direction only

    std::vector<Matrix33> Rs = {R1, R1, R2, R2};
    std::vector<Vector3D> ts = {t_candidate, t_candidate * (-1.0), t_candidate, t_candidate * (-1.0)};

    int best_count = -1;
    Matrix33 bestR;
    Vector3D bestt;
    std::vector<Vector3D> best_points3d;
    Matrix33 I_mat(1,0,0, 0,1,0, 0,0,1);
    Vector3D zero(0,0,0);
    Matrix34 P1 = constructProjectionMatrix(K, I_mat, zero);
    for (int i = 0; i < 4; ++i) {
        Matrix34 P2 = constructProjectionMatrix(K, Rs[i], ts[i]);
        std::vector<Vector3D> current_points3d;
        int count = computeCheirality(P1, P2, points_0, points_1, current_points3d);
        std::cout << "Candidate " << i+1 << ": cheirality count = " << count << std::endl;
        if (count > best_count) {
            best_count = count;
            bestR = Rs[i];
            bestt = ts[i];
            best_points3d = current_points3d;
        }
    }
    R = bestR;
    t = bestt;
    std::cout << "Recovered Rotation R:" << std::endl;
    std::cout << R << std::endl;
    std::cout << "Recovered Translation t:" << std::endl;
    std::cout << t << std::endl;

    // ----------------------------------------------------------------------
    // Step 5: Reconstruct all 3D points (triangulation)
    // ----------------------------------------------------------------------
    Matrix34 P2_final = constructProjectionMatrix(K, R, t);
    points_3d.clear();
    int N_pts = points_0.size();
    for (int i = 0; i < N_pts; ++i) {
        Vector4D Xh = triangulatePoint(P1, P2_final, points_0[i], points_1[i]);
        Vector3D X(Xh[0], Xh[1], Xh[2]);
        X = X / Xh[3];
        points_3d.push_back(X);
    }
    return points_3d.size() > 0;
}
