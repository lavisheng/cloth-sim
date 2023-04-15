#include "w_matrix.h"

void precompute_duv(Eigen::MatrixXd &global_duv, const Eigen::MatrixXi &F, Eigen::MatrixXd UV) {
  global_duv.resize(2 * F.rows(), 2);
  for (int i = 0; i < F.rows(); i++) {
    Eigen::RowVector3i tri = F.row(i);
    // we compute the duv we use to compute wuv here
    global_duv.block<2, 2>(i * 2, 0)
            << UV(i * 3 + tri(1), 0) - UV(i * 3 + tri(0), 0),
            UV(i * 3 + tri(2), 0) - UV(i * 3 + tri(0), 0),
            UV(i * 3 + tri(1), 1) - UV(i * 3 + tri(1), 1),
            UV(i * 3 + tri(2), 1) - UV(i * 3 + tri(0), 1);
  }
}

void wuv(Eigen::MatrixXd &wuv, Eigen::MatrixXi F, Eigen::MatrixXd V, int i, Eigen::Matrix2d duv) {
  Eigen::RowVector3i tri = F.row(i);
  Eigen::MatrixXd dx(3, 2);
  // construct dx
  dx.col(0) << V.row(tri(1)) - V.row(tri(0));
  dx.col(1) << V.row(tri(2)) - V.row(tri(0));
  // we calculate [Wu Wv] = dx * duv^-1
  // should be 3 x 2
  wuv = dx * duv.inverse();
}

void precompute_dwdx(Eigen::MatrixXd &dwudx, Eigen::MatrixXd &dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd duv) {
  dwudx.resize(F.rows(), 3);
  dwvdx.resize(F.rows(), 3);
  Eigen::Matrix2d local_duv(2, 2);
  double D;
  for (int i = 0; i < F.rows(); i++) {
    // basically use what we did in w_matrix to get duv
    local_duv = duv.block<2, 2>(i * 2, 0);
    // the denominator du1 * dv2 - du2 * dv1
    double du1 = local_duv(0, 0);
    double du2 = local_duv(0, 1);
    double dv1 = local_duv(1, 0);
    double dv2 = local_duv(1, 1);
    D = du1 * dv2 - du2 * dv1;
    dwudx.row(i) << (dv1 - dv2) / D, dv1 / D, -dv2 / D;
    dwvdx.row(i) << (du1 - du2) / D, -du1 / D, du2 / D;
  }
}

void dwdxi(Eigen::Vector3d &dwudxi, Eigen::Vector3d &dwvdxi, Eigen::MatrixXd wuv, double dwudx, double dwvdx) {
  // we use precomputed to form the appropriate matrices
  Eigen::Matrix3d ucomp, vcomp;
  ucomp = dwudx * Eigen::Matrix3d::Identity();
  vcomp = dwvdx * Eigen::Matrix3d::Identity();
  dwudxi = ucomp * wuv.col(0);
  dwvdxi = vcomp * wuv.col(1);
}
