#include "w_matrix.h"

void precompute_duv(Eigen::MatrixXd &global_duv, Eigen::MatrixXi F, Eigen::MatrixXd UV){
  global_duv.resize(2 * F.rows(), 2);
  for(int i = 0; i < F.rows(); i++){
    Eigen::RowVector3i tri = F.row(i);
    // we compute the duv we use to compute wuv here
    global_duv.block<2,2>(i * 2, 0) << UV(i * 3 + tri(1), 0) - UV(i * 3 + tri(0), 0), UV(i * 3 + tri(2), 0) - UV(i * 3 + tri(0), 0),
    UV(i * 3 + tri(1),1) - UV(i * 3 + tri(1), 1), UV(i * 3 + tri(2), 1) - UV(i * 3 + tri(0), 1);
  }
}

void wuv(Eigen::MatrixXd &wuv, Eigen::MatrixXi F, Eigen::MatrixXd V, int i, Eigen::Matrix2d duv){
  Eigen::RowVector3i tri = F.row(i);
  Eigen::MatrixXd dx(3,2);
  // construct dx
  dx.col(0) << V.row(tri(1)) - V.row(tri(0));
  dx.col(1) << V.row(tri(2)) - V.row(tri(0));
  // we calculate [Wu Wv] = dx * duv^-1
  // should be 3 x 2
  wuv = dx * duv.inverse();
}
