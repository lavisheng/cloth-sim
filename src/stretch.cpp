#include "stretch.h"
#Include "w_matrix.h"
void stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, double a, Eigen::Vector2d uv, Eigen::RowVector3d coords){
  cond(0) = a * (wuv.col(0).norm() - b(0));
  cond(1) = a * (wuv.col(1).norm() - b(1));
}

void stretch_dcdxi(Eigen::MatrixXd &dcdxi, int i, Eigen::MatrixXd wuv){
  dcdxi.resize(3,2);
  // we can precompute the matrices we need to calculate, and then choose and store them
  // in ucomp, vcomp
  // TODO: do the precomputation steps.
  Eigen::MatrixXd ucomp, vcomp;
  dcdxi.col(0) << ucomp * wuv.col(0);
  dcdxi.col(1) << vcomp * wuv.col(1);
}
