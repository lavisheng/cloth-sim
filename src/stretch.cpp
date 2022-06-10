#include "stretch.h"
#include "w_matrix.h"

void stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, double a, Eigen::Vector2d uv, Eigen::RowVector3d coords){
  cond(0) = a * (wuv.col(0).norm() - b(0));
  cond(1) = a * (wuv.col(1).norm() - b(1));
}

void stretch_precompute(Eigen::MatrixXd &stretch_dwudx, Eigen::MatrixXd &stretch_dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd UV){
  stretch_dwudx.resize(F.rows(), 3);
  stretch_dwvdx.resize(F.rows(), 3);
  Eigen::RowVector3i tri;
  Eigen::MatrixXd duv(2,2);
  double D;
  for(int i = 0; i < F.rows(); i++){
    tri = F.row(i);
    duv << UV(i * 3 + tri(1), 0) - UV(i * 3 + tri(0), 0), UV(i * 3 + tri(2), 0) - UV(i * 3 + tri(0), 0),
      UV(i * 3 + tri(1),1) - UV(i * 3 + tri(1), 1), UV(i * 3 + tri(2), 1) - UV(i * 3 + tri(0), 1);
    D = duv(0,0) * duv(1,1) - duv(0,1) * duv(1,0);
    stretch_dwudx.row(i) << (duv(1,0) - duv(1,1))/ D, duv(1,0)/ D, duv(1,1) / D;
    stretch_dwvdx.row(i) << (duv(0,0) - duv(0,1))/ D , duv(0,0)/ D, duv(0,1) / D;
  }
}
void stretch_dcdxi(Eigen::MatrixXd &dcdxi, int i0, int i1, Eigen::MatrixXd wuv, Eigen::MatrixXd stretch_dwudx, Eigen::MatrixXd stretch_dwvdx){
  dcdxi.resize(3,2);
  // we can precompute the matrices we need to calculate, and then choose and store them
  // in ucomp, vcomp
  // TODO: do the precomputation steps.
  Eigen::Matrix3d ucomp, vcomp;
  ucomp = stretch_dwudx(i0, i1) * Eigen::Matrix3d::Identity();
  vcomp = stretch_dwvdx(i0,i1) * Eigen::Matrix3d::Identity();
  dcdxi.col(0) << ucomp * wuv.col(0);
  dcdxi.col(1) << vcomp * wuv.col(1);
}
