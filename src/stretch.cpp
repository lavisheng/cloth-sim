#include "stretch.h"
#include "w_matrix.h"

void stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, double a){
  cond(0) = a * (wuv.col(0).norm() - b(0));
  cond(1) = a * (wuv.col(1).norm() - b(1));
}

void stretch_precompute(Eigen::MatrixXd &stretch_dwudx, Eigen::MatrixXd &stretch_dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd global_duv){
  stretch_dwudx.resize(F.rows(), 3);
  stretch_dwvdx.resize(F.rows(), 3);
  Eigen::MatrixXd duv(2,2);
  double D;
  for(int i = 0; i < F.rows(); i++){
    // basically use what we did in w_matrix to get duv
    duv = global_duv.block<2,2>(i * 2, 0);
    // the denominator du1 * dv2 - du2 * dv1
    D = duv(0,0) * duv(1,1) - duv(0,1) * duv(1,0);
    stretch_dwudx.row(i) << (duv(1,0) - duv(1,1))/ D, duv(1,0)/ D, -duv(1,1) / D;
    stretch_dwvdx.row(i) << (duv(0,0) - duv(0,1))/ D , -duv(0,0)/ D, duv(0,1) / D;
  }
}
void stretch_dcdxi(Eigen::MatrixXd &dcdxi,  Eigen::MatrixXd wuv, double dwudx, double dwvdx){
  dcdxi.resize(3,2);
  // we use precomputed to form the appropriate matrices
  Eigen::Matrix3d ucomp, vcomp;
  ucomp = dwudx * Eigen::Matrix3d::Identity();
  vcomp = dwvdx * Eigen::Matrix3d::Identity();
  dcdxi.col(0) << ucomp * wuv.col(0);
  dcdxi.col(1) << vcomp * wuv.col(1);
}

void stretch_d2cdxixj(Eigen::MatrixXd &d2cdxixj, double a, Eigen::MatrixXd wuv, double dwudx_prod, double dwvdx_prod){
  d2cdxixj.resize(3,2);
  // a / ||w_u|| dw/dxi dw/dxj (I - w_u * w_u^T) and vice versa for v
  d2cdxixj.col(0) << a / wuv.col(0).norm() * dwudx_prod * (Eigen::Matrix3d::Identity() - wuv.col(0) * wuv.col(0).transpose());
  d2cdxixj.col(1) << a / wuv.col(0).norm() * dwvdx_prod * (Eigen::Matrix3d::Identity() - wuv.col(1) * wuv.col(1).transpose());
  
}
