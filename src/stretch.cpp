#include "stretch.h"

void stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, double a){
  cond(0) = a * (wuv.col(0).norm() - b(0));
  cond(1) = a * (wuv.col(1).norm() - b(1));
}

void stretch_dcdxi(Eigen::VectorXd &dcudxi, Eigen::VectorXd &dcvdxi, Eigen::MatrixXd F, int V_size, Eigen::VectorXd a, Eigen::MatrixXd wuv, Eigen::VectorXd dwudx, Eigen::VectorXd dwvdx){  
  dcudxi = Eigen::VectorXd::Zero(V_size * 3);
  dcvdxi = Eigen::VectorXd::Zero(V_size * 3);
  // iterate through the mesh and all the triangles
  for(int i = 0; i < F.rows(); i++){
    for(int f_i = 0; f_i < 3; f_i++){
      int v_index = F(i, f_i);
      // calculates a * dwudx * wu
      dcudxi.segment(v_index * 3, (v_index + 1) * 3) = a[i] * dwudx.segment(v_index*3, (v_index + 1) * 3) * wuv.col(0).segment(v_index * 3, (v_index + 1) * 3);
      // calculates a * dwvdx * wv
      dcvdxi.segment(v_index * 3, (v_index + 1) * 3) = a[i] * dwvdx.segment(v_index*3, (v_index + 1) * 3) * wuv.col(1).segment(v_index * 3, (v_index + 1) * 3);
    }
  }
}

void stretch_d2cdxixj(Eigen::Matrix3d &d2cudxixj, Eigen::Matrix3d &d2cvdxixj, double a, Eigen::MatrixXd wuv, double dwudx_prod, double dwvdx_prod){
  //dwudxixj.resize(3,3);
  //dwvdxixj.resize(3,3);
  // a / ||w_u|| dw/dxi dw/dxj (I - w_u * w_u^T) and vice versa for v
  d2cudxixj = a / wuv.col(0).norm() * dwudx_prod * (Eigen::Matrix3d::Identity() - wuv.col(0) * wuv.col(0).transpose());
  d2cvdxixj = a / wuv.col(0).norm() * dwvdx_prod * (Eigen::Matrix3d::Identity() - wuv.col(1) * wuv.col(1).transpose());
  
}
