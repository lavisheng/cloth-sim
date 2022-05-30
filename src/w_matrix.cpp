#include "w_matrix.h"

void wuv(Eigen::MatrixXd &wuv, int i, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd UV){
  //double dx1,dx2, du1, du2, dv1, dv2;
  Eigen::MatrixXd dx(3,2);
  Eigen::MatrixXd duv(2,2);
  Eigen::RowVector3d tri = F.row(i);
  // construct duv
  duv << UV(i * 3 + tri(1), 0) - UV(i * 3 + tri(0), 0), UV(i * 3 + tri(2), 0) - UV(i * 3 + tri(0), 0),
    UV(i * 3 + tri(1),1) - UV(i * 3 + tri(1), 1), UV(i * 3 + tri(2), 1) - UV(i * 3 + tri(0), 1);
  // construct dx
  dx.col(0) << V.row(tri(1)) - V.row(tri(0), 0);
  dx.col(1) << V.row(tri(2)) - V.row(tri(0), 0);
  wuv = dx * duv.inverse();
}
