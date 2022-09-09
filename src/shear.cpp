#include "shear.h"
double shear_cond(Eigen::MatrixXd wuv, double a){
  // shear condition is a * wu^T * wv
  return a * wuv.col(0).transpose() * wuv.col(1);
}
