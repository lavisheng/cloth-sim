#include "shear.h"

double shear_cond(Eigen::MatrixXd wuv, double area){
  // shear condition is a * wu^T * wv
  return area * wuv.col(0).transpose() * wuv.col(1);
}
