#include <Eigen/Dense>

/**
   Calculates the C(x) for the shear condition
 * Inputs:
 *   wuv - wuv matrix
 *   a - area of triangle in (u, v) space (precomputable)
 * Out:
 *   cond - the condition value
 */
double shear_cond( Eigen::MatrixXd wuv, double a);
