#include <Eigen/Dense>
/**
 * Calculates C(x) for the stretch condition
 * Inputs:
 *   wuv - wuv matrix
 *   b - scalars controlling desired stretch amount in uv respectively
 *   a - area of triangle in (u,v) space (precomputable
 */
void stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, double a, Eigen::RowVector3d coords);

/**
 * Precomputes matrix values that are used to calculate dcdxi
 * Inputs:
 *   F - #F by 3 matrix storing indices to vertices
 *   global_duv - #F x 2 by 2 stack of duvs that were precomputed 
 * Out:
 *   stretch_dwudx - #F by 3 list of dwudx constants
 *   stretch_dwvdx - #F by 3 list of dwvdx constants
 */
void stretch_precompute(Eigen::MatrixXd &stretch_dwudx, Eigen::MatrixXd &stretch_dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd global_duv);

/**
 * Computes dcdxi for a vertex in a triangle
 * Inputs:
 *   wuv - wuv matrix
 *   dwudx - precomputed dwudx matrix value in dwudx
 *   dwvdx - precomputed dwvdx matrix value in dwvdx
 * Out:
 *   dcdxi - derivative of c
 */
void stretch_dcdxi(Eigen::MatrixXd & dcdxi, Eigen::MatrixXd wuv, double dwudx, double dwvdx);

/**
 * Computes the d2c/dxixj for a pair of vertices in a triangle
 * Inputs:
 *   a - area of triangle in (u,v) space, precomputable
 *   wuv - wuv matrix
 *   dwudx_prod - dwudxi * dwudxj for calculating the second derivative
 *   dwvdx_prod - dwvdxi * dwvdxj for calculating second derivative
 * Out:
 *   d2cdxixj - second derivative of c
 */
void stretch_d2cdxixj(Eigen::MatrixXd &d2cdxixj, double a, Eigen::MatrixXd wuv, double dwudx_prod, double dwvdx_prod);
