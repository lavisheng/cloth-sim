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
 * Computes dcdxi for a triangle
 * Inputs:
 *   i0 - triangle id
 *   i1 - vertex in triangle id
 *   wuv - wuv matrix
 *   stretch_dwudx - precomputed dwudx matrix
 *   stretch_dwvdx - precomputed dwvdx
 */
void stretch_dcdxi(Eigen::MatrixXd & dcdxi, int i0, int i1, Eigen::MatrixXd wuv, Eigen::MatrixXd stretch_dwudx, Eigen::MatrixXd stretch_dwvdx);
