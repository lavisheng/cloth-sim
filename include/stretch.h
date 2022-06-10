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
 *   F - #F x 3 matrix storing indices to vertices
 *   UV - #UV by 2 list of UV coordinates
 */
void stretch_precompute(Eigen::MatrixXd &stretch_dwudx, Eigen::MatrixXd &stretch_dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd UV);

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
