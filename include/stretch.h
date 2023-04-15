#include <Eigen/Dense>

/**
 * Calculates C(x) for the stretch condition
 * Inputs:
 *   wuv - wuv matrix
 *   b - scalars controlling desired stretch amount in uv respectively
 *   a - area of triangle in (u,v) space (precomputable
 */
void
stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, double area, Eigen::RowVector3d coords);

/**
 * Precomputes matrix values that are used to calculate dcdxi
 * Inputs:
 *   F - #F by 3 matrix storing indices to vertices
 *   global_duv - #F x 2 by 2 stack of duvs that were precomputed 
 * Out:
 *   stretch_dwudx - #F by 3 list of dwudx constants
 *   stretch_dwvdx - #F by 3 list of dwvdx constants
 */
void stretch_precompute(Eigen::MatrixXd &stretch_dwudx, Eigen::MatrixXd &stretch_dwvdx, Eigen::MatrixXi F,
                        Eigen::MatrixXd global_duv);

/**
 * Computes dcdx for the mesh
 * Inputs:
 *   F - #F by 3 list of triangle faces
 *   V_size - #V 
 *   a - #F list of triangle areas (precomputed)
 *   wuv - #F x 2 by 2 stacks of wuv matrices
 *   dwudx - #F x 3 stacks of dwudx per triangle
 *   dwvdx - #F x 3 stacks of dwvdx per triangle
 * Out:
 *   dcudxi - derivative of cu respect to xi
 *   dcvdxi - derivative of cv respect to xi
 */
void stretch_dcdxi(Eigen::VectorXd &dcudxi, Eigen::VectorXd &dcvdxi, Eigen::MatrixXd F, int V_size, Eigen::VectorXd a,
                   Eigen::MatrixXd wuv, Eigen::VectorXd dwudx, Eigen::VectorXd dwvdx);

/**
 * Computes the d2c/dxixj for a pair of vertices in a triangle
 * Inputs:
 *   a - area of triangle in (u,v) space, precomputable
 *   wuv - wuv matrix
 *   dwudx_prod - dwudxi * dwudxj for calculating the second derivative
 *   dwvdx_prod - dwvdxi * dwvdxj for calculating second derivative
 * Out:
 *   dwudxixj - second derivative of wu respect to xi, xj
 *   dwvdxixj - second derivative of wv respect to xi, xj
 */
void
stretch_d2cdxixj(Eigen::Matrix3d &dwudxixj, Eigen::Matrix3d &dwvdxixj, double a, Eigen::MatrixXd wuv, double dwudx_prod,
                 double dwvdx_prod);
