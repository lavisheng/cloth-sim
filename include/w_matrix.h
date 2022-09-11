#include <Eigen/Dense>
/**
 * Precomputes the duv we use in wuv and other places
 * Inputs:
 *   F - #F by 3 list of indices into V
 *   UV - #UV by 2 list of UV coordinates
 * Out:
 *   global_duv - a #F*2x2 matrix with all the duvs for each triangle stacked
 */
void precompute_duv(Eigen::MatrixXd &global_duv, Eigen::MatrixXi F, Eigen::MatrixXd UV);
/**
 * Calculates the W_{uv} matrix which describes the stretch and squash across a cloth's 
 * uv coordinates.
 * Inputs:
 *   F - #F by 3 list of triangle indices
 *   V - #V by 3 list of vertex positions
 *   i - index into F
 *   duv - 2 by 2 matrix used to calculate wuv, precomputed
 * Out:
 *   wuv - the w_{uv} matrix for a given triangle
 */
  void wuv(Eigen::MatrixXd &wuv, Eigen::MatrixXi F, Eigen::MatrixXd V, int i, Eigen::MatrixXd duv);

/**
 * Precomputes matrix values used to calculate dwu/dxi and dwvdxi
 * Inputs:
 *   F - #F by 3 matrix storing indices to vertices
 *   duv - #F x 2 by 2 stack of duvs that were precomputed
 * Out:
 *   dwudx - #F by 3 list of dwudx constants
 *   dwvdx - #F by 3 list of dwvdx constants
 */
void precompute_dwdx(Eigen::MatrixXd &dwudx, Eigen::MatrixXd dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd duv);

/**
 * Computes dwdxi for a vertex in a triangle
 * Inputs:
 *   wuv - wuv matrix
 *   dwudx - precomputed dwudx matrix value in dwudx
 *   dwvdx - precomputed dwvdx matrix value in dwvdx
 * Out:
 *   dwudxi - derivative of wu respect to xi
 *   dwvdxi - derivative of wv respect to xi
 */
void dxi(Eigen::Vector3d &dwudxi, Eigen::Vector3d &dwvdxi, Eigen::MatrixXd wuv, double dwudx, double dwvdx);
