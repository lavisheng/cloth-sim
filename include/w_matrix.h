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
