#include <Eigen/Dense>

/**
 * Calculates the W_{uv} matrix which describes the stretch and squash across a cloth's 
 * uv coordinates.
 * Inputs:
 *   i - triangle index
 *   V - #V by 3 list of vertex positions
 *   F - #F by 3 list of indieces into V
 *   UV - #UV by 2 list of UV coordinates
 */
void wuv(Eigen::MatrixXd &wuv, int i, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd UV);
