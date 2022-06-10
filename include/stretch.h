#include <Eigen/Dense>

void stretch_cond(Eigen::Vector2d &cond, Eigen::MatrixXd wuv, Eigen::Vector2d b, Eigen::Vector2d uv, Eigen::RowVector3d coords);

void stretch_precompute(Eigen::MatrixXd &stretch_dwudx, Eigen::MatrixXd &stretch_dwvdx, Eigen::MatrixXi F, Eigen::MatrixXd UV);

void stretch_dcdxi(Eigen::MatrixXd & dcdxi, int i0, int i1, Eigen::MatrixXd wuv, Eigen::MatrixXd stretch_dwudx, Eigen::MatrixXd stretch_dwvdx);
