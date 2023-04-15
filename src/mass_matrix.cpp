#include "mass_matrix.h"

void mass_matrix(Eigen::SparseMatrix<double> &M, Eigen::MatrixXd &V, double m0){
  M = Eigen::SparseMatrix<double>(V.rows() * 3, V.rows() * 3);
  M.setIdentity();
  M *= m0;
}
