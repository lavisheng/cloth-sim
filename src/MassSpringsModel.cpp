#include "MassSpringsModel.h"
#include <igl/edges.h>
#include <igl/edge_lengths.h>

MassSpringsModel::MassSpringsModel(Eigen::MatrixXd &V, Eigen::MatrixXi &F, double k) : Model(V, F), k(k) {}

bool MassSpringsModel::init_precompute() {
    igl::edges(F, E);
    igl::edge_lengths(V, F, l0);
    return true;
}

void MassSpringsModel::potential_energy(double &energy, const Eigen::VectorXd q) {
    double energy_temp = 0;
    for (int e = 0; e < E.rows(); e++) {
        Eigen::Vector3d q0 = q.segment(E(e, 0), 3);
        Eigen::Vector3d q1 = q.segment(E(e, 1), 3);
        energy_temp += 0.5 * k * pow((q1 - q0).norm() - l0(e), 2);
    }
    energy = energy_temp;
}

// TODO
void MassSpringsModel::force(Eigen::VectorXd &f, const Eigen::VectorXd q) {

}

// TODO
void MassSpringsModel::stiffness(Eigen::SparseMatrix<double> &K, const Eigen::VectorXd q) {

}
