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

void MassSpringsModel::force(Eigen::VectorXd &f, const Eigen::VectorXd q) {
    f = Eigen::VectorXd::Zero(q.rows());
    for (int e = 0; e < E.rows(); e++) {
        int v0 = E(e, 0);
        int v1 = E(e, 1);
        Eigen::Vector3d q0 = q.segment(E(e, 0), 3);
        Eigen::Vector3d q1 = q.segment(E(e, 1), 3);
        Eigen::Vector3d f0 = -k * ((q1 - q0).norm() - l0(e)) * (q1 - q0).normalized();
        f.segment<3>(v0 * 3) = f0;
        f.segment<3>(v1 * 3) = -f0;
    }
}

void MassSpringsModel::stiffness(Eigen::SparseMatrix<double> &K, const Eigen::VectorXd q) {
    std::list<Eigen::Triplet<double>> tl;
    for (int e = 0; e < E.rows(); e++) {
        int v0 = E(e, 0);
        int v1 = E(e, 1);
        Eigen::Vector3d q0 = q.segment(E(e, 0), 3);
        Eigen::Vector3d q1 = q.segment(E(e, 1), 3);
        double l = l0(e);
        double s1 = l * (q0(0) - q1(0)) * (q0(1) - q1(1));
        double s2 = l * (q0(0) - q1(0)) * (q0(2) - q1(2));
        double s3 = l * (q0(1) - q1(1)) * (q0(2) - q1(2));
        double t1 = l * (pow(q0(1) - q1(1), 2.) + pow(q0(2) - q1(2), 2.));
        double t2 = l * (pow(q0(0) - q1(0), 2.) + pow(q0(2) - q1(2), 2.));
        double t3 = l * (pow(q0(0) - q1(0), 2.) + pow(q0(1) - q1(1), 2.));
        double length_cubed = pow((q1 - q0).norm(), 3.);

        tl.emplace_back(v0 + 0, v0 + 0, -k * (-length_cubed + t1) / length_cubed);
        tl.emplace_back(v0 + 0, v0 + 1, k * s1 / length_cubed);
        tl.emplace_back(v0 + 0, v0 + 2, k * s2 / length_cubed);
        tl.emplace_back(v0 + 0, v1 + 0, k * (-length_cubed + t1) / length_cubed);
        tl.emplace_back(v0 + 0, v1 + 1, -k * s1 / length_cubed);
        tl.emplace_back(v0 + 0, v1 + 2, -k * s2 / length_cubed);
        tl.emplace_back(v0 + 1, v0 + 0, k * s1 / length_cubed);
        tl.emplace_back(v0 + 1, v0 + 1, -k * (-length_cubed + t2) / length_cubed);
        tl.emplace_back(v0 + 1, v0 + 2, k * s3 / length_cubed);
        tl.emplace_back(v0 + 1, v1 + 0, -k * s1 / length_cubed);
        tl.emplace_back(v0 + 1, v1 + 1, k * (-length_cubed + t2) / length_cubed);
        tl.emplace_back(v0 + 1, v1 + 2, -k * s3 / length_cubed);
        tl.emplace_back(v0 + 2, v0 + 0, k * s2 / length_cubed);
        tl.emplace_back(v0 + 2, v0 + 1, k * s3 / length_cubed);
        tl.emplace_back(v0 + 2, v0 + 2, -k * (-length_cubed + t3) / length_cubed);
        tl.emplace_back(v0 + 2, v1 + 0, -k * s2 / length_cubed);
        tl.emplace_back(v0 + 2, v1 + 1, -k * s3 / length_cubed);
        tl.emplace_back(v0 + 2, v1 + 2, k * (-length_cubed + t3) / length_cubed);
        tl.emplace_back(v1 + 0, v0 + 0, k * (-length_cubed + t1) / length_cubed);
        tl.emplace_back(v1 + 0, v0 + 1, -k * s1 / length_cubed);
        tl.emplace_back(v1 + 0, v0 + 2, -k * s2 / length_cubed);
        tl.emplace_back(v1 + 0, v1 + 0, -k * (-length_cubed + t1) / length_cubed);
        tl.emplace_back(v1 + 0, v1 + 1, k * s1 / length_cubed);
        tl.emplace_back(v1 + 0, v1 + 2, k * s2 / length_cubed);
        tl.emplace_back(v1 + 1, v0 + 0, -k * s1 / length_cubed);
        tl.emplace_back(v1 + 1, v0 + 1, k * (-length_cubed + t2) / length_cubed);
        tl.emplace_back(v1 + 1, v0 + 2, -k * s3 / length_cubed);
        tl.emplace_back(v1 + 1, v1 + 0, k * s1 / length_cubed);
        tl.emplace_back(v1 + 1, v1 + 1, -k * (-length_cubed + t2) / length_cubed);
        tl.emplace_back(v1 + 1, v1 + 2, k * s3 / length_cubed);
        tl.emplace_back(v1 + 2, v0 + 0, -k * s2 / length_cubed);
        tl.emplace_back(v1 + 2, v0 + 1, -k * s3 / length_cubed);
        tl.emplace_back(v1 + 2, v0 + 2, k * (-length_cubed + t3) / length_cubed);
        tl.emplace_back(v1 + 2, v1 + 0, k * s2 / length_cubed);
        tl.emplace_back(v1 + 2, v1 + 1, k * s3 / length_cubed);
        tl.emplace_back(v1 + 2, v1 + 2, -k * (-length_cubed + t3) / length_cubed);
    }
    K.setZero();
    K.setFromTriplets(tl.begin(), tl.end());
}
