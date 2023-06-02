#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifndef CLOTH_SIM_MODEL_H
#define CLOTH_SIM_MODEL_H

class Model {
public:
    Model(Eigen::MatrixXd &V, Eigen::MatrixXi &F) : V(V), F(F) {};

    virtual bool init_precompute() = 0;

    virtual void potential_energy(double &energy, const Eigen::VectorXd &q) = 0;

    virtual void force(Eigen::VectorXd &f, const Eigen::VectorXd &q) = 0;

    virtual void stiffness(Eigen::SparseMatrix<double> &K, const Eigen::VectorXd &q) = 0;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
};

#endif
