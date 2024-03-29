#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <stdio.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS>
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                                    const Eigen::SparseMatrix<double> &mass, FORCE &force, STIFFNESS &stiffness,
                                    Eigen::VectorXd &tmp_force, Eigen::SparseMatrix<double> &tmp_stiffness,
                                    Eigen::SparseMatrix<double> fp) {
    // eqn is A qdot = b, where A = M - dt^2 K and b = Mqdot + dt f(q)
    // solve for qdot by solving sparse system
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    stiffness(tmp_stiffness, q, qdot);
    solver.compute(fp * (mass - dt * dt * tmp_stiffness) * fp.transpose());
    force(tmp_force, q, qdot);
    qdot = fp.transpose() * solver.solve(fp * (mass * qdot + dt * tmp_force));
    // update q using updated qdot
    q = q + dt * qdot;
}
