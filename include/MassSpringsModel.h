#include "Model.h"

class MassSpringsModel : public Model {
public:
    MassSpringsModel(Eigen::MatrixXd &V, Eigen::MatrixXi &F, double k);

    bool init_precompute();

    void potential_energy(double &energy, const Eigen::VectorXd &q);

    void force(Eigen::VectorXd &f, const Eigen::VectorXd &q);

    void stiffness(Eigen::SparseMatrix<double> &K, const Eigen::VectorXd &q);

private:
    double k;  // stiffness
    Eigen::MatrixXi E;  // edges
    Eigen::VectorXd l0;  // edge lengths (precomputed)
};
