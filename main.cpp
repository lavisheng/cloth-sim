#include <igl/readOBJ.h>
#include <string>
#include "igl/opengl/glfw/Viewer.h"
#include "Model.h"
#include "MassSpringsModel.h"

void setup(std::string file, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    igl::readOBJ(file, V, F);
}

void animation_loop(Model *model, Eigen::VectorXd q, Eigen::VectorXd qdot, double dt) {
    Eigen::VectorXd f;
    Eigen::SparseMatrix<double> K;
    // TODO: fix animation loop
    model->force(f, q);
    model->stiffness(K, q);
    // linearly_implicit_euler(...)
}

int main(int argc, char **argv) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    auto filename = argc > 1 ? argv[1] : "../data/square_cloth.obj";
    setup(filename, V, F);

    Model *model = new MassSpringsModel(V, F, 1);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();
    return 0;
}
