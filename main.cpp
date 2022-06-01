#include <igl/readOBJ.h>
#include <string>
#include "igl/opengl/glfw/Viewer.h"

void setup(std::string file, Eigen::VectorXd &V, Eigen::MatrixXi &F) {
    igl::readOBJ(file, V, F);
}

int main(int argc, char **argv) {
    Eigen::VectorXd V;
    Eigen::MatrixXi F;
//    auto filename = argc > 1 ? argv[1] : "../data/square_cloth.obj";
//    setup(filename, V, F);
//
//    igl::opengl::glfw::Viewer viewer;
//    viewer.data().set_mesh(V, F);
//    viewer.launch();
    return 0;
}