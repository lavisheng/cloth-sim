#include <igl/readOBJ.h>
#include <string>

void setup(std::string file, Eigen::VectorXd &V, Eigen::MatrixXi &F) {
    igl::readOBJ(file, V, F);
}
