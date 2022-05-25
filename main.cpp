#include <igl/readMESH.h>
#include <string>

void setup(std::string file, Eigen::VectorXd &V, Eigen::VectorXd &F){
  igl::readOBJ(file, V, F);
}
