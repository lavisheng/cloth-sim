#include <igl/readOBJ.h>
#include <string>
#include "igl/opengl/glfw/Viewer.h"
#include "Model.h"
#include "MassSpringsModel.h"
#include "linearly_implicit_euler.h"


double t = 0;
double dt = 0.01;

Eigen::VectorXd flatten(const Eigen::MatrixXd &V) {
    long rows = V.rows(), cols = V.cols();
    Eigen::VectorXd out(rows * cols);
    for (long r = 0; r < rows; r++) {
        out.segment(r * cols, cols) = V.row(r);
    }
    return out;
}

Eigen::MatrixXd unflatten(const Eigen::VectorXd &v, int dim) {
    long rows = v.rows() / dim, cols = dim;
    Eigen::MatrixXd out(rows, cols);
    for (long r = 0; r < rows; r++) {
        out.row(r) = v.segment(r * cols, cols);
    }
    return out;
}

void setup(std::string file, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    igl::readOBJ(file, V, F);
}

void
animation_loop(Model *model, Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::SparseMatrix<double> &M, const double h,
               const Eigen::SparseMatrix<double> &fp) {
    Eigen::VectorXd f(q.rows());
    Eigen::SparseMatrix<double> K_temp(q.rows(), q.rows());
    // TODO: fix animation loop

    auto force = [&](Eigen::VectorXd &f,
                     Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot) {
        model->force(f, q);
        // add force of gravity
        Eigen::Vector3d g_single;
        g_single << 0, -9.8, 0;
        Eigen::VectorXd g = g_single.replicate(model->V.rows(), 1);
        f += M * g;
        f = f;
    };

    auto stiffness = [&](Eigen::SparseMatrix<double> &K,
                         Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot) {
        model->stiffness(K, q);
        K = K;
    };

    // time integrate
    linearly_implicit_euler(q, qdot, h, M, force, stiffness, f, K_temp, fp);
}

int main(int argc, char **argv) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    auto filename = argc > 1 ? argv[1] : "../data/square_cloth_rotated.obj";
    setup(filename, V, F);
    long n = V.rows();

    Eigen::VectorXd q = flatten(V);
    Eigen::VectorXd qdot = Eigen::VectorXd::Zero(n * 3);

    Eigen::SparseMatrix<double> M(n * 3, n * 3);
    M.setIdentity();

    Model *model = new MassSpringsModel(V, F, 1);
    model->init_precompute();

    std::vector<int> unfixed_point_indices;
    for (int v = 0; v < model->V.rows(); v++) {
        if (model->V(v, 1) < 0.433013) {
            unfixed_point_indices.emplace_back(v);
        }
    }
    std::list<Eigen::Triplet<double>> fixed_tl;
    for (int i = 0; i < unfixed_point_indices.size(); i++) {
        int v = unfixed_point_indices[i];
        fixed_tl.emplace_back(i * 3 + 0, v * 3 + 0, 1);
        fixed_tl.emplace_back(i * 3 + 1, v * 3 + 1, 1);
        fixed_tl.emplace_back(i * 3 + 2, v * 3 + 2, 1);
    }
    Eigen::SparseMatrix<double> fp(unfixed_point_indices.size() * 3, model->V.rows() * 3);
    fp.setFromTriplets(fixed_tl.begin(), fixed_tl.end());
    std::cout << fp.squaredNorm() / 3 << "\n";

    igl::opengl::glfw::Viewer viewer;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool {
        if (viewer.core().is_animating and t < 1) {
            animation_loop(model, q, qdot, M, dt, fp);
            t += dt;
            viewer.data(0).set_vertices(unflatten(q, 3));
        }
        return false;
    };
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod) {
        std::cout << key << "\n";
        switch (key) {
            case ' ':
                std::cout << "space" << "\n";
                viewer.core().is_animating = !viewer.core().is_animating;
                break;
            case 'c':
            case 'C':
                std::cout << "c" << "\n";
                break;
            default:
                break;
        }
        return false;
    };
    viewer.data().set_mesh(V, F);
    viewer.launch();
    return 0;
}
