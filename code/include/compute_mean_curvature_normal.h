#ifndef COMPUTE_MEAN_CURVATURE_NORMAL_HEADER_FILE
#define COMPUTE_MEAN_CURVATURE_NORMAL_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>

void compute_mean_curvature_normal(const Eigen::MatrixXd& V,
                                   const Eigen::MatrixXi& F,
                                   const Eigen::SparseMatrix<double>& L,
                                   const Eigen::VectorXd& vorAreas,
                                   Eigen::MatrixXd& Hn,
                                   Eigen::VectorXd& H){
    using namespace Eigen;
    using namespace std;
    Hn = MatrixXd::Zero(V.rows(), 3);
    H = VectorXd::Zero(V.rows());
}


#endif
