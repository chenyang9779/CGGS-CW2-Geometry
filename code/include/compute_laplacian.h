#ifndef COMPUTE_LAPLACIAN_HEADER_FILE
#define COMPUTE_LAPLACIAN_HEADER_FILE

#include <Eigen/Dense>

void compute_laplacian(const Eigen::MatrixXd& V,
                       const Eigen::MatrixXi& F,
                       const Eigen::MatrixXi& E,
                       const Eigen::MatrixXi& EF,
                       const Eigen::VectorXi& boundEMask,
                       Eigen::SparseMatrix<double>& d0,
                       Eigen::SparseMatrix<double>& W,
                       Eigen::VectorXd& vorAreas){
    
    using namespace Eigen;
    using namespace std;
    d0.resize(E.rows(), V.rows());
    W.resize(E.rows(), E.rows());
    vorAreas = VectorXd::Ones(V.rows());
   
    //TODO
}



#endif
