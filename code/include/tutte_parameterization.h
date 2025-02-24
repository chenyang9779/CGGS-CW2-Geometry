#ifndef COMPUTE_TUTTE_PARAMETERIZATION_HEADER_FILE
#define COMPUTE_TUTTE_PARAMETERIZATION_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "slice_columns_sparse.h"
#include "set_diff.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Eigen::MatrixXd compute_boundary_embedding(const Eigen::MatrixXd& V,
                                           const Eigen::VectorXi& boundVertices,
                                           const double r){
    
    //TODO
    return Eigen::MatrixXd::Zero(boundVertices.size(),2);
}

Eigen::MatrixXd compute_tutte_embedding(const Eigen::VectorXi& boundVertices,
                                        const Eigen::MatrixXd& UVBound,
                                        const Eigen::SparseMatrix<double>& d0,
                                        const Eigen::SparseMatrix<double>& W){
    
    //TODO
    return Eigen::MatrixXd::Zero(d0.cols(),2);
}


#endif
