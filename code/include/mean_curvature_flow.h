#ifndef COMPUTE_MEAN_CURVATURE_FLOW_HEADER_FILE
#define COMPUTE_MEAN_CURVATURE_FLOW_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "compute_areas_normals.h"

void mean_curvature_flow(const Eigen::MatrixXi& F,
                         const Eigen::SparseMatrix<double>& L,
                         const double timeStep,
                         const Eigen::SparseMatrix<double>& M,
                         const Eigen::SparseMatrix<double>& MInv,
                         const Eigen::VectorXi& boundVMask,
                         const bool isExplicit,
                         Eigen::MatrixXd& currV){
    
    using namespace Eigen;
    using namespace std;
    //TODO
}


#endif
