#ifndef COMPUTE_ANGLE_DEFECT_HEADER_FILE
#define COMPUTE_ANGLE_DEFECT_HEADER_FILE

#include <Eigen/Dense>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Eigen::VectorXd compute_angle_defect(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::VectorXi& boundVMask){
    
    //Stub
    Eigen::VectorXd G = Eigen::VectorXd::Zero(V.rows());
    //TODO
    return G;
}


#endif
