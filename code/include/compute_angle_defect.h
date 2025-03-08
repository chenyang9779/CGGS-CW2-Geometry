#ifndef COMPUTE_ANGLE_DEFECT_HEADER_FILE
#define COMPUTE_ANGLE_DEFECT_HEADER_FILE

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Eigen::VectorXd compute_angle_defect(const Eigen::MatrixXd &V,
                                     const Eigen::MatrixXi &F,
                                     const Eigen::VectorXi &boundVMask)
{

    // Stub
    Eigen::VectorXd angleSums = Eigen::VectorXd::Zero(V.rows());
    Eigen::VectorXd G = Eigen::VectorXd::Zero(V.rows());
    // TODO

    int nF = F.rows();
    for (int f = 0; f < nF; ++f) {
        int i0 = F(f, 0);
        int i1 = F(f, 1);
        int i2 = F(f, 2);
        Eigen::RowVector3d v0 = V.row(i0);
        Eigen::RowVector3d v1 = V.row(i1);
        Eigen::RowVector3d v2 = V.row(i2);
        

        {
            Eigen::RowVector3d a = v1 - v0;
            Eigen::RowVector3d b = v2 - v0;
            double cosAngle = a.dot(b) / (a.norm() * b.norm());
            double angle = std::acos(cosAngle);
            angleSums(i0) += angle;
        }

        {
            Eigen::RowVector3d a = v0 - v1;
            Eigen::RowVector3d b = v2 - v1;
            double cosAngle = a.dot(b) / (a.norm() * b.norm());
            double angle = std::acos(cosAngle);
            angleSums(i1) += angle;
        }

        {
            Eigen::RowVector3d a = v0 - v2;
            Eigen::RowVector3d b = v1 - v2;
            double cosAngle = a.dot(b) / (a.norm() * b.norm());
            double angle = std::acos(cosAngle);
            angleSums(i2) += angle;
        }
    }

    for (int i = 0; i < V.rows(); ++i) {
        if (boundVMask(i) == 1)
            G(i) = M_PI - angleSums(i);
        else
            G(i) = 2 * M_PI - angleSums(i);
    }
    return G;
}

#endif
