#ifndef COMPUTE_MEAN_CURVATURE_FLOW_HEADER_FILE
#define COMPUTE_MEAN_CURVATURE_FLOW_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "compute_areas_normals.h"

inline double compute_total_area(
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &V

)
{
    Eigen::VectorXd fAreas;
    Eigen::MatrixXd fNormals;

    compute_areas_normals(V, F, fAreas, fNormals);
    return fAreas.sum();
}

void mean_curvature_flow(const Eigen::MatrixXi &F,
                         const Eigen::SparseMatrix<double> &L,
                         const double timeStep,
                         const Eigen::SparseMatrix<double> &M,
                         const Eigen::SparseMatrix<double> &MInv,
                         const Eigen::VectorXi &boundVMask,
                         const bool isExplicit,
                         Eigen::MatrixXd &currV)
{

    using namespace Eigen;
    using namespace std;
    // TODO
    MatrixXd oldV = currV;

    double oldArea = compute_total_area(F, currV);

    if (isExplicit)
    {
        MatrixXd Lv = L * currV;
        MatrixXd update = MInv * Lv;
        currV -= timeStep * update;
    }
    else
    {
        SparseMatrix<double> A = M + timeStep * L;
        MatrixXd b = M * currV;

        SimplicialLDLT<SparseMatrix<double>> solver;
        solver.compute(A);
        currV = solver.solve(b);
    }

    for (int i = 0; i < boundVMask.size(); ++i)
    {
        if (boundVMask(i) == 1)
        {
            currV.row(i) = oldV.row(i);
        }
    }

    bool hasBoundary = (boundVMask.array() == 1).any();
    if (!hasBoundary)
    {

        double newArea = compute_total_area(F, currV);
        double scale = std::sqrt(oldArea / newArea);
        currV *= scale;

        RowVectorXd center = currV.colwise().mean();
        currV.rowwise() -= center;
    }
}

#endif
