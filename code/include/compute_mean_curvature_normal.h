#ifndef COMPUTE_MEAN_CURVATURE_NORMAL_HEADER_FILE
#define COMPUTE_MEAN_CURVATURE_NORMAL_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>

void compute_mean_curvature_normal(const Eigen::MatrixXd &V,
                                   const Eigen::MatrixXi &F,
                                   const Eigen::SparseMatrix<double> &L,
                                   const Eigen::VectorXd &vorAreas,
                                   Eigen::MatrixXd &Hn,
                                   Eigen::VectorXd &H)
{
    using namespace Eigen;
    int nV = V.rows();
    int nF = F.rows();

    MatrixXd Lv = L * V;

    Hn.resize(nV, 3);
    for (int i = 0; i < nV; ++i)
    {
        double area = vorAreas(i);
        Hn.row(i) = Lv.row(i) / (2.0 * area);
    }

    MatrixXd faceNormals(nF, 3);
    for (int f = 0; f < nF; ++f)
    {
        int i0 = F(f, 0), i1 = F(f, 1), i2 = F(f, 2);
        Vector3d v0 = V.row(i0);
        Vector3d v1 = V.row(i1);
        Vector3d v2 = V.row(i2);
        Vector3d n = (v1 - v0).cross(v2 - v0);
        double norm = n.norm();
        faceNormals.row(f) = n;
    }

    MatrixXd vertexNormals = MatrixXd::Zero(nV, 3);
    for (int f = 0; f < nF; ++f)
    {
        for (int j = 0; j < 3; ++j)
        {
            int v = F(f, j);
            vertexNormals.row(v) += faceNormals.row(f);
        }
    }

    for (int i = 0; i < nV; ++i)
    {
        double norm = vertexNormals.row(i).norm();
    }

    H.resize(nV);
    for (int i = 0; i < nV; ++i)
    {
        double Hn_norm = Hn.row(i).norm();
        double dot = Hn.row(i).dot(vertexNormals.row(i));
        H(i) = (dot < 0) ? -Hn_norm : Hn_norm;
    }
}

#endif
