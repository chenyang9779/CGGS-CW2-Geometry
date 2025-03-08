#ifndef COMPUTE_LAPLACIAN_HEADER_FILE
#define COMPUTE_LAPLACIAN_HEADER_FILE

#include <Eigen/Dense>

inline double cotangent(const Eigen::Vector3d &v0,
                        const Eigen::Vector3d &v1)
{

    double dot = v0.dot(v1);
    double crossNorm = (v0.cross(v1)).norm();

    return dot / crossNorm;
}

void compute_laplacian(const Eigen::MatrixXd &V,
                       const Eigen::MatrixXi &F,
                       const Eigen::MatrixXi &E,
                       const Eigen::MatrixXi &EF,
                       const Eigen::VectorXi &boundEMask,
                       Eigen::SparseMatrix<double> &d0,
                       Eigen::SparseMatrix<double> &W,
                       Eigen::VectorXd &vorAreas)
{

    using namespace Eigen;
    using namespace std;

    // TODO
    int nV = V.rows();
    int nE = E.rows();
    int nF = F.rows();

    d0.resize(nE, nV);
    W.resize(nE, nE);
    vorAreas = VectorXd::Zero(nV);

    for (int i = 0; i < nF; i++)
    {
        int i0 = F(i, 0);
        int i1 = F(i, 1);
        int i2 = F(i, 2);

        Vector3d v0 = V.row(i0);
        Vector3d v1 = V.row(i1);
        Vector3d v2 = V.row(i2);

        double area = 0.5 * ((v1 - v0).cross(v2 - v0)).norm();

        vorAreas(i0) += area / 3.0;
        vorAreas(i1) += area / 3.0;
        vorAreas(i2) += area / 3.0;
    }

    vector<Triplet<double>> d0Triplets;
    for (int e = 0; e < nE; e++)
    {
        int source = E(e, 0);
        int target = E(e, 1);
        d0Triplets.push_back(Triplet<double>(e, source, -1));
        d0Triplets.push_back(Triplet<double>(e, target, 1));
    }
    d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());

    vector<Triplet<double>> wTriplets;
    wTriplets.reserve(nE);
    for (int e = 0; e < nE; e++)
    {
        int v0 = E(e, 0);
        int v1 = E(e, 1);
        double cotSum = 0.0;

        int fIdx = EF(e, 0);
        if (fIdx >= 0)
        {
            int fv0 = F(fIdx, 0);
            int fv1 = F(fIdx, 1);
            int fv2 = F(fIdx, 2);
            int other;
            if (fv0 != v0 && fv0 != v1)
                other = fv0;
            else if (fv1 != v0 && fv1 != v1)
                other = fv1;
            else
                other = fv2;

            Vector3d vec0 = V.row(v0) - V.row(other);
            Vector3d vec1 = V.row(v1) - V.row(other);
            double c = cotangent(vec0, vec1);
            cotSum += c;
        }

        fIdx = EF(e, 2);
        if (fIdx >= 0)
        {
            int fv0 = F(fIdx, 0);
            int fv1 = F(fIdx, 1);
            int fv2 = F(fIdx, 2);
            int other;
            if (fv0 != v0 && fv0 != v1)
                other = fv0;
            else if (fv1 != v0 && fv1 != v1)
                other = fv1;
            else
                other = fv2;

            Vector3d vec0 = V.row(v0) - V.row(other);
            Vector3d vec1 = V.row(v1) - V.row(other);
            double c = cotangent(vec0, vec1);
            cotSum += c;
        }
        double weight = 0.5 * cotSum;

        wTriplets.push_back(Triplet<double>(e, e, weight));
    }
    W.resize(nE, nE);
    W.setFromTriplets(wTriplets.begin(), wTriplets.end());
}

#endif
