#ifndef COMPUTE_TUTTE_PARAMETERIZATION_HEADER_FILE
#define COMPUTE_TUTTE_PARAMETERIZATION_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "slice_columns_sparse.h"
#include "set_diff.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Eigen::MatrixXd compute_boundary_embedding(const Eigen::MatrixXd &V,
                                           const Eigen::VectorXi &boundVertices,
                                           const double r)
{

    // TODO

    using namespace Eigen;
    using namespace std;

    int nbV = boundVertices.size();
    MatrixXd UVBound(nbV, 2);

    if (nbV == 0)
    {
        return UVBound;
    }

    if (nbV == 1)
    {
        UVBound(0, 0) = r;
        UVBound(0, 1) = 0.0;
        return UVBound;
    }

    vector<double> edgeLength(nbV);
    double totalLength = 0.0;
    for (int i = 0; i < nbV; i++)
    {
        int iNext = (i + 1) % nbV;

        RowVector3d p0 = V.row(boundVertices[i]);
        RowVector3d p1 = V.row(boundVertices[iNext]);
        double len = (p1 - p0).norm();
        edgeLength[i] = len;
        totalLength += len;
    }

    vector<double> sectorAngles(nbV);
    for (int i = 0; i < nbV; i++)
    {
        sectorAngles[i] = 2 * M_PI * (edgeLength[i] / totalLength);
    }

    double sumAngle = 0.0;
    for (int i = 0; i < nbV; i++)
    {
        UVBound(i, 0) = r * cos(sumAngle);
        UVBound(i, 1) = r * sin(sumAngle);

        sumAngle += sectorAngles[i];
    }
    return UVBound;
}

Eigen::MatrixXd compute_tutte_embedding(const Eigen::VectorXi &boundVertices,
                                        const Eigen::MatrixXd &UVBound,
                                        const Eigen::SparseMatrix<double> &d0,
                                        const Eigen::SparseMatrix<double> &W)
{

    // TODO
    using namespace Eigen;
    using namespace std;

    const int bD0 = d0.cols();
    const int nbV = boundVertices.size();
    const int nI = bD0 - nbV;

    vector<bool> isBoundary(bD0, false);
    for (int i = 0; i < nbV; ++i)
    {
        isBoundary[boundVertices[i]] = true;
    }

    vector<int> interiorVerts;
    interiorVerts.reserve(nI);
    for (int v = 0; v < bD0; ++v)
    {
        if (!isBoundary[v])
        {
            interiorVerts.push_back(v);
        }
    }

    unordered_map<int, int> interiorMap;
    interiorMap.reserve(nI);
    for (int i = 0; i < nI; ++i)
    {
        interiorMap[interiorVerts[i]] = i;
    }

    unordered_map<int, int> boundaryMap;
    boundaryMap.reserve(nbV);
    for (int i = 0; i < nbV; ++i)
    {
        boundaryMap[boundVertices[i]] = i;
    }

    vector<Triplet<double>> tripI;
    vector<Triplet<double>> tripB;
    tripI.reserve(d0.nonZeros());
    tripB.reserve(d0.nonZeros());

    for (int k = 0; k < d0.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(d0, k); it; ++it)
        {
            const int e = it.row();
            const int v = it.col();
            const double val = it.value();

            if (isBoundary[v])
            {
                int newCol = boundaryMap[v];
                tripB.emplace_back(e, newCol, val);
            }
            else
            {
                int newCol = interiorMap[v];
                tripI.emplace_back(e, newCol, val);
            }
        }
    }

    SparseMatrix<double> d0I(d0.rows(), nI), d0B(d0.rows(), nbV);
    d0I.setFromTriplets(tripI.begin(), tripI.end());
    d0B.setFromTriplets(tripB.begin(), tripB.end());

    SparseMatrix<double> WD0I = W * d0I;
    SparseMatrix<double> A = d0I.transpose() * WD0I;
    SparseMatrix<double> WD0B = W * d0B;
    SparseMatrix<double> M = d0I.transpose() * WD0B;
    MatrixXd M_UVBound = M * UVBound;
    MatrixXd RHS = -M_UVBound;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(A);

    MatrixXd UVI(nI, 2);
    for (int c = 0; c < 2; ++c)
    {
        VectorXd rhsCol = RHS.col(c);
        VectorXd solCol = solver.solve(rhsCol);
        UVI.col(c) = solCol;
    }

    MatrixXd UV(bD0, 2);

    for (int i = 0; i < nbV; ++i)
    {
        int vB = boundVertices[i];
        UV(vB, 0) = UVBound(i, 0);
        UV(vB, 1) = UVBound(i, 1);
    }

    for (int i = 0; i < nI; ++i)
    {
        int vI = interiorVerts[i];
        UV(vI, 0) = UVI(i, 0);
        UV(vI, 1) = UVI(i, 1);
    }

    return UV;
}

#endif
