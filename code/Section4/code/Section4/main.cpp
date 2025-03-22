#include <iostream>
#include <vector>
#include <unordered_set>
#include <Eigen/Dense>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include "readOFF.h"

using namespace Eigen;
using namespace std;

struct Mesh
{
    vector<Vector3d> vertices;
    vector<Vector3i> faces;
    vector<Vector3d> normals;
};

Mesh loadMesh(const string &filename)
{
    Mesh mesh;
    MatrixXd V;
    MatrixXi F;
    bool success = readOFF(filename, V, F);
    if (!success)
    {
        cerr << "Error reading OFF file: " << filename << "\n";
        exit(1);
    }
    mesh.vertices.resize(V.rows());
    for (int i = 0; i < V.rows(); ++i)
    {
        mesh.vertices[i] = Vector3d(V(i, 0), V(i, 1), V(i, 2));
    }
    mesh.faces.resize(F.rows());
    for (int i = 0; i < F.rows(); ++i)
    {
        mesh.faces[i] = Vector3i(F(i, 0), F(i, 1), F(i, 2));
    }
    return mesh;
}

void computeVertexNormals(Mesh &mesh)
{
    mesh.normals.resize(mesh.vertices.size(), Vector3d::Zero());
    for (auto &f : mesh.faces)
    {
        Vector3d e1 = mesh.vertices[f[1]] - mesh.vertices[f[0]];
        Vector3d e2 = mesh.vertices[f[2]] - mesh.vertices[f[0]];
        Vector3d faceNormal = e1.cross(e2);
        for (int i = 0; i < 3; i++)
        {
            mesh.normals[f[i]] += faceNormal;
        }
    }
    for (auto &n : mesh.normals)
    {
        n.normalize();
    }
}

void computePrincipalCurvatures(const Mesh &mesh, vector<Vector2d> &curvatures, vector<Vector3d> &dir1, vector<Vector3d> &dir2)
{
    size_t n = mesh.vertices.size();
    curvatures.resize(n);
    dir1.resize(n);
    dir2.resize(n);

    for (size_t i = 0; i < n; ++i)
    {   
        // Current vertex position & normal
        const Vector3d &v = mesh.vertices[i];
        const Vector3d &nm = mesh.normals[i];

        // Construct local frame R with x_axis, y_axis in tangent plane, nm as z-axis
        Vector3d x_axis = nm.cross(Vector3d::UnitY());

        x_axis = nm.cross(Vector3d::UnitX());

        x_axis.normalize();
        Vector3d y_axis = nm.cross(x_axis);
        Matrix3d R;
        R.col(0) = x_axis;
        R.col(1) = y_axis;
        R.col(2) = nm;

        // Collect 1-ring neighbor positions
        unordered_set<int> neighborIndices;
        for (const auto &f : mesh.faces)
        {
            if (f[0] == (int)i || f[1] == (int)i || f[2] == (int)i)
            {
                neighborIndices.insert(f[0]);
                neighborIndices.insert(f[1]);
                neighborIndices.insert(f[2]);
            }
        }
        neighborIndices.erase(i);

        // Build the system A * coeffs = b
        vector<Vector3d> neighbors;
        for (int idx : neighborIndices)
        {
            neighbors.push_back(mesh.vertices[idx]);
        }
        MatrixXd A(neighbors.size(), 6);
        VectorXd b(neighbors.size());
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
            // Local coords relative to v
            Vector3d p_local = R.transpose() * (neighbors[j] - v);
            double x = p_local.x();
            double y = p_local.y();
            double z = p_local.z();

            // Fill row in A and b
            A.row(j) << x * x, y * y, x * y, x, y, 1.0;
            b(j) = z;
        }

        // Solve via SVD
        JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
        VectorXd coeffs = svd.solve(b);

        // Hessian
        Matrix2d H;
        H << 2.0 * coeffs(0), coeffs(2), coeffs(2), 2.0 * coeffs(1);

        // Diagonalize H
        SelfAdjointEigenSolver<Matrix2d> es(H);
        curvatures[i] = es.eigenvalues();

        // Rotate eigenvectors back to global
        Eigen::Vector2d e1 = es.eigenvectors().col(0);
        Eigen::Vector2d e2 = es.eigenvectors().col(1);

        // Construct global directions
        dir1[i] = e1(0) * x_axis + e1(1) * y_axis;
        dir2[i] = e2(0) * x_axis + e2(1) * y_axis;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <mesh.off>" << endl;
        return 1;
    }

    Mesh mesh = loadMesh(argv[1]);
    computeVertexNormals(mesh);

    vector<Vector2d> curvatures;
    vector<Vector3d> dir1, dir2;
    computePrincipalCurvatures(mesh, curvatures, dir1, dir2);

    vector<double> c1(mesh.vertices.size()), c2(mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
    {
        c1[i] = curvatures[i](0);
        c2[i] = curvatures[i](1);
    }

    polyscope::init();
    auto *psMesh = polyscope::registerSurfaceMesh("Mesh", mesh.vertices, mesh.faces);

    psMesh->addVertexScalarQuantity("Curvature 1", c1);
    psMesh->addVertexScalarQuantity("Curvature 2", c2);

    psMesh->addVertexVectorQuantity("Principal Direction 1", dir1);
    psMesh->addVertexVectorQuantity("Principal Direction 2", dir2);

    polyscope::show();
    return 0;
}