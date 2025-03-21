#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <queue>
#include "readOFF.h"
#include "create_edge_list.h"
#include "compute_angle_defect.h"
#include "compute_laplacian.h"
#include "compute_mean_curvature_normal.h"
#include "compute_areas_normals.h"
#include "mean_curvature_flow.h"
#include <chrono>
#include <filesystem>
#include <cstdlib>
#include <string>

using namespace Eigen;
using namespace std;

// Global parameters with default values:
bool isFlowing;
bool isExplicit;
double timeStepRatio;
double timeStep;

polyscope::SurfaceMesh* psMesh;

MatrixXi F, E, EF;
VectorXi boundEMask, boundVMask, boundVertices;
MatrixXd origV, currV, Hn, faceNormals;
SparseMatrix<double> d0, W, M, MInv, L;
VectorXd vorAreas, H, faceAreas;

void callback_function() {
    ImGui::PushItemWidth(50);
    ImGui::TextUnformatted("Flow Parameters");
    ImGui::Separator();
    ImGui::Checkbox("Flow", &isFlowing);
    ImGui::PopItemWidth();
    if (!isFlowing)
        return;
    
    mean_curvature_flow(F, L, timeStep, M, MInv, boundVMask, isExplicit, currV);
    psMesh->updateVertexPositions(currV);
}

int main(int argc, char** argv)
{
    if (argc != 5 ) {
        cerr << "Usage: " << argv[0] << "isFlowing isExplicit timeStepRatio <mesh.off>" << endl;
        return 1;
    }

    string arg1(argv[1]);
    isFlowing = (arg1 == "true" || arg1 == "1");
    string arg2(argv[2]);
    timeStepRatio = atof(argv[3]);
    isExplicit = (arg2 == "true" || arg2 == "1");

    readOFF(argv[4], origV, F);
    currV = origV.rowwise() - origV.colwise().mean();
    create_edge_list(F, E, EF, boundEMask, boundVMask, boundVertices);
    
    polyscope::init();
    psMesh = polyscope::registerSurfaceMesh("Mesh", currV, F);
    compute_areas_normals(origV, F, faceAreas, faceNormals);
    timeStep = timeStepRatio * faceAreas.minCoeff();
    
    compute_laplacian(origV, F, E, EF, boundEMask, d0, W, vorAreas);
    L = d0.transpose() * W * d0;
    vector<Triplet<double>> MTris, MInvTris;
    for (int i = 0; i < vorAreas.size(); i++){
        MTris.push_back(Triplet<double>(i, i, vorAreas(i)));
        MInvTris.push_back(Triplet<double>(i, i, 1.0 / vorAreas(i)));
    }
    M.resize(vorAreas.size(), vorAreas.size());
    MInv.resize(vorAreas.size(), vorAreas.size());
    M.setFromTriplets(MTris.begin(), MTris.end());
    MInv.setFromTriplets(MInvTris.begin(), MInvTris.end());
    
    /*psMesh->addVertexScalarQuantity("Gaussian Curvature", G.array()/vorAreas.array())->setEnabled(true);
    psMesh->addVertexScalarQuantity("Gaussian Regions", G.unaryExpr([](double x) { return (x > 0) - (x < 0); }));
    psMesh->addVertexScalarQuantity("Mean Curvature", H);
    psMesh->addVertexVectorQuantity("Mean Curvature Normal", Hn);*/
    
    polyscope::state::userCallback = callback_function;
    polyscope::show();
}
