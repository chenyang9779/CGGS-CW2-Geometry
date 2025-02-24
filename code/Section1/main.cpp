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
#include <chrono>
#include <filesystem>


using namespace Eigen;
using namespace std;

MatrixXi F, E, EF;
VectorXi boundEMask, boundVMask, boundVertices;
MatrixXd V, Hn;
SparseMatrix<double> d0, W;
VectorXd vorAreas, H;


int main()
{
    readOFF(DATA_PATH "/cheburashka.off",V, F);
    create_edge_list(F, E, EF, boundEMask, boundVMask, boundVertices);

    polyscope::init();
    polyscope::SurfaceMesh* psMesh = polyscope::registerSurfaceMesh("Mesh", V, F);
    
    auto start = std::chrono::high_resolution_clock::now();
    VectorXd G = compute_angle_defect(V, F, boundVMask);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Function took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    compute_laplacian(V, F, E, EF, boundEMask, d0, W, vorAreas);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      std::cout << "Function took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
    
    
    start = std::chrono::high_resolution_clock::now();
    compute_mean_curvature_normal(V, F, d0.transpose()*W*d0, vorAreas, Hn, H);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Function took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
   
    psMesh->addVertexScalarQuantity("Gaussian Curvature", G.array()/vorAreas.array())->setEnabled(true);
    psMesh->addVertexScalarQuantity("Gaussian Regions", G.unaryExpr([](double x) { return (x > 0) - (x < 0); }));
    psMesh->addVertexScalarQuantity("Mean Curvature", H);
    psMesh->addVertexVectorQuantity("Mean Curvature Normal", Hn);
    
    polyscope::show();
    
}

