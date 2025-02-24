#include "readOFF.h"
#include "tutte_parameterization.h"
#include "compute_areas_normals.h"
#include "compute_laplacian.h"
#include "create_edge_list.h"
#include "serialization.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <filesystem>

using namespace Eigen;
using namespace std;


double tolerance = 1e-3;

namespace fs = std::filesystem;


int main()
{
    double section2Points = 15.0;
    int pointGain=0;
    int pointSum=0;
    
    std::string folderPath(DATA_PATH "/param"); // Replace with your folder path
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".off") {
            cout<<"Working on file "<<entry.path().filename()<<endl;
            std::string dataName = entry.path().string();
            dataName.erase(dataName.size() - 4, 4);
            std::ifstream ifs(dataName+"-section3.data", std::ofstream::binary);
            
            MatrixXi F, E, EF;
            VectorXi boundEMask, boundVMask, boundVertices;
            MatrixXd V;
            SparseMatrix<double> d0, W;
            VectorXd vorAreas;
            MatrixXd UVBound, UV;
            
            readOFF(entry.path().string(), V, F);
            create_edge_list(F, E, EF, boundEMask, boundVMask, boundVertices, true);
            compute_laplacian(V, F, E, EF, boundEMask, d0, W, vorAreas);
            double r = sqrt(vorAreas.sum()/M_PI);
            
            auto start = std::chrono::high_resolution_clock::now();
            UVBound = compute_boundary_embedding(V, boundVertices, r);
            UV = compute_tutte_embedding(boundVertices, UVBound, d0, W);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Tutte Parameterization took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
            VectorXd durVector(1); durVector<<(double)(duration.count())/1000.0;
            
            MatrixXd UVBoundGT, UVGT;
            VectorXd durVectorGT;
            deserializeMatrix(UVBoundGT, ifs);
            deserializeMatrix(UVGT, ifs);
            deserializeVector(durVectorGT, ifs);
            
            pointSum+=2;
            if ((durVectorGT(0)*10.0 < (double)(duration.count())/1000.0)&&(durVectorGT(0)>1000.0)){
                cout<<"Running took too long! "<<endl;
            } else {
                int rowIndex, colIndex;
                if ((UVBoundGT-UVBound).cwiseAbs().maxCoeff(&rowIndex, &colIndex)<=tolerance){
                    cout<<"UVBound is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"UVBound("<<rowIndex<<","<<colIndex<<")="<<UVBound(rowIndex, colIndex)<<", Ground-truth UVBound("<<rowIndex<<","<<colIndex<<")="<<UVBoundGT(rowIndex, colIndex)<<endl;
                }
                if ((UVGT-UV).cwiseAbs().maxCoeff(&rowIndex, &colIndex)<=tolerance){
                    cout<<"Full UV is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"UV("<<rowIndex<<","<<colIndex<<")="<<UV(rowIndex, colIndex)<<", Ground-truth UV("<<rowIndex<<","<<colIndex<<")="<<UVGT(rowIndex, colIndex)<<endl;
                }
            }
        }
    }
    cout<<"Total point gained: "<<pointGain<<"/"<<pointSum<<endl;
    cout<<"Grade for Section 1: "<<round((double)pointGain*section2Points/(double)pointSum)<<endl;
    
}

