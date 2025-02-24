#include "readOFF.h"
#include "mean_curvature_flow.h"
#include "compute_areas_normals.h"
#include "compute_laplacian.h"
#include "create_edge_list.h"
#include "serialization.h"
#include <Eigen/Dense>
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
    
    std::string folderPath(DATA_PATH); // Replace with your folder path
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".off") {
            cout<<"Working on file "<<entry.path().filename()<<endl;
            std::string dataName = entry.path().string();
            dataName.erase(dataName.size() - 4, 4);
            std::ifstream ifs(dataName+"-section2.data", std::ofstream::binary);
            
            bool isExplicit;
            double timeStepRatio = 0.05;
            double timeStep;
            
            MatrixXi F, E, EF;
            VectorXi boundEMask, boundVMask, boundVertices;
            MatrixXd origV, currV, Hn, faceNormals;
            SparseMatrix<double> d0, W, M, MInv, L;
            VectorXd vorAreas, H, faceAreas;
            
            readOFF(entry.path().string(), origV, F);
            currV = origV;
            create_edge_list(F, E, EF, boundEMask, boundVMask, boundVertices);
            compute_areas_normals(origV, F, faceAreas, faceNormals);
            
            compute_laplacian(origV, F, E, EF, boundEMask, d0, W, vorAreas);
            L = d0.transpose()*W*d0;
            vector<Triplet<double>> MTris, MInvTris;
            for (int i=0;i<vorAreas.size();i++){
                MTris.push_back(Triplet<double>(i,i,vorAreas(i)));
                MInvTris.push_back(Triplet<double>(i,i,1.0/vorAreas(i)));
            }
            M.resize(vorAreas.size(), vorAreas.size());
            MInv.resize(vorAreas.size(), vorAreas.size());
            M.setFromTriplets(MTris.begin(), MTris.end());
            MInv.setFromTriplets(MInvTris.begin(), MInvTris.end());
            
            //Explicit flow
            isExplicit = true;
            timeStep = timeStepRatio*faceAreas.minCoeff();
            
            auto start = std::chrono::high_resolution_clock::now();
            for (int i=0;i<200;i++)
                mean_curvature_flow(F, L, timeStep, M, MInv, boundVMask, isExplicit, currV);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Explicit flow took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
            VectorXd durVector(1); durVector<<(double)(duration.count())/1000.0;
            
            MatrixXd currVGT;
            VectorXd durVectorGT;
            deserializeMatrix(currVGT, ifs);
            deserializeVector(durVectorGT, ifs);
            
            pointSum++;
            if ((durVectorGT(0)*10.0 < (double)(duration.count())/1000.0)&&(durVectorGT(0)>1000.0)){
                cout<<"Running took too long! "<<endl;
            } else {
                int rowIndex, colIndex;
                if ((currVGT-currV).cwiseAbs().maxCoeff(&rowIndex, &colIndex)<=tolerance){
                    cout<<"Explicit currV is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"currV("<<rowIndex<<","<<colIndex<<")="<<currV(rowIndex, colIndex)<<", Ground-truth currV("<<rowIndex<<","<<colIndex<<")="<<currVGT(rowIndex, colIndex)<<endl;
                }
            }
            
            //Implicit flow
            currV = origV;
            isExplicit = false;
            timeStep = 50*timeStepRatio*faceAreas.minCoeff();
            
            start = std::chrono::high_resolution_clock::now();
            for (int i=0;i<15;i++)
                mean_curvature_flow(F, L, timeStep, M, MInv, boundVMask, isExplicit, currV);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "Implicit flow took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
            
            deserializeMatrix(currVGT, ifs);
            deserializeVector(durVectorGT, ifs);
            
            pointSum++;
            if ((durVectorGT(0)*10.0 < (double)(duration.count())/1000.0)&&(durVectorGT(0)>1000.0)){
                cout<<"Running took too long! "<<endl;
            } else {
                int rowIndex, colIndex;
                if ((currVGT-currV).cwiseAbs().maxCoeff(&rowIndex, &colIndex)<=tolerance){
                    cout<<"Implicit currV is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"currV("<<rowIndex<<","<<colIndex<<")="<<currV(rowIndex, colIndex)<<", Ground-truth currV("<<rowIndex<<","<<colIndex<<")="<<currVGT(rowIndex, colIndex)<<endl;
                }
            }
        }
    }
    cout<<"Total point gained: "<<pointGain<<"/"<<pointSum<<endl;
    cout<<"Grade for Section 1: "<<round((double)pointGain*section2Points/(double)pointSum)<<endl;
    
}

