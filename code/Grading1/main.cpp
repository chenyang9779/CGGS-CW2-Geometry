#include "readOFF.h"
#include "compute_angle_defect.h"
#include "compute_mean_curvature_normal.h"
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
    double section1Points = 25.0;
    int pointGain=0;
    int pointSum=0;
    
    std::string folderPath(DATA_PATH); // Replace with your folder path
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".off") {
            cout<<"Working on file "<<entry.path().filename()<<endl;
            std::string dataName = entry.path().string();
            dataName.erase(dataName.size() - 4, 4);
            std::ifstream ifs(dataName+"-section1.data", std::ofstream::binary);
            
            MatrixXi F, E, EF;
            VectorXi boundEMask, boundVMask, boundVertices;
            MatrixXd V, Hn;
            SparseMatrix<double> d0, W;
            VectorXd vorAreas, H;
            
            readOFF(entry.path().string(), V, F);
            create_edge_list(F, E, EF, boundEMask, boundVMask, boundVertices);
            
            auto start = std::chrono::high_resolution_clock::now();
            VectorXd G = compute_angle_defect(V, F, boundVMask);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "compute_angle_defect() took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
            VectorXd durVector(1); durVector<<(double)(duration.count())/1000.0;
            VectorXd GGT;
            VectorXd durVectorGT;
            deserializeVector(GGT, ifs);
            deserializeVector(durVectorGT, ifs);
            
            pointSum++;
            if ((durVectorGT(0)*10.0 < (double)(duration.count())/1000.0)&&(durVectorGT(0)>1000.0)){
                cout<<"Running took too long! "<<endl;
            } else {
                int where;
                if ((GGT-G).cwiseAbs().maxCoeff(&where)<=tolerance){
                    cout<<"G is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"G("<<where<<")="<<G(where)<<", Ground-truth G("<<where<<")="<<GGT(where)<<endl;
                }
            }
            
            start = std::chrono::high_resolution_clock::now();
            compute_laplacian(V, F, E, EF, boundEMask, d0, W, vorAreas);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "compute_laplacian() took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
            VectorXd vorAreasGT;
            deserializeVector(vorAreasGT, ifs);
            deserializeVector(durVectorGT, ifs);
            pointSum++;
            if ((durVectorGT(0)*10.0 < (double)(duration.count())/1000.0)&&(durVectorGT(0)>1000.0)){
                cout<<"Running took too long! "<<endl;
            } else {
                int where;
                if ((vorAreasGT-vorAreas).cwiseAbs().maxCoeff(&where)<=tolerance){
                    cout<<"vorAreas is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"vorAreas("<<where<<")="<<vorAreas(where)<<", Ground-truth vorAreas("<<where<<")="<<vorAreasGT(where)<<endl;
                }
            }
            
            start = std::chrono::high_resolution_clock::now();
            compute_mean_curvature_normal(V, F, d0.transpose()*W*d0, vorAreas, Hn, H);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "compute_mean_curvature_normal() took " << (double)(duration.count())/1000.0 << " seconds to execute." << std::endl;
            MatrixXd HnGT;
            VectorXd HGT;
            deserializeMatrix(HnGT, ifs);
            deserializeVector(durVectorGT, ifs);
            pointSum++;
            if ((durVectorGT(0)*10.0 < (double)(duration.count())/1000.0)&&(durVectorGT(0)>1000.0)){
                cout<<"Running took too long! "<<endl;
            } else {
                int rowIndex, colIndex;
                if ((HnGT-Hn).cwiseAbs().maxCoeff(&rowIndex, &colIndex)<=tolerance){
                    cout<<"Hn is good!"<<endl;
                    pointGain++;
                } else {
                    cout<<"Hn("<<rowIndex<<","<<colIndex<<")="<<HnGT(rowIndex, colIndex)<<", Ground-truth Hn("<<rowIndex<<","<<colIndex<<")="<<Hn(rowIndex, colIndex)<<endl;
                }
            }
        }
    }
    cout<<"Total point gained: "<<pointGain<<"/"<<pointSum<<endl;
    cout<<"Grade for Section 1: "<<round((double)pointGain*section1Points/(double)pointSum)<<endl;
}

