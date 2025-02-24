#ifndef CREATE_EDGE_LIST_HEADER_FILE
#define CREATE_EDGE_LIST_HEADER_FILE

#include <Eigen/Dense>
#include "unique.h"
#include "sort_rows.h"


#include <iostream>
#include <Eigen/Dense>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>

Eigen::VectorXi extract_boundary_loop(const Eigen::MatrixXi& E, const Eigen::VectorXi& boundEMask) {
    // Step 1: Extract boundary edges into an ordered adjacency map (std::map ensures deterministic order)
    std::map<int, std::vector<int>> adjacency;
    for (int i = 0; i < E.rows(); ++i) {
        if (boundEMask[i] == 1) { // This edge is a boundary edge
            int u = E(i, 0), v = E(i, 1);
            adjacency[u].push_back(v);
            adjacency[v].push_back(u);
        }
    }

    // Step 2: Sort neighbors to ensure deterministic traversal order
    for (auto& [_, neighbors] : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
    }

    // Step 3: Find a valid starting vertex
    int start = -1;
    for (const auto& [vertex, neighbors] : adjacency) {
        if (neighbors.size() == 1) { // Prefer an endpoint if available
            start = vertex;
            break;
        }
    }
    if (start == -1) { // If no endpoints exist (closed loop), pick the smallest boundary vertex
        start = adjacency.begin()->first;
    }

    // Step 4: Traverse the boundary loop in order
    std::vector<int> loop;
    std::set<int> visited;
    int current = start;
    int prev = -1;  // Previous vertex to avoid stepping backward

    while (true) {
        loop.push_back(current);
        visited.insert(current);

        int next = -1;
        for (int neighbor : adjacency[current]) {
            if (neighbor != prev) { // Move forward, avoiding the previous step
                next = neighbor;
                break;
            }
        }

        if (next == -1 || visited.count(next)) break; // Stop if we can't continue

        prev = current;
        current = next;

        if (current == start && loop.size() > 1) break; // Ensure full loop
    }

    // Step 5: Convert result to Eigen VectorXi
    Eigen::VectorXi orderedLoop(loop.size());
    for (size_t i = 0; i < loop.size(); ++i) {
        orderedLoop(i) = loop[i];
    }

    return orderedLoop;
}


void create_edge_list(const Eigen::MatrixXi& F,
                      Eigen::MatrixXi& E,
                      Eigen::MatrixXi& EF,
                      Eigen::VectorXi& boundEMask,
                      Eigen::VectorXi& boundVMask,
                      Eigen::VectorXi& boundVertices,
                      const bool computeBoundaryLoop = false){
    
    MatrixXi H(3*F.rows(),2);  //halfedges
    VectorXi twinH = VectorXi::Constant(3*F.rows(), -1);
    for (int i=0;i<F.rows();i++){
        H.row(3*i)<<F(i,0), F(i,1);
        H.row(3*i+1)<<F(i,1), F(i,2);
        H.row(3*i+2)<<F(i,2), F(i,0);
    }
    /*Eigen::MatrixXi HSorted=H;
    sort_rows(HSorted);
    std::vector<int> uniqueIndices,counts,inverse;
    unique(HSorted, uniqueIndices, counts, inverse);*/
    
    //finding edge twins
    struct ComparePairs {
        bool operator()(const std::pair<std::pair<int, int>, int>& a, const std::pair<std::pair<int, int>, int>& b) const {
            if (a.first.first == b.first.first) {
                return a.first.second < b.first.second;
            } else {
                return a.first.first < b.first.first;
            }
        }
    };
    
    //finding twins
    typedef std::pair<std::pair<int, int>, int> pairPlusOne;
    std::set<pairPlusOne, ComparePairs> edgeSet;
    std::vector<int> EHList;
    for (int i=0;i<H.rows();i++){
        std::pair<int,int> oppEdge(H(i,1), H(i,0));
        pairPlusOne oppEdgePlus(oppEdge, -1);
        std::set<pairPlusOne>::iterator si = edgeSet.find(oppEdgePlus);
        if (si == edgeSet.end()) {
            edgeSet.insert(pairPlusOne(std::pair<int, int>(H(i,0), H(i,1)), i));
            EHList.push_back(i);
        } else {  //found matching twin
            twinH[si->second] = i;
            twinH[i] = si->second;
        }
    }
    
    E.resize(EHList.size(),2);
    EF = MatrixXi::Constant(EHList.size(),4,-1);
    boundEMask.resize(E.rows());
    boundVMask=VectorXi::Zero(F.maxCoeff()+1);  
    for (int i=0;i<EHList.size();i++){
        E.row(i)<<H(EHList[i],0), H(EHList[i],1);
        EF(i,0)=EHList[i] / 3;
        EF(i,1)=EHList[i] % 3;
        boundEMask(i) = (twinH(EHList[i]) == -1);
        if (boundEMask(i)==1){
            boundVMask(E(i,0))=boundEMask(i);
            boundVMask(E(i,1))=boundEMask(i);
        }
        if (!boundEMask(i)){
            EF(i,2)=twinH(EHList[i]) / 3;
            EF(i,3)=twinH(EHList[i]) % 3;
        }
    }
    if (computeBoundaryLoop)
        boundVertices = extract_boundary_loop(E, boundEMask);
    
    //std::cout<<"boundVertices: "<<boundVertices.head(5)<<std::endl;
    //std::cout<<"sum boundEMask :"<<boundEMask.sum()<<std::endl;
    //std::cout<<"sum boundVMask :"<<boundVMask.sum()<<std::endl;
}



#endif
