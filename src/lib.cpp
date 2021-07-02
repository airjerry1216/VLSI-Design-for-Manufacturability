#include <iostream>
#include <climits>
#include <cmath>
#include "lib.h"
using namespace std;

CuttingSet::CuttingSet()
{

}
Rect::Rect()
{
    this->coord.resize(4);
}
Polygon::Polygon()
{
    this->boundingRect.coord[0] = INT_MAX;
    this->boundingRect.coord[2] = 0;
    this->boundingRect.coord[1] = INT_MAX;
    this->boundingRect.coord[3] = 0;
    this->assign = 0;
}
void Polygon::checkBoundingRect()
{
    int index = this->rectArray.size() - 1;
    if (boundingRect.coord[0] > this->rectArray[index].coord[0])
        boundingRect.coord[0] = this->rectArray[index].coord[0];
    if (boundingRect.coord[2] < this->rectArray[index].coord[2])
        boundingRect.coord[2] = this->rectArray[index].coord[2];
    if (boundingRect.coord[1] > this->rectArray[index].coord[1])
        boundingRect.coord[1] = this->rectArray[index].coord[1];
    if (boundingRect.coord[3] < this->rectArray[index].coord[3])
        boundingRect.coord[3] = this->rectArray[index].coord[3];
}
Row::Row()
{
    this->coord.resize(2);
    this->polygonNum = 0;
}
void Row::dummyExtension()
{   
    //#pragma omp parallel for
    for (int i = 0; i < this->polygonNum; ++i) {
        int j = 0;
        for (j = this->polygonNum - 1; j > i; --j) {
            if (this->conflictGraph[i][j])
                break;
        }
        if (j > i) {
            double newBoundary = this->polygonArray[j].boundingRect.coord[0] - 0.000001;
            if (this->polygonArray[i].boundingRect.coord[2] < newBoundary) {
                this->polygonArray[i].boundingRect.coord[2] = newBoundary;
            }
        }
    }

}
void Row::createCuttingSet()
{
    int index = 0;
    for (int i = 0; i < this->polygonNum; ++i) {
        double x = this->polygonArray[i].boundingRect.coord[0]; // cutting line
        CuttingSet cuttingSet;
        for (int j = 0; j < this->polygonNum; ++j) {
            if (this->polygonArray[j].boundingRect.coord[0] > this->polygonArray[i].boundingRect.coord[0])
                break;
            if (this->polygonArray[j].boundingRect.coord[0] <= x && x <= this->polygonArray[j].boundingRect.coord[2]) {
                cuttingSet.polygons.push_back(j);
            }
        }
        //cuttingSet.nodesAll.resize(pow(4, cuttingSet.polygons.size()));
        //cuttingSet.nodesAllInvalid.resize(cuttingSet.nodesAll.size(), 0);
        vector<vector<int> > nodesAll;
        vector<bool> nodesAllInvalid;
        nodesAll.resize(pow(4, cuttingSet.polygons.size()));
        nodesAllInvalid.resize(nodesAll.size(), 0);
        for (int j = 0; j < nodesAll.size(); ++j) {
            nodesAll[j].resize(cuttingSet.polygons.size());
        }
        /******  create permutation  ******/
        int solution = 1;
        int step = 1, cnt = 0;
        for (int j = cuttingSet.polygons.size() - 1; j >= 0; --j) {
            int cnt = 0;
            for (int k = 0; k < nodesAll.size(); ++k) {
                nodesAll[k][j] = solution;
                ++cnt;
                if (cnt >= step) {
                    cnt = 0;
                    if (solution == 4)
                        solution = 1;
                    else
                        solution += 1;
                }
            }
            step *= 4;
        }
        /******  remove conflict node  ******/
        for (int j = 0; j < nodesAll.size(); ++j) {
            for (int k = 0; k < cuttingSet.polygons.size() - 1; ++k) {
                if (nodesAllInvalid[j])
                    break;
                for (int l = k + 1; l < cuttingSet.polygons.size(); ++l) {
                    // conflict polygons but same color
                    if (this->conflictGraph[cuttingSet.polygons[k]][cuttingSet.polygons[l]] && nodesAll[j][k] == nodesAll[j][l] && nodesAll[j][k] != 4) {
                        nodesAllInvalid[j] = 1;
                        break;
                    }
                }
            }
            if (!nodesAllInvalid[j])
                cuttingSet.nodes.push_back(nodesAll[j]);
        }
        cuttingSet.nodeStartIndex = index + 1;
        index += cuttingSet.nodes.size();
        this->cuttingSetArray.push_back(cuttingSet);
    }
}
void Row::createSolutionGraph() {
    int SGSize = this->cuttingSetArray[this->cuttingSetArray.size()-1].nodeStartIndex + this->cuttingSetArray[this->cuttingSetArray.size()-1].nodes.size() + 2 - 1;
    this->solutionGraph.resize(SGSize);
    for (int i = 0; i < SGSize; ++i)
        this->solutionGraph[i].resize(SGSize, INT_MAX);
    //cout << "SG size: " << this->solutionGraph.size() << endl;
    /******  create edges  ******/
    //#pragma omp parallel for
    for (int k = 0; k < this->cuttingSetArray[0].nodes.size(); ++k) {
        int VSB = 0;
        for (int n = 0; n < this->cuttingSetArray[0].polygons.size(); ++n) {
            if (this->cuttingSetArray[0].nodes[k][n] == 4) {
                //VSB += this->polygonArray[this->cuttingSetArray[0].polygons[n]].rectArray.size();
                ++VSB;
            }
        }
        this->solutionGraph[0][this->cuttingSetArray[0].nodeStartIndex+k] = VSB;
    }
    //#pragma omp parallel for
    for (int j = 0; j < this->cuttingSetArray[this->cuttingSetArray.size()-1].nodes.size(); ++j) {
        this->solutionGraph[this->cuttingSetArray[this->cuttingSetArray.size()-1].nodeStartIndex+j][SGSize-1] = 0;
    }
    for (int i = 1; i < this->cuttingSetArray.size(); ++i) {
        // j ,k : nodes pair
        //cout << "------ i: " << i << " ------" << endl;
        for (int j = 0; j < this->cuttingSetArray[i-1].nodes.size(); ++j) {
            for (int k = 0; k < this->cuttingSetArray[i].nodes.size(); ++k) {
                //cout << "---- j: " << j << " k: " << k << " ----" << endl;
                bool validEdge = 1;
                // m, n : polygons pair
                for (int m = 0; m < this->cuttingSetArray[i-1].polygons.size(); ++m) {
                    if (!validEdge)
                        break;
                    for (int n = 0; n < this->cuttingSetArray[i].polygons.size(); ++n) {
                        /*if (this->conflictGraph[this->cuttingSetArray[i-1].polygons[m]][this->cuttingSetArray[i].polygons[n]]) {
                            if (this->cuttingSetArray[i-1].nodes[j][m] == this->cuttingSetArray[i].nodes[k][n] && this->cuttingSetArray[i-1].nodes[j][m] != 4) {
                                validEdge = 0;
                                break;
                            }
                        }
                        if (this->cuttingSetArray[i-1].polygons[m] == this->cuttingSetArray[i].polygons[n]) {
                            if (this->cuttingSetArray[i-1].nodes[j][m] != this->cuttingSetArray[i].nodes[k][n] && this->cuttingSetArray[i-1].nodes[j][m] != 4 && this->cuttingSetArray[i].nodes[k][n] != 4) {
                                validEdge = 0;
                                break;
                            }
                            if (this->cuttingSetArray[i-1].nodes[j][m] == 4 && this->cuttingSetArray[i].nodes[k][n] != 4) {
                                validEdge = 0;
                                break;
                            }
                        }*/

                        // same polygon
                        if (this->cuttingSetArray[i-1].polygons[m] == this->cuttingSetArray[i].polygons[n]) {
                            // differnt color
                            if (this->cuttingSetArray[i-1].nodes[j][m] != this->cuttingSetArray[i].nodes[k][n]) {
                                validEdge = 0;
                                break;
                            }        
                        }
                        // different polygon
                        else {
                            // same color and both not ebeam
                            if (this->cuttingSetArray[i-1].nodes[j][m] == this->cuttingSetArray[i].nodes[k][n] && this->cuttingSetArray[i-1].nodes[j][m] != 4) {
                                // conflict polygons
                                if (this->conflictGraph[this->cuttingSetArray[i-1].polygons[m]][this->cuttingSetArray[i].polygons[n]]) {
                                    validEdge = 0;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (validEdge) {
                    //cout << "*i: " << i << " j: " << j << " k: " << k << endl;
                    //cout << this->cuttingSetArray[i-1].nodeStartIndex+j << " " << this->cuttingSetArray[i].nodeStartIndex+k << endl;
                    int VSB = 0;
                    /*for (int n = 0; n < this->cuttingSetArray[i].polygons.size(); ++n) {
                        if (this->cuttingSetArray[i].nodes[k][n] == 4)
                            VSB += this->polygonArray[this->cuttingSetArray[i].polygons[n]].rectArray.size();
                    }*/
                    for (int n = 0; n < this->cuttingSetArray[i].polygons.size(); ++n) {
                        if (this->cuttingSetArray[i].nodes[k][n] == 4 && this->cuttingSetArray[i].polygons[n] == i) {
                            //VSB += this->polygonArray[this->cuttingSetArray[i].polygons[n]].rectArray.size();
                            ++VSB;
                        }            
                    }
                    this->solutionGraph[this->cuttingSetArray[i-1].nodeStartIndex+j][this->cuttingSetArray[i].nodeStartIndex+k] = VSB;
                }
            }
        }
    }
}