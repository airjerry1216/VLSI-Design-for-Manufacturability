#ifndef LIB_H_
#define LIB_H_
#include <string>
#include <vector>
using namespace std;
class CuttingSet {
public:
    vector<int> polygons;
    vector<vector<int> > nodes;
    //vector<vector<int> > nodesAll;
    //vector<bool> nodesAllInvalid;
    int nodeStartIndex;
    CuttingSet();
};
class Rect {
public:
    vector<double> coord;
    Rect();
};
class Polygon {
public:
    string name;
    string mask;
    //string layer;
    vector<Rect> rectArray;
    Rect boundingRect;
    bool assign;
    Polygon();
    void checkBoundingRect();
};
class Row {
public:
    string name;
    vector<double> coord;
    vector<Polygon> polygonArray;
    //vector<Polygon> polygonArray; // sorted by x
    vector<int> mapToIndex;
    vector<vector<bool> > conflictGraph;
    vector<CuttingSet> cuttingSetArray;
    vector<vector<int> > solutionGraph;
    vector<int> prev;
    int polygonNum;
    int prevPolygonCnt;
    Row();
    void dummyExtension();
    void createCuttingSet();
    void createSolutionGraph();
};

#endif