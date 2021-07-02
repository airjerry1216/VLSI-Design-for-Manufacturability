#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <thread>
#include <ctime>
#include <chrono>
//#include <omp.h>
#include "lib.h"
using namespace std;

bool compare(const Polygon &a, const Polygon &b)
{
    return a.boundingRect.coord[0] < b.boundingRect.coord[0];
}
double EuclideanDis(double x1, double y1, double x2, double y2)
{
    double xDif = x2 - x1;
    double yDif = y2 - y1;
    return sqrt(xDif * xDif + yDif * yDif);
}
double rectDis(Rect &a, Rect &b)
{
    bool left = b.coord[2] < a.coord[0];
    bool right = a.coord[2] < b.coord[0];
    bool bottom = b.coord[3] < a.coord[1];
    bool top = a.coord[3] < b.coord[1];
    double recDis = 0;

    if (left && top)
    {
        recDis = EuclideanDis(a.coord[0], a.coord[3], b.coord[2], b.coord[1]);
    }
    else if (left && bottom)
    {
        recDis = EuclideanDis(a.coord[0], a.coord[1], b.coord[2], b.coord[3]);
    }
    else if (right && bottom)
    {
        recDis = EuclideanDis(a.coord[2], a.coord[1], b.coord[0], b.coord[3]);
    }
    else if (right && top)
    {
        recDis = EuclideanDis(a.coord[2], a.coord[3], b.coord[0], b.coord[1]);
    }
    else if (left)
    {
        recDis = a.coord[0] - b.coord[2];
    }
    else if (right)
    {
        recDis = b.coord[0] - a.coord[2];
    }
    else if (bottom)
    {
        recDis = a.coord[1] - b.coord[3];
    }
    else if (top)
    {
        recDis = b.coord[1] - a.coord[3];
    }
    return recDis;
}
double polygonDis(Polygon &A, Polygon &B)
{
    double minDis = INT_MAX, recDis = 0;
    for (auto a : A.rectArray)
    {
        for (auto b : B.rectArray)
        {
            recDis = rectDis(a, b);
            if (minDis > recDis)
            {
                minDis = recDis;
            }
        }
    }
    return minDis;
}
void createConflictGraph(Row &row, const double &dis)
{
    int graphSize = row.polygonNum;
    row.conflictGraph.resize(graphSize);
    row.mapToIndex.resize(graphSize);
    for (int i = 0; i < graphSize; ++i)
    {
        row.conflictGraph[i].resize(graphSize, 0);
        row.mapToIndex[stoi(row.polygonArray[i].name.substr(1)) - row.prevPolygonCnt] = i;
    }
    for (int i = 0; i < graphSize - 1; ++i)
    {
        for (int j = i + 1; j < graphSize; ++j)
        {
            if (rectDis(row.polygonArray[i].boundingRect, row.polygonArray[j].boundingRect) > dis)
                continue;
            double polyDis = polygonDis(row.polygonArray[i], row.polygonArray[j]);
            //cout << row.polygonArray[i].name << " " << row.polygonArray[j].name << " " << polyDis << endl;
            if (polyDis <= dis + 0.000001)
            {
                row.conflictGraph[i][j] = 1;
                row.conflictGraph[j][i] = 1;
            }
        }
    }
}
void Dijkstra(Row &row)
{
    vector<int> dist(row.solutionGraph.size(), INT_MAX);
    row.prev.resize(row.solutionGraph.size(), 0);
    typedef pair<int, int> pi;
    priority_queue<pi, vector<pi>, greater<pi> > pq;
    dist[0] = 0;
    pq.push(make_pair(0, 0));

    while (!pq.empty())
    {
        pi top = pq.top();
        pq.pop();
        int u = top.second;
        for (int v = 0; v < row.solutionGraph.size(); ++v)
        {
            if (row.solutionGraph[u][v] < dist[v] - dist[u])
            {
                dist[v] = dist[u] + row.solutionGraph[u][v];
                pq.push(make_pair(dist[v], v));
                row.prev[v] = u;
            }
        }
    }
    //cout << "final dis: " << dist[row.solutionGraph.size() - 1] << endl;
}
void constructPath(Row &row)
{
    int t = row.solutionGraph.size() - 1, column = row.cuttingSetArray.size() - 1;
    for (int i = 0; i < row.cuttingSetArray.size(); ++i)
    {
        //cout << "i: " << i << " " << t << " ";
        //cout << row.prev[t] << " ";
        t = row.prev[t];
        //cout << t - row.cuttingSetArray[column].nodeStartIndex << endl;
        int nodeIndex = t - row.cuttingSetArray[column].nodeStartIndex;
        for (int j = 0; j < row.cuttingSetArray[column].polygons.size(); ++j)
        {
            if (!row.polygonArray[row.cuttingSetArray[column].polygons[j]].assign)
            {
                if (row.cuttingSetArray[column].nodes[nodeIndex][j] == 1)
                    row.polygonArray[row.cuttingSetArray[column].polygons[j]].mask = "A";
                else if (row.cuttingSetArray[column].nodes[nodeIndex][j] == 2)
                    row.polygonArray[row.cuttingSetArray[column].polygons[j]].mask = "B";
                else if (row.cuttingSetArray[column].nodes[nodeIndex][j] == 3)
                    row.polygonArray[row.cuttingSetArray[column].polygons[j]].mask = "C";
                row.polygonArray[row.cuttingSetArray[column].polygons[j]].assign = 1;
            }
        }
        --column;
    }
}
void TPL(Row &row, double dis)
{
    //cout << row.name << endl;
    sort(row.polygonArray.begin(), row.polygonArray.end(), compare);
    createConflictGraph(row, dis);
    row.dummyExtension();
    row.createCuttingSet();
    row.createSolutionGraph();
    Dijkstra(row);
    constructPath(row);
}
int main(int argc, char *argv[])
{

    //time_t start;
    chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    /*start = clock();
    const auto processor_count = thread::hardware_concurrency();
    vector<thread> threads(processor_count);
    cout << processor_count << endl;
    clock_t initTime = clock();*/
    double dis = atof(argv[1]);
    ifstream input_file1(argv[2]);
    ifstream input_file2(argv[2]);
    ofstream output_file(argv[3]);
    stringstream ss;
    string line1, tmp, line2;
    int polygonCnt;
    //vector<Row> rows;

    /********** For Row **********/
    while (getline(input_file1, line1))
    {
        Row row;
        ss << line1;
        ss >> tmp >> row.name >> row.coord[0] >> row.coord[1];
        ss.str("");
        ss.clear();
        /********** For Polygon **********/
        while (getline(input_file1, line1))
        {
            ss << line1;
            ss >> tmp;
            if (tmp != "POLYGON")
            {
                ss.str("");
                ss.clear();
                break;
            }
            Polygon polygon;
            ss >> polygon.name;
            ss.str("");
            ss.clear();
            getline(input_file1, line1);
            ss << line1;
            ss >> tmp >> polygon.mask;
            ss.str("");
            ss.clear();
            getline(input_file1, line1);
            ss << line1;
            ss >> tmp >> tmp;
            ss.str("");
            ss.clear();
            /********** For Rect **********/
            while (getline(input_file1, line1))
            {
                ss << line1;
                ss >> tmp;
                if (tmp != "RECT")
                {
                    ss.str("");
                    ss.clear();
                    break;
                }
                Rect rect;
                ss >> rect.coord[0] >> rect.coord[1] >> rect.coord[2] >> rect.coord[3];
                polygon.rectArray.push_back(rect);
                polygon.checkBoundingRect();
                ss.str("");
                ss.clear();
            }
            getline(input_file1, line1);
            row.polygonArray.push_back(polygon);
        }
        row.polygonNum = row.polygonArray.size();
        row.prevPolygonCnt = polygonCnt;
        //rows.push_back(row);
        TPL(row, dis);
        /*for (int i = 0; i < row.polygonNum; ++i)
            cout << i << " " << row.polygonArray[i].name << " " << row.polygonArray[i].mask << endl;*/
        //int r = 0;
        while (getline(input_file2, line2)) {
            //cout << line2 << endl;
            output_file << line2 << "\n";
            while (getline(input_file2, line2)) {
                output_file << line2 << "\n";
                ss << line2;
                ss >> tmp;
                if (tmp != "POLYGON") {
                    ss.str("");
                    ss.clear();
                    break;     
                }
                string polygonName;
                ss >> polygonName;
                getline(input_file2, line2);
                //cout << polygonName << " " << rows[r].polygonArray[rows[r].mapToIndex[stoi(polygonName.substr(1)) - rows[r].prevPolygonCnt]].mask << endl;
                //output_file << "\t\tMASK " << rows[r].polygonArray[rows[r].mapToIndex[stoi(polygonName.substr(1)) - rows[r].prevPolygonCnt]].mask << "\n";
                output_file << "\t\tMASK " << row.polygonArray[row.mapToIndex[stoi(polygonName.substr(1)) - row.prevPolygonCnt]].mask << "\n";
                ss.str("");
                ss.clear();
                getline(input_file2, line2);
                output_file << line2 << "\n";
                while (getline(input_file2, line2)) {
                    output_file << line2 << "\n";
                    ss << line2;
                    ss >> tmp;
                    if (tmp != "RECT") {
                        ss.str("");
                        ss.clear();
                        break;
                    }
                    ss.str("");
                    ss.clear();
                }
                getline(input_file2, line2);
                output_file << line2 << "\n";
            }
            //++r;
            break;
        }
        //break; // only test ROW R0
        polygonCnt += row.polygonNum;
    }
    /*for (size_t t = 0; t < processor_count; ++t) {
        cout << "t: " << t << endl;
        threads[t] = thread(bind([&](int a, int b, int t) {
            for (int i = a; i < b; ++i) {
                cout << "i: " << i << endl;
                TPL(rows[i], dis);
            }
        }, t * rows.size() / processor_count, (t + 1) == processor_count ? rows.size() : (t + 1) * rows.size() / processor_count, t));
    }*/
    /*for(int i = 0; i < rows.size(); i++){
        threads.push_back(thread(TPL, ref(rows[i]), dis));
    }*/
    /*for(int i = 0; i < threads.size(); i++){
        if (threads[i].joinable()) {
        threads[i].join();
      }
    }*/
    //TPL(rows[0], dis);
    /*int n = rows.size(), m = 8;
    for (int i = 0; i < n; i += m)
    {
        #pragma omp parallel for
        for (int j = 0; j < m; ++j)
        {
            if (i + j < n)
                TPL(rows[i + j], dis);
        }
    }*/
    /*for (int i = 0; i < processor_count; ++i) {
        thread t1(TPL, ref(rows[i]), dis);
        t1.join();
    }*/
    //thread t1(TPL, ref(rows[0]), dis);
    input_file1.close();
    input_file2.close();
    output_file.close();
    /*clock_t endTime = clock();
    cout << "time: " << ((double)(endTime - initTime)) / CLOCKS_PER_SEC << endl;*/
    cout << "time(chrono): " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - t1).count() << "\n";
    return 0;
}