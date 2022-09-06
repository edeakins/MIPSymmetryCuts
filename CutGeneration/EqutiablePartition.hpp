#ifndef EQUITABLE_H
#define EQUITABLE_H

#include <string>
#include <vector>
#include <algorithm>
#include <stack>
#include <set>
#include <list>
#include <map>
#include <tuple>
#include <numeric>
#include <functional>
using namespace std;

class EquitablePartition{
public:
    EquitablePartition(const int nR, const int nC, const vector<double>& cC, const vector<double>& cL,
                            const vector<double>& cU, const vector<double>& rL, const vector<double>& rU,
                            const vector<double>& Av, const vector<int>& Ai, const vector<int>& As,
                            const vector<double>& ARv, const vector<int>& ARi, const vector<int>& ARs);
    void initRefinement();
    void handleNegatives();
    void refine();
    void splitColor(int s);
    void isolate(int i);
    void findTarget();
    void packVectors();
    bool isDiscrete();

    // Original LP info
    int nRows;
    int nCols; 
    int numTot;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper; 
    vector<double> rowLower;
    vector<double> rowUpper;
    vector<double> Avalue; 
    vector<int> Aindex; 
    vector<int> Astart;
    vector<double> AvaluePos;
    vector<double> ARvalue;
    vector<int> ARindex;
    vector<int> ARstart;
    vector<double> ARvaluePos;
    vector<int> AindexP;

    // EP info
    int vCol, cCol, refinements = 0;
    vector<vector<int> > C;
    vector<vector<int> > A;
    vector<int> color;
    vector<bool> SCheck;
    vector<double> cdeg;
	vector<double> mincdeg;
	vector<double> maxcdeg;
    vector<bool> isAdj;
    vector<int> colorsAdj;
    vector<int> colorsToSplit;
    vector<bool> isolates;
    stack<int> S;
    vector<int> Csize;
    vector<int> Asize;
    vector<int> initialParts;
    vector<int> numEdges;
};

#endif