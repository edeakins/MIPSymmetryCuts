#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "EquitablePartition.hpp"
#include <limits>

class AggregateLp{
public:
    AggregateLp(EquitablePartition& ep);
    void updateMasterLpAndEp(EquitablePartition& ep, int _nC, int _nR,
                            int _nnz, vector<int>& As, vector<int>& Ai,
                            vector<double>& Av, vector<double>& rL, vector<double>& rU);
    void clear();
    void findDimensions();
    void pairCutWithIndex();
    void pairOrderConstraintWithIndex();
    void aggregateColBounds();
    void aggregateRowBounds();
    void aggregateAMatrix();
    void addCutsToAggregate();
    void addOrderConstraints();
    void aggregateCostVector();
    void aggregate();
    int getNumCol();
    int getNumRowAfterCuts();
    int getNumRowAfterOrderConstraints();
    int getNumRow();
    vector<double>& getColUpper();
    vector<double>& getColLower();
    vector<double>& getColCost();
    vector<double>& getRowUpper();
    vector<double>& getRowLower();
    vector<double>& getAvalue();
    vector<int>& getAindex();
    vector<int>& getAstart();

    // Infinity
    const double inf = std::numeric_limits<double>::max();

    // Original LP info
    int nRows = 0;
    int _nRows = 0;
    int nCols = 0; 
    int _nCols = 0;
    int numTot = 0;
    int _numTot = 0;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper; 
    vector<double> rowLower;
    vector<double> _rowLower;
    vector<double> rowUpper;
    vector<double> _rowUpper;
    vector<double> Avalue;
    vector<double> _Avalue; 
    vector<int> Aindex; 
    vector<int> _Aindex;
    vector<int> Astart;
    vector<int> _Astart;
    vector<int> AindexP;

    // Ordering constraints to break symmetries
    vector<int> parentPartition;

    // Reduced LP 
    int nRows_ = 0;
    int tempNRows_ = 0;
    int nCols_ = 0; 
    int tempNCols_ = 0;
    int numTot_ = 0;
    int nCuts = 0;
    int nOrdConstraints = 0;
    vector<double> colCost_;
    vector<double> colLower_;
    vector<double> colUpper_; 
    vector<double> rowLower_;
    vector<double> rowUpper_;
    vector<double> Avalue_; 
    vector<int> Aindex_; 
    vector<int> Astart_;
    set<int> cut;
    map<int, int> cutIdx;
    map<int, int> orderConstraintIdx;

    // EP info
    vector<vector<int> > C;
    vector<int> Csize;
    vector<int> color;

};

#endif 