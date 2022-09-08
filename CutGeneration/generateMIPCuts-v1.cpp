#include "generateMIPCuts.hpp"
using namespace std;

int iter = 0;
int nCuts = 0;
int nCutsAdded = 0;

void update(int& nCols, int& nRows, int nnz, vector<int>& Astart,
            vector<int>& Aindex, vector<double>& Avalue, vector<int>& ARstart,
            vector<int>& ARindex, vector<double>& ARvalue, vector<double>& colLower,
            vector<double>& colUpper, vector<double>& rowLower, vector<double>& rowUpper,
            EquitablePartition& ep, AggregateLp& alp, OsiClpSolverInterface& aggSi, OsiClpSolverInterface& si){
  // Grab info we will need for update
  int nAggRowsWCuts = alp.getNumRowAfterOrderConstraints() + nCutsAdded;
  nCutsAdded = 0;
  int nAggColsWCuts = aggSi.getNumCols();
  int nnzAggWCuts = aggSi.getMatrixByRow()->getVectorStarts()[nAggRowsWCuts];
  int nAggRows = alp.getNumRowAfterOrderConstraints();
  int nAggCols = alp.nCols_;
  vector<double> rowLowerAgg(aggSi.getRowLower(), aggSi.getRowLower() + nAggRowsWCuts);
  vector<double> rowUpperAgg(aggSi.getRowUpper(), aggSi.getRowUpper() + nAggRowsWCuts);
  vector<double> ARvalueAgg(aggSi.getMatrixByRow()->getElements(), 
                                aggSi.getMatrixByRow()->getElements() + nnzAggWCuts);
  vector<int> ARindexAgg(aggSi.getMatrixByRow()->getIndices(), 
                                aggSi.getMatrixByRow()->getIndices() + nnzAggWCuts);
  vector<int> ARstartAgg(aggSi.getMatrixByRow()->getVectorStarts(), 
                                aggSi.getMatrixByRow()->getVectorStarts() + nAggRowsWCuts + 1);
  vector<vector<int> > C = ep.C; 
  // Update info for original LP
  int nRowsAdded = 0;
  for (int i = nAggRows; i < nAggRowsWCuts; ++i){
    double lb = rowLowerAgg[i];
    double ub = rowUpperAgg[i];
    vector<int> indices;
    vector<double> values;
    // Do row wise update (This is easy because cuts are in row format)
    for (int j = ARstartAgg[i]; j < ARstartAgg[i + 1]; ++j){
      int color = ARindexAgg[j];
      double alpha = ARvalueAgg[j];
      vector<int> tempIndices(C[color].begin(), C[color].end());
      indices.insert(indices.end(), tempIndices.begin(), tempIndices.end());
      vector<double> tempValues(tempIndices.size(), alpha);
      values.insert(values.end(), tempValues.begin(), tempValues.end());
      // cout << "one col finished in this row" << endl;
    }
    int numEl = indices.size();
    si.addRow(numEl, &indices[0], &values[0], lb, ub);
  } 
  for (int i = 0; i < si.getNumCols(); ++i){
    si.setInteger(i);
  }
  si.writeLp("testSi");
  // Columns will be automagically updated using clp model class in background
  // Pull new LP info with added cuts in the highest space
  int nRowsWCuts = si.getNumRows();
  int nColsWCuts = si.getNumCols();
  int nnzWCuts = si.getMatrixByCol()->getVectorStarts()[nColsWCuts];
  rowLower.insert(rowLower.end(), si.getRowLower() + nRows, si.getRowLower() + nRowsWCuts);
  rowUpper.insert(rowUpper.end(), si.getRowUpper() + nRows, si.getRowUpper() + nRowsWCuts);
  //Avalue.clear();
  Avalue.assign(si.getMatrixByCol()->getElements(), 
                si.getMatrixByCol()->getElements() + nnzWCuts);
  //Aindex.clear();
  Aindex.assign(si.getMatrixByCol()->getIndices(), 
                si.getMatrixByCol()->getIndices() + nnzWCuts);
  //Astart.clear();
  Astart.assign(si.getMatrixByCol()->getVectorStarts(),
                si.getMatrixByCol()->getVectorStarts() + nColsWCuts + 1);
  ARvalue.insert(ARvalue.end(), si.getMatrixByRow()->getElements() + nnz, 
                 si.getMatrixByRow()->getElements() + nnzWCuts);
  ARindex.insert(ARindex.end(), si.getMatrixByRow()->getIndices() + nnz, 
                 si.getMatrixByRow()->getIndices() + nnzWCuts);
  ARstart.insert(ARstart.end(), si.getMatrixByRow()->getVectorStarts() + nRows + 1, 
                 si.getMatrixByRow()->getVectorStarts() + nRowsWCuts + 1);
  nRows = nRowsWCuts;
  nCols = nColsWCuts;  
  nnz = nnzWCuts;   
}

int main(int argc, const char *argv[]){
    // File name for problem being read in 
    string mpsFileName = argv[1];
    string orbit_file_name = argv[2];
    // Intialize solvers and cut generators
    OsiClpSolverInterface si;
    OsiClpSolverInterface aggSi;
    CglGomory gomoryCuts;
    CglGMI gmiCuts;
    CglTwomir twomirCuts;
    CglMixedIntegerRounding mIRoundingCuts;
    OsiCuts cuts;
    OsiSolverInterface::ApplyCutsReturnCode acRc;

    // Grab info for equitable partition scheme;
    si.readMps(mpsFileName.c_str());
    int nRows = si.getNumRows();
    int _nRows = si.getNumRows();
    int nCols = si.getNumCols();
    int _nCols = si.getNumCols();
    int nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    int _nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    vector<double> colCost(si.getObjCoefficients(), 
                                  si.getObjCoefficients() + nCols);
    vector<double> colLower(si.getColLower(), si.getColLower() + nCols);
    vector<double> colUpper(si.getColUpper(), si.getColUpper() + nCols);
    vector<double> rowLower(si.getRowLower(), si.getRowLower() + nRows);
    vector<double> rowUpper(si.getRowUpper(), si.getRowUpper() + nRows);
    vector<double> Avalue(si.getMatrixByCol()->getElements(), 
                                  si.getMatrixByCol()->getElements() + nnz);
    vector<int> Aindex(si.getMatrixByCol()->getIndices(), 
                                  si.getMatrixByCol()->getIndices() + nnz);
    vector<int> Astart(si.getMatrixByCol()->getVectorStarts(), 
                                  si.getMatrixByCol()->getVectorStarts() + nCols + 1);
    vector<double> ARvalue(si.getMatrixByRow()->getElements(), 
                                  si.getMatrixByRow()->getElements() + nnz);
    vector<int> ARindex(si.getMatrixByRow()->getIndices(), 
                                  si.getMatrixByRow()->getIndices() + nnz);
    vector<int> ARstart(si.getMatrixByRow()->getVectorStarts(), 
                                  si.getMatrixByRow()->getVectorStarts() + nRows + 1);
    orbit_reader orbit_read(orbit_file_name.c_str());
    orbital_partition op = orbit_read.get_orbital_partition();
    // // Pass LP info off to EP algorithm to compute an
    // // EP of the LP
    // EquitablePartition ep(nRows, nCols, colCost, colLower,mv 
    //                         colUpper, rowLower, rowUpper,
    //                         Avalue, Aindex, Astart,
    //                         ARvalue, ARindex, ARstart);
    // // Pass EP to aggregator for first aggregation
    // AggregateLp alp(ep);
    // alp.aggregate();
    
    // // Grab aggregate Lp dimensions
    // int nCols_ = alp.getNumCol();
    // int nRows_ = alp.getNumRowAfterOrderConstraints();
    // // Load problem into a native Coin clpmodel type
    // aggSi.loadProblem(nCols_, nRows_, &alp.getAstart()[0], &alp.getAindex()[0],
    //                     &alp.getAvalue()[0], &alp.getColLower()[0], &alp.getColUpper()[0],
    //                     &alp.getColCost()[0], &alp.getRowLower()[0], &alp.getRowUpper()[0]);
    // for (int i = 0; i < nCols_; ++i){
    //   aggSi.setInteger(i);
    // }
    // aggSi.writeLp("AggBeforeCutsAddedFromCbc");
    // // // Solve initial aggregate
    // double prevObj = numeric_limits<double>::infinity();
    // aggSi.initialSolve();
    // double obj = aggSi.getObjValue();
    // while (fabs(obj - prevObj) > 1e-4){
    //   prevObj = obj;
    //   cuts = OsiCuts();
    //   // Generate l & p cuts for this model
    //   gomoryCuts.generateCuts(aggSi, cuts);
    //   // gmiCuts.generateCuts(aggSi, cuts);
    //   // twomirCuts.generateCuts(aggSi, cuts);
    //   // mIRoundingCuts.generateCuts(aggSi, cuts);
    //   acRc = aggSi.applyCuts(cuts,0.0);
    //   // Print applyCuts return code
    //   cout <<endl <<endl;
    //   cout <<cuts.sizeCuts() <<" cuts were generated" <<endl;
    //   cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<endl;
    //   cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel() 
    //         <<" were inconsistent for this problem" <<endl;
    //   cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<endl;
    //   cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<endl;
    //   cout <<"  " <<acRc.getNumApplied() <<" were applied" <<endl;
    //   cout <<endl <<endl;
    //   aggSi.writeLp("AggAfterCuts");
    //   aggSi.resolve();
    //   obj = aggSi.getObjValue();
    //   iter++;
    //   const char* fName = argv[1];
    //   char iteration[] = "_iter_";
    //   char outPut[100];
    //   // sprintf(outPut, "%s%s%d", fName, iteration, iter);
    //   nCuts += acRc.getNumApplied();
    //   nCutsAdded += acRc.getNumApplied();
    // }
    // update(_nCols, _nRows, _nnz, Astart, Aindex,
    //         Avalue, ARstart, ARindex, ARvalue,
    //         colLower, colUpper, rowLower, rowUpper,
    //         ep, alp, aggSi, si);
    // while(!ep.isDiscrete()){
    //   ep.refine();
    //   alp.updateMasterLpAndEp(ep, _nCols, _nRows, 
    //                     _nnz, Astart, Aindex,
    //                     Avalue, rowLower, rowUpper);
    //   alp.aggregate();
    //   // Grab aggregate Lp dimensions
    //   nCols_ = alp.getNumCol();
    //   nRows_ = alp.getNumRowAfterOrderConstraints();
    //   // Load problem into a native Coin clpmodel type
    //   aggSi = OsiClpSolverInterface();
    //   aggSi.loadProblem(nCols_, nRows_, &alp.getAstart()[0], &alp.getAindex()[0],
    //                       &alp.getAvalue()[0], &alp.getColLower()[0], &alp.getColUpper()[0],
    //                       &alp.getColCost()[0], &alp.getRowLower()[0], &alp.getRowUpper()[0]);
    //   for (int i = 0; i < nCols_; ++i){
    //     aggSi.setInteger(i);
    //   }
    //   aggSi.writeLp("AggBeforeCutsAddedFromCbc");
    //   // Solve 

    //   aggSi.initialSolve();
    //   double obj = aggSi.getObjValue();
    //   while (fabs(obj - prevObj) > 1e-4){
    //     prevObj = obj;
    //     // Generate l & p cuts for this model
    //     cuts = OsiCuts();
    //     gomoryCuts.generateCuts(aggSi, cuts);
    //     // gmiCuts.generateCuts(aggSi, cuts);
    //     // twomirCuts.generateCuts(aggSi, cuts);
    //     // mIRoundingCuts.generateCuts(aggSi, cuts);
    //     acRc = aggSi.applyCuts(cuts,0.0);
    //     // Print applyCuts return code
    //     cout <<endl <<endl;
    //     cout <<cuts.sizeCuts() <<" cuts were generated" <<endl;
    //     cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<endl;
    //     cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel() 
    //           <<" were inconsistent for this problem" <<endl;
    //     cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<endl;
    //     cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<endl;
    //     cout <<"  " <<acRc.getNumApplied() <<" were applied" <<endl;
    //     cout <<endl <<endl;
    //     aggSi.writeLp("AggAfterCuts");
    //     aggSi.resolve();
    //     obj = aggSi.getObjValue();
    //     iter++;
    //     const char* fName = argv[1];
    //     char iteration[] = "_iter_";
    //     char outPut[100];
    //     // sprintf(outPut, "%s%s%d", fName, iteration, iter);
    //     nCuts += acRc.getNumApplied();
    //     nCutsAdded += acRc.getNumApplied();
    //   }
    //   update(_nCols, _nRows, _nnz, Astart, Aindex,
    //           Avalue, ARstart, ARindex, ARvalue,
    //           colLower, colUpper, rowLower, rowUpper,
    //           ep, alp, aggSi, si);
    // }
    // aggSi.initialSolve();
    // cout << "Symmetry Cuts Generated: " << nCuts << endl;
}
