#include "generateMIPCuts.hpp"
using namespace std;

int iter = 0;
int nCuts = 0;
int nCutsAdded = 0;

void update(OsiClpSolverInterface& aggSi, OsiClpSolverInterface& si, orbital_partition& op){
    // Grab info we will need for update
    int num_aggregate_cols = aggSi.getNumCols();
    int num_aggregate_rows = aggSi.getNumRows() - nCutsAdded;
    int nnz_cut_start = aggSi.getMatrixByRow()->getVectorStarts()[num_aggregate_rows];
    int nnz_aggregate = aggSi.getMatrixByRow()->getVectorStarts()[aggSi.getNumRows()];
    vector<double> rowLowerAgg(aggSi.getRowLower() + num_aggregate_rows, 
                                aggSi.getRowLower() + aggSi.getNumRows());
    vector<double> rowUpperAgg(aggSi.getRowUpper() + num_aggregate_rows, 
                                aggSi.getRowUpper() + aggSi.getNumRows());
    vector<double> rowRhsAgg(aggSi.getRightHandSide() + num_aggregate_rows,
                             aggSi.getRightHandSide() + aggSi.getNumRows());
    vector<char> rowSenseAgg(aggSi.getRowSense() + num_aggregate_rows,
                            aggSi.getRowSense() + aggSi.getNumRows());
    vector<double> ARvalueAgg(aggSi.getMatrixByRow()->getElements() + nnz_cut_start, 
                                aggSi.getMatrixByRow()->getElements() + nnz_aggregate);
    vector<int> ARindexAgg(aggSi.getMatrixByRow()->getIndices() + nnz_cut_start, 
                                aggSi.getMatrixByRow()->getIndices() + nnz_aggregate);
    vector<int> ARstartAgg(aggSi.getMatrixByRow()->getVectorStarts() + num_aggregate_rows, 
                                aggSi.getMatrixByRow()->getVectorStarts() + aggSi.getNumRows() + 1); 
    int offset = ARstartAgg.at(0);
    ofstream write_cuts("./CutGeneration/cut_txt_files/cuts.txt");
    // int iter = std::ceil(nCutsAdded/2.0);
    int iter = 10; 
    if (nCutsAdded < 10) iter = nCutsAdded;
    for (int i = 0; i < iter; ++i){
    // for (int i = 0; i < 1; ++i){
        double rhs = rowRhsAgg[i];
        double lb = rowLowerAgg[i];
        double ub = rowUpperAgg[i];
        char sense = rowSenseAgg[i];
        vector<int> indices;
        vector<double> values;
        // Do row wise update (This is easy because cuts are in row format)
        for (int j = ARstartAgg[i]; j < ARstartAgg[i + 1]; ++j){
            int orbit = ARindexAgg[j - offset];
            int start = op.orbit_start.at(orbit);
            int end = op.orbit_start.at(orbit + 1);
            double alpha = ARvalueAgg[j - offset];
            vector<int> tempIndices(op.element.begin() + start, op.element.begin() + end);
            indices.insert(indices.end(), tempIndices.begin(), tempIndices.end());
            vector<double> tempValues(tempIndices.size(), alpha);
            values.insert(values.end(), tempValues.begin(), tempValues.end());
            // cout << "one col finished in this row" << endl;
        }
        int numEl = indices.size();
        vector<double> expanded_values(si.getNumCols(), 0);
        for (int i = 0; i < indices.size(); ++i){
            int idx = indices.at(i);
            double val = values.at(i);
            expanded_values.at(idx) = val;
        }
        // double rhs = std::min(std::fabs(lb), std::fabs(ub));
        for (int i = 0; i < si.getNumCols(); ++i)
            write_cuts << expanded_values.at(i) << ",";
        write_cuts << sense << "," << rhs << "\n";
        si.addRow(numEl, &indices[0], &values[0], lb, ub);
    }   
    write_cuts.close();        
    // si.writeLp("./CutGeneration/full_lps_with_cuts/full_after_cuts");    
    si.initialSolve();           
    nCutsAdded = 0;
    // std::cout << "break point" << std::endl;
//   int nAggRowsWCuts = alp.getNumRowAfterOrderConstraints() + nCutsAdded;
//   nCutsAdded = 0;
//   int nAggColsWCuts = aggSi.getNumCols();
//   int nnzAggWCuts = aggSi.getMatrixByRow()->getVectorStarts()[nAggRowsWCuts];
//   int nAggRows = alp.getNumRowAfterOrderConstraints();
//   int nAggCols = alp.nCols_;
//   
//   // Update info for original LP
//   int nRowsAdded = 0;

//   for (int i = 0; i < si.getNumCols(); ++i){
//     si.setInteger(i);
//   }
//   si.writeLp("testSi");   
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
    CglZeroHalf zeroHalf;
    CglMixedIntegerRounding mIRoundingCuts;
    OsiCuts cuts;
    OsiSolverInterface::ApplyCutsReturnCode acRc;
    // Lp to house original lp instance
    linear_program lp;

    // Grab info for equitable partition scheme;
    // si.setHintParam(OsiDoPresolveInInitial, false, OsiForceDo);
    // si.setHintParam(OsiDoPresolveInResolve, false, OsiForceDo);
    si.readMps(mpsFileName.c_str());
    int nRows = si.getNumRows();
    int _nRows = si.getNumRows();
    int nCols = si.getNumCols();
    int _nCols = si.getNumCols();
    int nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    int _nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    lp.num_col = si.getNumCols();
    lp.num_row = si.getNumRows();
    lp.num_total = lp.num_col + lp.num_row;
    lp.col_cost.assign(si.getObjCoefficients(), 
                                  si.getObjCoefficients() + nCols);
    lp.col_lower.assign(si.getColLower(), si.getColLower() + nCols);
    lp.col_upper.assign(si.getColUpper(), si.getColUpper() + nCols);
    lp.row_lower.assign(si.getRowLower(), si.getRowLower() + nRows);
    lp.row_upper.assign(si.getRowUpper(), si.getRowUpper() + nRows);
    lp.a_value.assign(si.getMatrixByCol()->getElements(), 
                                  si.getMatrixByCol()->getElements() + nnz);
    lp.a_index.assign(si.getMatrixByCol()->getIndices(), 
                                  si.getMatrixByCol()->getIndices() + nnz);
    lp.a_start.assign(si.getMatrixByCol()->getVectorStarts(), 
                                  si.getMatrixByCol()->getVectorStarts() + nCols + 1);
    lp.ar_value.assign(si.getMatrixByRow()->getElements(), 
                                  si.getMatrixByRow()->getElements() + nnz);
    lp.ar_index.assign(si.getMatrixByRow()->getIndices(), 
                                  si.getMatrixByRow()->getIndices() + nnz);
    lp.ar_start.assign(si.getMatrixByRow()->getVectorStarts(), 
                                  si.getMatrixByRow()->getVectorStarts() + nRows + 1);
    // Get orbital partition from file
    orbit_reader orbit_read(orbit_file_name.c_str(), lp.num_total);
    orbital_partition op = orbit_read.get_orbital_partition();
    // Pass the orbital partition of the variables and constraints to the lp aggregator
    // This is a revision from the first time I wrote this code.  Since this implementation will 
    // be more general we will need to use the symmetry group generators and orbital partitions
    // and not fractional symmetries and equitable partitions.
    Aggregate_lp alp(op, lp);
    alp.aggregate();
    // Grab aggregate Lp dimensions
    int nCols_ = alp.getNumCol();
    int nRows_ = alp.getNumRowAfterOrderConstraints();
    // Load problem into a native Coin clpmodel type
    // aggSi.setHintParam(OsiDoPresolveInInitial, false, OsiForceDo);
    // aggSi.setHintParam(OsiDoPresolveInResolve, false, OsiForceDo);
    aggSi.loadProblem(nCols_, nRows_, &alp.getAstart()[0], &alp.getAindex()[0],
                        &alp.getAvalue()[0], &alp.getColLower()[0], &alp.getColUpper()[0],
                        &alp.getColCost()[0], &alp.getRowLower()[0], &alp.getRowUpper()[0]);
    for (int i = 0; i < nCols_; ++i){
      aggSi.setInteger(i);
    }
    // aggSi.writeLp("./CutGeneration/aggregate_lps/agg_before_cuts");
    // Solve initial aggregate
    double prevObj = numeric_limits<double>::infinity();
    aggSi.initialSolve();
    double obj = aggSi.getObjValue();
    cuts = OsiCuts();
    gmiCuts.generateCuts(aggSi, cuts);
    acRc = aggSi.applyCuts(cuts, 0.0);
    cout << cuts.sizeCuts() << " cuts were generated" <<endl;
    // aggSi.writeLp("./CutGeneration/aggregate_lps_with_cuts/agg_after_cuts");
    aggSi.resolve();
    nCutsAdded += acRc.getNumApplied();
    update(aggSi, si, op);
    
    
    // while (fabs(obj - prevObj) > 1e-4){
    //   prevObj = obj;
    //   cuts = OsiCuts();
    //   // Generate l & p cuts for this model
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
