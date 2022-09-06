#include "Aggregate.hpp"

AggregateLp::AggregateLp(EquitablePartition& ep){
    // Original Lp info
    nRows = ep.nRows;
    _nRows = nRows;
	nCols = ep.nCols;
    _nCols = nCols;
	numTot = nRows + nCols;
    _numTot = numTot;
	colCost.assign(ep.colCost.begin(), ep.colCost.end());
    colLower.assign(ep.colLower.begin(), ep.colLower.end());
    colUpper.assign(ep.colUpper.begin(), ep.colUpper.end());
    rowLower.assign(ep.rowLower.begin(), ep.rowLower.end());
    rowUpper.assign(ep.rowUpper.begin(), ep.rowUpper.end());
    Avalue.assign(ep.Avalue.begin(), ep.Avalue.end());
    Aindex.assign(ep.Aindex.begin(), ep.Aindex.end());
    AindexP.assign(ep.AindexP.begin(), ep.AindexP.end());
    Astart.assign(ep.Astart.begin(), ep.Astart.end());

    // EP info
    C.assign(ep.C.begin(), ep.C.end());
    Csize.assign(ep.Csize.begin(), ep.Csize.end());
    color.assign(ep.color.begin(), ep.color.end());
}

void AggregateLp::updateMasterLpAndEp(EquitablePartition& ep, int _nC, int _nR,
                            int _nnz, vector<int>& As, vector<int>& Ai,
                            vector<double>& Av, vector<double>& rL, vector<double>& rU){
   // update original Lp info
    _nRows = _nR;
    _nCols = _nC;
	_numTot = _nRows + _nCols;
    _rowLower.assign(rL.begin(), rL.end());
    _rowUpper.assign(rU.begin(), rU.end());
    _Avalue.assign(Av.begin(), Av.end());
    _Aindex.assign(Ai.begin(), Ai.end());
    _Astart.assign(As.begin(), As.end());

    // update EP info
    C.assign(ep.C.begin(), ep.C.end());
    Csize.assign(ep.Csize.begin(), ep.Csize.end());
    color.assign(ep.color.begin(), ep.color.end());
    AindexP.assign(ep.AindexP.begin(), ep.AindexP.end());
    parentPartition.assign(ep.parentPartition.begin(), ep.parentPartition.end());
}

void AggregateLp::clear(){
    Avalue_.clear();
    Astart_.clear();
    Aindex_.clear();
}

void AggregateLp::findDimensions(){
    for (int i = tempNCols_; i < nCols; ++i)
        if (Csize[i])
            nCols_++;
    for (int i = tempNRows_; i < nRows; ++i)
        if (Csize[i + nCols])
            nRows_++;
    colLower_.resize(nCols_);
    colUpper_.resize(nCols_);
    colCost_.resize(nCols_); 
}

void AggregateLp::pairCutWithIndex(){
    cut.clear();
    cutIdx.clear();
    nCuts = nRows_;
    int start = nRows;
    int finish = _nRows;
    for (int i = start; i < finish; ++i){
        cutIdx[i] = nCuts++;
    } 
}

void AggregateLp::pairOrderConstraintWithIndex(){
    orderConstraintIdx.clear();
    nOrdConstraints = nCuts;
    int orderIdx = 0;
    for (int i = 0; i < parentPartition.size(); ++i)
        if (parentPartition[i] != -1)
            orderConstraintIdx[orderIdx++] = nOrdConstraints++;
    rowLower_.resize(nOrdConstraints);
    rowUpper_.resize(nOrdConstraints);
}

void AggregateLp::aggregateColBounds(){
    for (int i = 0; i < nCols_; ++i){
        int rep = C[i].front();
        int orbitSize = Csize[color[rep]];
        colLower_[i] = colLower[rep] * orbitSize;
        colUpper_[i] = colUpper[rep] * orbitSize;
    }
}

void AggregateLp::aggregateRowBounds(){
    for (int i = 0; i < nRows_; ++i){
        int rep = C[i + nCols].front() - nCols;
        int orbitSize = Csize[color[rep + nCols]];
        if (rowLower[rep] == inf)
            rowLower_[i] = inf;
        else
            rowLower_[i] = rowLower[rep] * orbitSize;
        if (rowUpper[rep] == inf)
            rowUpper_[i] = inf;
        else
            rowUpper_[i] = rowUpper[rep] * orbitSize;
    }
    for (int i = nRows_; i < nCuts; ++i){
        if (_rowLower[i - nRows_ + nRows] == inf)
            rowLower_[i] = inf;
        else{
            if (_rowLower[i - nRows_ + nRows] < 0)
                rowLower_[i] = _rowLower[i - nRows_ + nRows];
            else
                rowLower_[i] = _rowLower[i - nRows_ + nRows];  
        }
        if (_rowUpper[i - nRows_ + nRows] == inf)
            rowUpper_[i] = inf;
        else{
            if (_rowUpper[i - nRows_ + nRows] < 0)
                rowUpper_[i] = _rowUpper[i - nRows_ + nRows];
            else
                rowUpper_[i] = _rowUpper[i - nRows_ + nRows];  
        }
    }
    for (int i = nCuts; i < nOrdConstraints; ++i){
        rowLower_[i] = 0;
        rowUpper_[i] = inf;
    }
}

void AggregateLp::aggregateAMatrix(){
	Astart_.push_back(0);
	for (int i = 0; i < nCols_; ++i){
        int rep = C[i].front();
        vector<double> coeff(nRows_, 0);
        for (int j = Astart[rep]; j < Astart[rep + 1]; ++j){
            int rowIdx = AindexP[j] - nCols;
            coeff[rowIdx] += Avalue[j];
        }
        for (int j = 0; j < coeff.size(); ++j){
            if (coeff[j]){
                Avalue_.push_back(coeff[j]);
                Aindex_.push_back(j);
            }
        }
        Astart_.push_back(Avalue_.size());
	}
}

void AggregateLp::addCutsToAggregate(){
    if (!_Astart.size())
        return;
    for (int i = 0; i < nCols_; ++i){
        int rep = C[i].front();
        // vector<double> coeff(nRows_ + _nRows - nRows);
        for (int j = _Astart[rep]; j < _Astart[rep + 1]; ++j){
            if (_Aindex[j] >= nRows){
                int aggColEnd = Astart_[i + 1];
                int cIdx = cutIdx[_Aindex[j]];
                double value = _Avalue[j];
                Aindex_.insert(Aindex_.begin() + aggColEnd, cIdx);
                Avalue_.insert(Avalue_.begin() + aggColEnd, value);
                for (int k = i + 1; k <= nCols_; ++k){
                    Astart_[k]++;
                }
            }
        }
        //std::cout << "column finished for generated cuts" << std::endl;
    }
}

void AggregateLp::addOrderConstraints(){
    if (!_Astart.size()) return;
    int orderIdx = 0;
    for (int i  = 0; i < parentPartition.size(); ++i){
        if (parentPartition[i] != -1){
            int parent = parentPartition[i];
            int child = i;
            // add row to parent node
            int aggColEnd = Astart_[parent + 1];
            int rowIdx = orderConstraintIdx[orderIdx++];
            double value = Csize[child];
            Aindex_.insert(Aindex_.begin() + aggColEnd, rowIdx);
            Avalue_.insert(Avalue_.begin() + aggColEnd, value);
            for (int k = parent + 1; k <= nCols_; ++k){
                Astart_[k]++;
            }
            // add row to childe node
            aggColEnd = Astart_[child + 1];
            value = -1;
            Aindex_.insert(Aindex_.begin() + aggColEnd, rowIdx);
            Avalue_.insert(Avalue_.begin() + aggColEnd, value);
            for (int k = child + 1; k <= nCols_; ++k){
                Astart_[k]++;
            }
        }
    }
}

void AggregateLp::aggregateCostVector(){
    colCost_.assign(nCols_, 0);
	for (int i = 0; i < nCols_; ++i){
        int rep = C[i].front();
        colCost_[i] = colCost[rep];
    }
}

void AggregateLp::aggregate(){
    clear();
    findDimensions();
    pairCutWithIndex();
    pairOrderConstraintWithIndex();
    aggregateColBounds();
    aggregateRowBounds();
    aggregateAMatrix();
    addCutsToAggregate();
    addOrderConstraints();
    aggregateCostVector();
    tempNCols_ = nCols_;
    tempNRows_ = nRows_;
}

int AggregateLp::getNumCol(){
    return nCols_;
}

int AggregateLp::getNumRowAfterCuts(){
    if (!_Astart.size()) return nRows_;
    return nCuts;
}
int AggregateLp::getNumRowAfterOrderConstraints(){
    if (!_Astart.size()) return nRows_;
    return nOrdConstraints;
}

int AggregateLp::getNumRow(){
    return nRows_;
}

vector<double>& AggregateLp::getColUpper(){
    return colUpper_;
}

vector<double>& AggregateLp::getColLower(){
    return colLower_;
}

vector<double>& AggregateLp::getColCost(){
    return colCost_;
}

vector<double>& AggregateLp::getRowUpper(){
    return rowUpper_;
}

vector<double>& AggregateLp::getRowLower(){
    return rowLower_;
}
vector<double>& AggregateLp::getAvalue(){
    return Avalue_;
}

vector<int>& AggregateLp::getAindex(){
    return Aindex_;
}   

vector<int>& AggregateLp::getAstart(){
    return Astart_;
}
