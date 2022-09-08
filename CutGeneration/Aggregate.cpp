#include "Aggregate.hpp"

Aggregate_lp::Aggregate_lp(orbital_partition& op, linear_program& lp){
    // Orbital partition info
    orbits = op;
    // Original Lp info
    nRows = lp.num_row;
    _nRows = nRows;
	nCols = lp.num_col;
    _nCols = nCols;
	numTot = nRows + nCols;
    _numTot = numTot;
	colCost.assign(lp.col_cost.begin(), lp.col_cost.end());
    colLower.assign(lp.col_lower.begin(), lp.col_lower.end());
    colUpper.assign(lp.col_upper.begin(), lp.col_upper.end());
    rowLower.assign(lp.row_lower.begin(), lp.row_lower.end());
    rowUpper.assign(lp.row_upper.begin(), lp.row_upper.end());
    Avalue.assign(lp.a_value.begin(), lp.a_value.end());
    Aindex.assign(lp.a_index.begin(), lp.a_index.end());
    // AindexP.assign(lp.a_indexP.begin(), lp.AindexP.end());
    Astart.assign(lp.a_start.begin(), lp.a_start.end());

    // // EP info
    // C.assign(ep.C.begin(), ep.C.end());
    // Csize.assign(ep.Csize.begin(), ep.Csize.end());
    // color.assign(ep.color.begin(), ep.color.end());
}

void Aggregate_lp::updateMasterLpAndEp(orbital_partition& op, int _nC, int _nR,
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
    // // update EP info
    // C.assign(ep.C.begin(), ep.C.end());
    // Csize.assign(ep.Csize.begin(), ep.Csize.end());
    // color.assign(ep.color.begin(), ep.color.end());
    // AindexP.assign(ep.AindexP.begin(), ep.AindexP.end());
    // parentPartition.assign(ep.parentPartition.begin(), ep.parentPartition.end());
}

void Aggregate_lp::clear(){
    Avalue_.clear();
    Astart_.clear();
    Aindex_.clear();
}

void Aggregate_lp::findDimensions(){
    int i_part, el, el_idx;
    for (i_part = 0; i_part < orbits.orbit_start.size() - 1; ++i_part){
        el_idx = orbits.orbit_start.at(i_part);
        el = orbits.element.at(el_idx);
        if (el < _nCols)
            nCols_++;
        else
            nRows_++;
    }
    colLower_.resize(nCols_);
    colUpper_.resize(nCols_);
    colCost_.resize(nCols_); 
}

void Aggregate_lp::pairCutWithIndex(){
    cut.clear();
    cutIdx.clear();
    nCuts = nRows_;
    int start = nRows;
    int finish = _nRows;
    for (int i = start; i < finish; ++i){
        cutIdx[i] = nCuts++;
    }
    rowLower_.resize(nCuts);
    rowUpper_.resize(nCuts); 
}

void Aggregate_lp::pairOrderConstraintWithIndex(){
    orderConstraintIdx.clear();
    nOrdConstraints = nCuts;
    int orderIdx = 0;
    for (int i = 0; i < parentPartition.size(); ++i)
        if (parentPartition[i] != -1)
            orderConstraintIdx[orderIdx++] = nOrdConstraints++;
    rowLower_.resize(nOrdConstraints);
    rowUpper_.resize(nOrdConstraints);
}

void Aggregate_lp::aggregateColBounds(){
    for (int i = 0; i < nCols_; ++i){
        int el_idx = orbits.orbit_start.at(i);
        int el = orbits.element.at(el_idx);
        int orbitSize = orbits.orbit_start.at(i + 1) - el_idx;
        colLower_[i] = colLower[el] * orbitSize;
        colUpper_[i] = colUpper[el] * orbitSize;
    }
}

void Aggregate_lp::aggregateRowBounds(){
    for (int i = 0; i < nRows_; ++i){
        int el_idx = orbits.orbit_start.at(i + nCols_);
        int el = orbits.element.at(el_idx) - nCols;
        int orbitSize = orbits.orbit_start.at(i + nCols_ + 1) - el_idx;
        if (rowLower[el] == inf)
            rowLower_[i] = inf;
        else
            rowLower_[i] = rowLower[el] * orbitSize;
        if (rowUpper[el] == inf)
            rowUpper_[i] = inf;
        else
            rowUpper_[i] = rowUpper[el] * orbitSize;
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

void Aggregate_lp::aggregateAMatrix(){
	Astart_.push_back(0);
	for (int i = 0; i < nCols_; ++i){
        int el_idx = orbits.orbit_start.at(i);
        int el = orbits.element.at(el_idx);
        vector<double> coeff(nRows_, 0);
        for (int j = Astart[el]; j < Astart[el + 1]; ++j){
            int a_idx = Aindex.at(j);
            int rowIdx = orbits.orbit.at(a_idx + nCols) - nCols_;
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

void Aggregate_lp::addCutsToAggregate(){
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

void Aggregate_lp::addOrderConstraints(){
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

void Aggregate_lp::aggregateCostVector(){
    colCost_.assign(nCols_, 0);
	for (int i = 0; i < nCols_; ++i){
        int el_idx = orbits.orbit_start.at(i);
        int el = orbits.element.at(el_idx);
        colCost_[i] = colCost[el];
    }
}

void Aggregate_lp::aggregate(){
    clear();
    findDimensions();
    pairCutWithIndex();
    // pairOrderConstraintWithIndex();
    aggregateColBounds();
    aggregateRowBounds();
    aggregateAMatrix();
    addCutsToAggregate();
    // addOrderConstraints();
    aggregateCostVector();
    tempNCols_ = nCols_;
    tempNRows_ = nRows_;
}

int Aggregate_lp::getNumCol(){
    return nCols_;
}

int Aggregate_lp::getNumRowAfterCuts(){
    if (!_Astart.size()) return nRows_;
    return nCuts;
}
int Aggregate_lp::getNumRowAfterOrderConstraints(){
    if (!_Astart.size()) return nRows_;
    return nOrdConstraints;
}

int Aggregate_lp::getNumRow(){
    return nRows_;
}

vector<double>& Aggregate_lp::getColUpper(){
    return colUpper_;
}

vector<double>& Aggregate_lp::getColLower(){
    return colLower_;
}

vector<double>& Aggregate_lp::getColCost(){
    return colCost_;
}

vector<double>& Aggregate_lp::getRowUpper(){
    return rowUpper_;
}

vector<double>& Aggregate_lp::getRowLower(){
    return rowLower_;
}
vector<double>& Aggregate_lp::getAvalue(){
    return Avalue_;
}

vector<int>& Aggregate_lp::getAindex(){
    return Aindex_;
}   

vector<int>& Aggregate_lp::getAstart(){
    return Astart_;
}
