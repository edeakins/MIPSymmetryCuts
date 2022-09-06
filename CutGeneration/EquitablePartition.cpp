#include "EquitablePartition.hpp"

using namespace std;

EquitablePartition::EquitablePartition( int nR,  int nC,  vector<double>& cC,  vector<double>& cL,
                                             vector<double>& cU,  vector<double>& rL,  vector<double>& rU,
                                             vector<double>& Av,  vector<int>& Ai,  vector<int>& As,
                                             vector<double>& ARv,  vector<int>& ARi,  vector<int>& ARs){
    // Grab LP data
    nRows = nR;
    nCols = nC;
    numTot = nR + nC;
    colCost.assign(cC.begin(), cC.end());
    colLower.assign(cL.begin(), cL.end());
    colUpper.assign(cU.begin(), cU.end());
    rowLower.assign(rL.begin(), rL.end());
    rowUpper.assign(rU.begin(), rU.end());
    Avalue.assign(Av.begin(), Av.end());
    Aindex.assign(Ai.begin(), Ai.end());
    Astart.assign(As.begin(), As.end());
    ARvalue.assign(ARv.begin(), ARv.end());
    ARindex.assign(ARi.begin(), ARi.end());
    ARstart.assign(ARs.begin(), ARs.end());

    // Pre allocate arrays for EP algorithm
	AindexP.resize(Aindex.size());
    AvaluePos.resize(Avalue.size());
    ARvaluePos.resize(ARvalue.size());
    SCheck.assign(numTot, false);
	mincdeg.assign(numTot, 0);
	maxcdeg.assign(numTot, 0);
	cdeg.assign(numTot, 0);
    isAdj.assign(numTot, false);
    numEdges.assign(numTot, 0);
	isolates.assign(numTot, false);
    color.assign(numTot, 0);
	Csize.assign(numTot, 0);
	Asize.assign(numTot, 0);
    C.resize(numTot);
    A.resize(numTot);
	handleNegatives();
    initRefinement();
}

void EquitablePartition::updateForCuts(int nR,  int nC,  vector<double>& cC,  vector<double>& cL,
                                    	vector<double>& cU,  vector<double>& rL,  vector<double>& rU,
                                        vector<double>& Av,  vector<int>& Ai,  vector<int>& As,
                                        vector<double>& ARv,  vector<int>& ARi,  vector<int>& ARs){
	// Original Lp info but edited for cuts
    nCols = nC;
    numTot = nR + nC;
    colCost.assign(cC.begin(), cC.end());
    colLower.assign(cL.begin(), cL.end());
    colUpper.assign(cU.begin(), cU.end());
    rowLower.assign(rL.begin(), rL.end());
    rowUpper.assign(rU.begin(), rU.end());
    Avalue.assign(Av.begin(), Av.end());
    Aindex.assign(Ai.begin(), Ai.end());
    Astart.assign(As.begin(), As.end());
    ARvalue.assign(ARv.begin(), ARv.end());
    ARindex.assign(ARi.begin(), ARi.end());
    ARstart.assign(ARs.begin(), ARs.end());
	// Make color scheme edits for new cut rows
	// These rows are singletons and therefore do not effect anything 
	// other than that theyu exist.
	C.resize(numTot);
	Csize.resize(numTot);
	A.resize(numTot);
	Asize.resize(numTot);
	color.resize(numTot);
	AindexP.resize(Avalue.size());
	for (int i = nRows; i < nR; ++i){
		C[cCol].push_back(nRows + nCols);
		Csize[cCol]++;
		color[nRows + nCols] = cCol;
		cCol++;
		nRows++;
	}
}

void EquitablePartition::initRefinement(){
    int i,j;
	int numParts = 0;
	int varColor = 0;
	int conColor = nCols;
	set<double> Rhs;
	vector<double> Rhs_;
	set<tuple<double, double, double> > rhs;
	set<tuple<double, double, double, double> > objBounds;
    initialParts.resize(numTot);

    // Partition the rows
	pair<set<tuple<double, double, double> >::iterator, bool> retRow;
	int rowSum = 0;
	for (i = ARstart[0]; i < ARstart[1]; ++i)
		rowSum += ARvaluePos[i];
	rhs.insert(make_tuple(rowLower[0], rowUpper[0], rowSum));
	initialParts[nCols] = conColor;
	S.push(conColor);
	SCheck[conColor] = true;
	rowSum = 0;
    for (i = 1; i < nRows; ++i){
		for (j = ARstart[i]; j < ARstart[i + 1]; ++j)
			rowSum += ARvaluePos[j];
		retRow = rhs.insert(make_tuple(rowLower[i], rowUpper[i], rowSum));
		if (retRow.second){
			initialParts[i + nCols] = ++conColor;
			S.push(conColor);
			SCheck[conColor] = true;
		}
		else{
			initialParts[i + nCols] = conColor;
		}
		rowSum = 0;
	}
	conColor++;

    // Partition the columns
	pair<set<tuple<double, double, double, double> >::iterator, bool> retCol;
	int colSum = 0;
	for (i = Astart[0]; i < Astart[1]; ++i){
		colSum += AvaluePos[i];
	}
	objBounds.insert(make_tuple(colCost[0], colLower[0], colUpper[0], colSum));
	initialParts[0] = varColor;
	S.push(varColor);
	SCheck[varColor] = true;
	colSum = 0;
    for (i = 1; i < nCols; ++i){
		for (j = Astart[i]; j < Astart[i + 1]; ++j)
			colSum += AvaluePos[j];
		retCol = objBounds.insert(make_tuple(colCost[i], colLower[i], colUpper[i], colSum));
		if (retCol.second){
			initialParts[i] = ++varColor;
			S.push(varColor);
			SCheck[varColor] = true;;
		}
		else{
			initialParts[i] = varColor;
		}
		colSum = 0;
	}
	varColor++;

    // Define number of colors that have been used so far
	vCol = varColor;
	cCol = conColor;

	// Store initial refinement
	for (i = 0; i < numTot; ++i){
		C[initialParts[i]].push_back(i);
		Csize[initialParts[i]]++;
		color[i] = initialParts[i];
	}
	refine();
}

void EquitablePartition::refine(){
	if (refinements) findTarget();
	refinements++;
	double weight;
	colorsAdj.clear();
	colorsToSplit.clear();
	while (!S.empty()){
		int r = S.top();
		S.pop();
		SCheck[r] = false;
		bool var = r < nCols ? 1 : 0;
		for (int i = 0; i < Csize[r]; ++i){
			int v = C[r][i];
			if (var){
				for (int j = Astart[v]; j < Astart[v + 1]; ++j){
					int w = Aindex[j] + nCols;
					weight = AvaluePos[j];
					cdeg[w] += weight;
					numEdges[w]++;
					if (numEdges[w] == 1){
						A[color[w]].push_back(w);
						Asize[color[w]]++;
					}
					if (!isAdj[color[w]]){
						isAdj[color[w]] = true;
						colorsAdj.push_back(color[w]);
					}
					if (cdeg[w] > maxcdeg[color[w]])
						maxcdeg[color[w]] = cdeg[w];			
				}
			}
			else{
				v -= nCols;
				for (int j = ARstart[v]; j < ARstart[v + 1]; ++j){
					int w = ARindex[j];
					weight = ARvaluePos[j];
					cdeg[w] += weight;
					numEdges[w]++;
					if (numEdges[w] == 1){
						A[color[w]].push_back(w);
						Asize[color[w]]++;
					} 
					if (!isAdj[color[w]]){
						isAdj[color[w]] = true;
						colorsAdj.push_back(color[w]);
					}
					if (cdeg[w] > maxcdeg[color[w]])
						maxcdeg[color[w]] = cdeg[w];			
				}
			}
		}
		for (int i = 0; i < colorsAdj.size(); ++i){
			int c = colorsAdj[i];
			if (Csize[c] != Asize[c])
				mincdeg[c] = 0;
			else{
				mincdeg[c] = maxcdeg[c];
				for (int j = 0; j < Asize[c]; ++j){
					int w = A[c][j];
					if (cdeg[w] < mincdeg[c]) mincdeg[c] = cdeg[w];
				}
			}
		}
		colorsToSplit.clear();
		for (int i = 0; i < colorsAdj.size(); ++i){
			int c = colorsAdj[i];
			if (mincdeg[c] < maxcdeg[c])
				colorsToSplit.push_back(c);
		}
		sort(colorsToSplit.begin(), colorsToSplit.end());
		for (int i = 0; i < colorsToSplit.size(); ++i){
			int s = colorsToSplit[i];
			splitColor(s);
		}
		for (int i = 0; i < colorsAdj.size(); ++i){
			int c = colorsAdj[i];
			for (int j = 0; j < Asize[c]; ++j){
				int w = A[c][j];
				cdeg[w] = 0;
				numEdges[w] = 0;
			}
			maxcdeg[c] = 0;
			isAdj[c] = false;
			A[c].clear();
			Asize[c] = 0;
		}
		colorsAdj.clear();
	}			
	// cout << "\n\n Partition \n\n" << endl;
	for (int i = 0; i < C.size(); ++i){
		if (Csize[i] == 1) 
			isolates[C[i].front()] = true;
	}
	// for (int i = 0; i < coeff.size(); ++i){
	// 	cout << "node: " << i << ": ";
	// 	for (int j = 0; j < coeff[i].size(); ++j){
	// 		cout << coeff[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	packVectors();
}

void EquitablePartition::splitColor(int s){
	bool var = (s < nCols) ? true : false;
	set<double> cdegCopy;
	vector<int> colorFreq(numTot, 0);
	map<double, int> degSumColor;
	pair<map<double, int>::iterator, bool> ret;
	degSumColor.insert(pair<double, int>(mincdeg[s], s));
	colorFreq[0] = Csize[s] - Asize[s];
	for (int i = 0; i < Asize[s]; ++i){
		int w = A[s][i];
		if (var){
			ret = degSumColor.insert(pair<double, int>(cdeg[w], vCol));
			if (ret.second){
				if (parentPartition[s] != -1)
					parentPartition[vCol] = parentPartition[s];
				vCol++;
			}
		}
		else{
			ret = degSumColor.insert(pair<double, int>(cdeg[w], cCol));
			cCol += ret.second == 1;
		}
		colorFreq[cdeg[w]]++;
	}
	int b = distance(colorFreq.begin(), max_element(colorFreq.begin(), colorFreq.end()));
	int instack = (SCheck[s]) ? 1 : 0;
	for(map<double, int>::iterator it = degSumColor.begin(); it != degSumColor.end(); ++it){
		//coeff[r].push_back(it->second);
		if (it->first == mincdeg[s]){
			if (!instack && it->first != b){
				S.push(it->second);
				SCheck[s] = true;
			}
		}
		else{
			if (it->first != b){
				S.push(it->second);
				SCheck[it->second] = true;
			}
		}
	}
	for (int i = 0; i < Asize[s]; ++i){
		int w = A[s][i];
		if (degSumColor[cdeg[w]] != s){
			C[s].erase(remove(C[s].begin(), C[s].end(), w), C[s].end());
			Csize[s]--;
			C[degSumColor[cdeg[w]]].push_back(w);
			Csize[degSumColor[cdeg[w]]]++;
			color[w] = degSumColor[cdeg[w]];
		}
	}
}

void EquitablePartition::handleNegatives(){
    vector<double>::iterator min = min_element(Avalue.begin(), Avalue.end());
    if (*min < 0){
        for (int i = 0; i < Avalue.size(); ++i){
            AvaluePos[i] = Avalue[i] + (-*min) + 1;
            ARvaluePos[i] = ARvalue[i] + (-*min) + 1;
        }
    }
    else{
        for (int i = 0; i < Avalue.size(); ++i){
            AvaluePos[i] = Avalue[i];
            ARvaluePos[i] = ARvalue[i];
        }
    }
}

void EquitablePartition::findTarget(){
    for (int col = 0; col < numTot; ++col){
		if (Csize[col] > 1){
			// int isolated = i;
			isolate(col);
			return;
		}
	}
}

void EquitablePartition::isolate(int s){
	parentPartition.assign(nCols, -1);
	vector<int> temp1;
	vector<int> temp2;
	for (int i = 0; i < C[s].size(); ++i)
		if (!i) {temp1.push_back(C[s][i]); color[C[s][i]] = s;}
		else {temp2.push_back(C[s][i]); color[temp2[i - 1]] = vCol; Csize[vCol]++;}
	C[s] = temp1;
	C[vCol] = temp2;
	Csize[s] = 1;
	for (int i = 0; i < C.size(); ++i){
		if (Csize[i] == 1)
			isolates[C[i][0]] = true;
	}
	parentPartition[vCol] = s;
	SCheck[vCol] = true;
	S.push(vCol);
	vCol++;
}

bool EquitablePartition::isDiscrete(){
	for (int i = 0; i < numTot; ++i){
		if (Csize[i] > 1)
			return false;
	}
	return true;
}

void EquitablePartition::packVectors(){
	for (int i = 0; i < Aindex.size(); ++i)
		AindexP[i] = color[Aindex[i] + nCols];
}