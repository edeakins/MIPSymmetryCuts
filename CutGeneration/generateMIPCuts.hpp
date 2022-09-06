#ifndef GENERATEMIPCUTS_H_
#define GENERATEMIPCUTS_H_

// Necessary coinor headers
#include "ClpSimplex.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include "CoinModel.hpp"
#include "OsiCuts.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
// Combinatorial cuts
#include "CglAllDifferent.hpp"
// #include "CglBKClique.hpp"
#include "CglClique.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
// #include "CglOddWheel.hpp"
#include "CglZeroHalf.hpp"
#include "CglFlowCover.hpp"
// Gomory cuts and variants 
#include "CglGomory.hpp"
#include "CglGMI.hpp"
#include "CglRedSplit.hpp"
#include "CglRedSplit2.hpp"
// Lift and project cuts
#include "CglLiftAndProject.hpp"
#include "CglLandP.hpp"
// Mixed integer rounding cuts
#include "CglMixedIntegerRounding.hpp"
#include "CglMixedIntegerRounding2.hpp" 
#include "CglTwomir.hpp"
#include "CglResidualCapacity.hpp"
// Strengthening cuts
// #include "CglCliqueStrengthening.hpp"
#include "CglDuplicateRow.hpp"
// #include "CglPreprocess.hpp"
#include "CglProbing.hpp"
#include "CglSimpleRounding.hpp"
// Color refinement and aggregation 
#include "EquitablePartition.hpp"
#include "Aggregate.hpp"

#include <iomanip>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <array>

#endif