//
//  GenealogicalNetworkProbStore.cpp
//  
//
//  Created by Yufeng Wu on 2/5/19.
//
//

#include "GenealogicalNetworkProbStore.h"
//#include "GenealogicalNetworkProb.h"
#include "GenealogicalNetwork.h"
#include "PhylogenyTreeBasic.h"
#include "MarginalTree.h"
#include "Utils3.h"
#include "Utils4.h"
#include <cmath>
#include "UtilsNumerical.h"

const int DEF_SIZE_BR_LEN_GRID = 43;
const double LIST_DEF_BR_LEN_GRID_VALS[DEF_SIZE_BR_LEN_GRID] = {0.0001, 0.0005, 0.001,0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};
const int DEF_SIZE_MR_GRID = 13;
const double LIST_DEF_MR_GRID_VALS[DEF_SIZE_MR_GRID] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};

//***********************************************************************************
// Length grid

GTBranchLengthGrid & GTBranchLengthGrid :: Instance()
{
    static GTBranchLengthGrid instance;
    return instance;
}

GTBranchLengthGrid :: GTBranchLengthGrid()
{
    Init();
}

double GTBranchLengthGrid :: SnapToGrid(double brOrig) const
{
return brOrig;
    int posClose = GetClosestTo(listBrLenGridVals, brOrig);
    YW_ASSERT_INFO( posClose >= 0, "Fail to find" );
    return listBrLenGridVals[posClose];
}

double GTBranchLengthGrid :: SnapMRToGrid(double mr) const
{
return mr;
    int posClose = GetClosestTo(listMRGridVals, mr);
    YW_ASSERT_INFO( posClose >= 0, "Fail to find" );
    return listMRGridVals[posClose];
}

void GTBranchLengthGrid :: Init()
{
    listBrLenGridVals.clear();
    for(int i=0; i<DEF_SIZE_BR_LEN_GRID; ++i)
    {
        listBrLenGridVals.push_back(LIST_DEF_BR_LEN_GRID_VALS[i]);
    }
    listMRGridVals.clear();
    for(int i=0; i<DEF_SIZE_MR_GRID; ++i)
    {
        listMRGridVals.push_back( LIST_DEF_MR_GRID_VALS[i] );
    }
}

//***********************************************************************************
// Gene tree prob store


GTTreeProbStore & GTTreeProbStore :: Instance()
{
    static GTTreeProbStore instance;
    return instance;
}

void GTTreeProbStore :: AddTreeProb(const MarginalTree &mtree, const PhylogenyTreeBasic &gtree, double prob)
{
return;
    //
    string strmTree = mtree.GetNewickSorted(true);
    string strgtree;
    PhylogenyTreeBasic &gtreeInst = const_cast<PhylogenyTreeBasic &>(gtree);
    gtreeInst.ConsNewickSorted(strgtree);
    mapPreComputedProbForTrees[strmTree][strgtree] = prob;
}

bool GTTreeProbStore :: IsProbAlreadyComputed(const MarginalTree &mtree, const PhylogenyTreeBasic &gtree, double &probPre) const
{
return false;
    //
    string strmTree = mtree.GetNewickSorted(true);
    string strgtree;
    PhylogenyTreeBasic &gtreeInst = const_cast<PhylogenyTreeBasic &>(gtree);
    gtreeInst.ConsNewickSorted(strgtree);
    map<string, map<string,double> > :: const_iterator it = mapPreComputedProbForTrees.find(strmTree);
    if( it == mapPreComputedProbForTrees.end() )
    {
        return false;
    }
    else
    {
        map<string,double> :: const_iterator it2 = it->second.find(strgtree);
        if( it2 == it->second.end() )
        {
            return false;
        }
        probPre = it2->second;
        return true;
    }
}

GTTreeProbStore :: GTTreeProbStore()
{
}


