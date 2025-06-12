//
//  GenealogicalNetworkProbStore.h
//  
//
//  Created by Yufeng Wu on 2/5/19.
//  Efficient store and retrievel of gene network probability
//

#ifndef ____GenealogicalNetworkProbStore__
#define ____GenealogicalNetworkProbStore__

#include <vector>
#include <map>
using namespace std;

#include "GenealogicalNetwork.h"

class PhylogenyTreeBasic;
class MarginalTree;
class GenericGeneSpeciesTreeProb;

//***********************************************************************************
// Length grid

class GTBranchLengthGrid
{
public:
    static GTBranchLengthGrid &Instance();
    double SnapToGrid(double brOrig) const;
    double SnapMRToGrid(double mr) const;

private:
    GTBranchLengthGrid();
    void Init();
    
    vector<double> listBrLenGridVals;
    vector<double> listMRGridVals;
};

//***********************************************************************************
// Gene tree prob store

class GTTreeProbStore
{
public:
    static GTTreeProbStore &Instance();
    void AddTreeProb(const MarginalTree &mtree, const PhylogenyTreeBasic &gtree, double prob);
    bool IsProbAlreadyComputed(const MarginalTree &mtree, const PhylogenyTreeBasic &gtree, double &probPre) const;
    
private:
    GTTreeProbStore();
    
    map<string, map<string,double> > mapPreComputedProbForTrees;
};


#endif /* defined(____GenealogicalNetworkProbStore__) */
