//
//  AGHeuristicSearch.hpp
//  
//
//  Created by Yufeng Wu on 2/4/25.
//

#ifndef AGHeuristicSearch_hpp
#define AGHeuristicSearch_hpp

#include <vector>
#include <set>
#include <map>

#include "GenealogicalNetwork.h"

class GenealogicalNetwork;
class GenealogicalNetworkBranch;
class PhylogenyTreeBasic;

//***********************************************************************************
// Evaluate how a gene tree fits AG

class AGHeuristicSearchEval
{
public:
    AGHeuristicSearchEval(GenealogicalNetwork &agCurr);
    
    double EvalTree( PhylogenyTreeBasic *pTree );
    
private:
    // evaluate the AG
    void ProcessAG();
    
    GenealogicalNetwork &agCurr;
    // <AG branch, <set of taxa, score/prob> >
    std::map<GenealogicalNetworkBranch *,std::map<std::set<int>, double> > mapAGInfo;
};

//***********************************************************************************
// Evaluate how a gene tree fits AG based on trees contained in AG

class AGHeuristicSearchTreeBaseEval
{
public:
    AGHeuristicSearchTreeBaseEval(GenealogicalNetwork &agCurr);
    virtual ~AGHeuristicSearchTreeBaseEval();
    
    double EvalTree( PhylogenyTreeBasic *pTree );
    
private:
    // evaluate the AG
    void ProcessAG();
    
    GenealogicalNetwork &agCurr;
    
    // index: tree clade index, <set of taxa> >
    std::vector< std::set<std::set<int> > > listContainedAGTreeInfo;
};

//***********************************************************************************
// construct network from trees
class AGHeuristicBuilder
{
public:
    AGHeuristicBuilder(std::vector<PhylogenyTreeBasic *> &listGeneTreePtrs, int og);
    virtual ~AGHeuristicBuilder();
    void Build(int numMixEvts);
    double GetCurrScore() const { return scoreCurrNet; }
    GenealogicalNetwork *GetNet() const { return pnetCurr; }
    
private:
    // init
    void Init();
    // score the curr net
    double ScoreNet(GenealogicalNetwork &net);
    // search for one best net by making one tree branch to be reticulate and then regraft this edge to somewhere else
    double FindBestNetAddOneMix();
    // Find best neighboring net by local search
    double FindOptNetLocalSearch();
    
    std::vector<PhylogenyTreeBasic *> &listGeneTreePtrs;
    int og;
    GenealogicalNetwork *pnetCurr;
    double scoreCurrNet;
};

#endif /* AGHeuristicSearch_hpp */
