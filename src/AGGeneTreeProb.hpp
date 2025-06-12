//
//  AGGeneTreeProb.hpp
//  
//
//  Created by Yufeng Wu on 5/15/24.
//

#ifndef AGGeneTreeProb_hpp
#define AGGeneTreeProb_hpp

#include <vector>
#include <map>
#include "GenealogicalNetwork.h"

class PhylogenyTreeBasic;
class MarginalTree;
class GenericGeneSpeciesTreeProb;
class STApproxGeneTreeProb;
class STApproxGeneTreeProb2;
class ApproxGeneTreeProbHeu2;
class SGSTApproxGeneTreeProb;

//***********************************************************************************
// storing computed prob of networks

class AGGeneTreeProbDepot
{
public:
    static AGGeneTreeProbDepot &Instance();
    
    // query whether this tree has been computed or not
    bool QueryProbForNetTree(const std::string &netTreeNW, std::vector<double> &listGTProbs) const;
    
    // store a newly computed tree
    void AppendProb( const std::string &netTreeNW, const std::vector<double> &listGTProbs );
    
    void DumpStats() const;
    
    // clear out
    void Clear();
    
private:
    AGGeneTreeProbDepot() : numTotQuery(0), numQueryCached(0), numStores(0) {}
    
    // for each tree (w/ branch length), prob of each given GT
    std::map<std::string, std::vector<double> > mapNetTreeGTProbs;
    
    // statistics of usage
    int numTotQuery;
    int numQueryCached;
    int numStores;
};

//***********************************************************************************
// Calc prob of genealogical networks: enhanced version

class AGGeneTreeProb2
{
public:
    AGGeneTreeProb2(GenealogicalNetwork &netIn, std::vector<PhylogenyTreeBasic *> &listGeneTreesIn, int threadInd = 0);
    virtual ~AGGeneTreeProb2();
    double CalcProb();
    void UpdateBranch( GenealogicalNetworkBranch *pBrChange, double lenDiff );
    void UpdateWts();
    //void SetMultithread(int nt) { numThreads = nt; }
    
private:
    void Clear();
    void Init();
    void CreateProbCalc();
    int GetIndexForNetcode( const GenealogicalNetworkMTreeCode &ncode );
    double GetCachedProbValFor(int gtIndex, int mtIndex);
    void SetCachedProbValFor(int gtIndex, int mtIndex, double val);
    double CalcProbFromCached();
    double CalcProbMultithread();
    bool IsFresh() const { return fFresh; }
    void SetFresh(bool f ) { fFresh = f; }
    void SnapMargTreeLens();
    void SnapMargTreeLensFor(int gtIndex, int mtIndex, std::vector<double> &listOrigLens);
    void RestoreBrLenFor(int gtIndex, int mtIndex, const std::vector<double> &listOrigLens);
    int GetNumMargTreeEdges() const;
    int GetNumTaxa() const;
    
    GenealogicalNetwork &network;
    std::vector<PhylogenyTreeBasic *> &listGeneTrees;
    //std::vector< STApproxGeneTreeProb *> listGTProbPtrs;
    //std::vector< STApproxGeneTreeProb2 *> listGTProbPtrs;
    std::vector< ApproxGeneTreeProbHeu2 *> listGTProbPtrs;
    //std::vector< SGSTApproxGeneTreeProb *> listGTProbPtrs;
    std::vector<MarginalTree *> listMargTreesInNet;
    std::vector<double> listGenetreeProbWts;
    std::map<GenealogicalNetworkMTreeCode, int> mapNetCodeIndex;
    std::map<GenealogicalNetworkMTreeCode, MarginalTree *> mapNodeToMargTree;
    std::map<GenealogicalNetworkMTreeCode, std::map<GenealogicalNetworkBranch *, int> > mapMargTreeNetBr;    // for each marginal tree, which net branch affects which mtree branch
    std::vector<bool> listMargTreeFresh;             // is marginal tree changed or not
    std::vector<std::vector<double> > listCachedProbs;    // to avoid recomputing
    bool fFresh;
    bool fWeightUnchanged;
    double probLastComputed;
    //int numThreads;
    int threadIndex;
};


#endif /* AGGeneTreeProb_hpp */
