//
//  AGFastCoalBuilder.hpp
//  
//
//  Created by Yufeng Wu on 3/18/24.
//

#ifndef AGFastCoalBuilder_hpp
#define AGFastCoalBuilder_hpp

#include "PhylogenyTreeBasic.h"
#include "GenealogicalNetwork.h"
#include "MarginalTree.h"
#include <thread>
#include <mutex>

//////////////////////////////////////////////////////////////////////////////////
// Store processed network topologies

class AGProcessedNetsDepot
{
public:
    static AGProcessedNetsDepot& Instance();
    bool IsNetProcessed( GenealogicalNetwork *pNet);
    void MarkNetProcessed(GenealogicalNetwork *pNet);
    void SkipNet(GenealogicalNetwork *pNet);
    void DumpStats() const;
    
private:
    AGProcessedNetsDepot();
    void ConsNetSignature(GenealogicalNetwork *pNet, map<int,int> &sigNet);
    
    // next id to use
    //int idNextToUse;
    
    // for each string, assign an id
    set<map<int,int> > mapNettoId;
    
    // skipped counts
    int numSkipped;
    
    // for concurrency
    std::mutex mut;
};

//***********************************************************************************
// Fast inference of admixture graph from gene trees

class AGFastCoalBuilder
{
public:
    AGFastCoalBuilder(vector<PhylogenyTreeBasic *> &listGeneTreesIn, GenealogicalNetwork &networkInitIn, TaxaMapper &mapperTaxaIdsIn);
    ~AGFastCoalBuilder();
    double Search(int numMixNodes=1);
    double Search2(int numMixNodes=1);
    void InfMixPops( set<string> &setMixTaxaUser );
    GenealogicalNetwork *GetBestNet();
    void SetMaxNumMixNodes(int mn) { maxNumMixNodes = mn; }
    int GetMaxNumMixNodes() const { return maxNumMixNodes; }
    void SetHeuSearch(bool fHeu) { fHeuSearch = fHeu; }
    void MarkNetProcessed(GenealogicalNetwork *pNet) const;
    void SetInitNetOnly(bool f) { fInitNetOnly = f; }
    void SetInitMixTaxa(const set<int> &setInitMix) { taxaMixInit = setInitMix; }
    static void SetOutgroup(int og);
    static int GetOutgroup();
    
private:
    // implementation
    void Init();
    void FindMixTaxaImp( vector<int> &setMixTaxa );
    int FindOneMixTaxonFrom( const set<int> &setTaxaCurr );
    double ScoreForNetwork(const vector<PhylogenyTreeBasic *> &listTrees, GenealogicalNetwork *pNetCurr);
    double ScoreForNetwork2(const vector<PhylogenyTreeBasic *> &listTrees, GenealogicalNetwork *pNetCurr);
    double ScoreForNetworkMDC(const vector<PhylogenyTreeBasic *> &listTrees, GenealogicalNetwork *pNetCurr);
    double ExploreNgbrs();
    bool IsNetProcBefore(GenealogicalNetwork *pNet) const;
    
    // helpers
    void SetCurLogProb(double lp) { logprobBestCurr = lp; }
    double GetCurrBestLogprob() const { return logprobBestCurr; }
    void GetAllTaxa(set<int> &setTaxa) const;
    void GetReducedGenetreesForSubsetTaxa( const set<int> &setTaxaSub, vector<PhylogenyTreeBasic *> &listGeneTreesSub ) const;
    bool IsNetworkOGGood( const GenealogicalNetwork &netTest );
    static int ScoreMDCGTandST( PhylogenyTreeBasic *pGeneTree, MarginalTree *pST );
    
    GenealogicalNetwork *pnetworkOpt;
    vector<PhylogenyTreeBasic *> &listGeneTrees;
    GenealogicalNetwork &networkInit;
    TaxaMapper &mapperTaxaIds;
    int maxNumMixNodes;
    set<int> taxaMixInit;
    double logprobBestCurr;
    bool fHeuSearch;
    bool fInitNetOnly;
    static int taxonOutgroup;
};


#endif /* AGFastCoalBuilder_hpp */
