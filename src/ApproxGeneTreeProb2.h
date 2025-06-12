//
//  ApproxGeneTreeProb2.h
//  
//
//  Created by Yufeng Wu on 9/26/24.
//  Approximate gene tree probability by dropping the coalescent coefficient terms
//  Designed to handle larger number of trees
//  Different from the original version, this
//  approximation has polynomial-time running time

#ifndef ____ApproxGeneTreeProb2__
#define ____ApproxGeneTreeProb2__

#include "GeneSpeciesTreeProb.h"
#include "Utils4.h"
#include <cmath>
#include <string>
#include <unordered_map>

///////////////////////////////////////////////////////////////////////////////////////
void AGTPSetVerbose(bool f);

///////////////////////////////////////////////////////////////////////////////////////
// for benchmark run time

class ApproxGTPStats
{
public:
    static ApproxGTPStats &Instance();
    void DumpStats() const;
    void RecordProbComputeStart();
    void RecordProbComputeEnd();
    void RecordMaxACSize(int szACList);
    void RecordBranchLenOpt();
    
private:
    ApproxGTPStats();
    
    vector<long> listProbCalcTimeStart;
    vector<long> listProbCalcTimeEnd;
    int maxSzACList;
    int brOptNums;
};


///////////////////////////////////////////////////////////////////////////////////////
// supporting multi-threading

class ApproxGTPCache
{
public:
    static ApproxGTPCache &Instance();
    void Init();
    //map< pair<int,int>, double> &GetdbValCache(int tid);
    //map< pair<int,int>, double> &GetwbValCache(int tid);
    //map< pair<pair<int,int>,int>, double > &GetBrProbCache(int tid);
    //double GetFac(int k) const;
    double GetProdChoose2(int ub, int cb) const;
    void DumpStats() const;
    //void NotedbValCacheHit() { ++numCacheddVals; }
    //void NotewbValCacheHit() { ++numCachedwbVals; }
    
private:
    ApproxGTPCache();
    ~ApproxGTPCache();

    //vector<map< pair<int,int>, double> *> listMapdbValsThreadSafe;
    //vector<map< pair<int,int>, double> * > listMapwbValsThreadSafe;
    //vector<map< pair<pair<int,int>,int>, double > *> listCacheBranchProbThreadSafe;
    
    // additional
    //vector<double> listFactorials;
    vector<vector<double> > listProdChoose2;    // [ub][cb]: ub lineages coalesce cb times; ub >=cb+1
    
    // stats for cache performance
    //int numTotdVals;
    //int numCacheddVals;
    //int numTotwbVals;
    //int numCachedwbVals;
};

///////////////////////////////////////////////////////////////////////////////////////
struct KeyHasher
{
  std::size_t operator()(const pair<pair<int,int>,int>& k) const
  {
    using std::size_t;
    using std::hash;
    using std::string;

    return ((hash<int>()(k.first.first)
             ^ (hash<int>()(k.first.second) << 1)) >> 1)
             ^ (hash<int>()(k.second) << 1);
  }
};

///////////////////////////////////////////////////////////////////////////////////////
// pre-compute some commonly used coalescent terms

class ApproxGTPCache2
{
public:
    static ApproxGTPCache2 &Instance();
    void Init(int maxNumLins);
    void InitNumThreads(int nt);
    double GetdbValCache(int ub, int cb);
    double GetwbValCache(int ub, int cb);
    //double GetBrProbCache(int u, int v, double brLen);
    int GetMaxBrProbPreValu() const;
    //map< pair<pair<int,int>,int>, double > &GetBrProbCache(int tid);
    //double GetFac(int k) const;
    //double GetProdChoose2(int ub, int cb) const;
    //double GetCoalProbBtwCfgLins( int u, int v, double brLen );
    double GetCoalProbBtwCfgLins2( int u, int v, double brLen, int tid );
    double GetCoalProbBtwCfgLins2ConvertedLen( int u, int v, double brLen, int grLen, int tid );     // grLen: grid of length
    double GetCoalProbBtwCfgLins2ConvertedLenLRU( int u, int v, double brLen, int grLen, int tid );     // grLen: grid of length
    int ConvDistToGrid(double dist) const;
    double ConvGridToDist(int grid) const;
    void DumpStats() const;
    
private:
    ApproxGTPCache2();
    ~ApproxGTPCache2();
    //bool IsCacheBrProbInit() const { return listCacheBranchProb.size() > 0; }

    vector<vector<double> > listdbVals;
    //vector<vector<double> > listwbVals;
    vector<vector<vector<double> > > listCacheBranchProb;
    //vector<map< pair<pair<int,int>,int>, double > *> listCacheBranchProbThreadSafe;
    // for thread-safe access
    vector<map<pair<pair<int,int>,int>, double> > listCachedCoalProb;
    
    // another version of cache based LRU
    vector< lru_cache2< pair<pair<int,int>,int>, double, KeyHasher > * > listCachedCoalProbLRU;
    
    
    // additional
    //vector<double> listFactorials;
    //vector<vector<double> > listProdChoose2;    // [ub][cb]: ub lineages coalesce cb times; ub >=cb+1
};

//////////////////////////////////////////////////////////////////////////////////
// Approx gene tree prob clustser

class ApproxGeneTreeProbCluster
{
public:
    ApproxGeneTreeProbCluster(TreeNode *proot) : pClusterRoot(proot)
    {
        InitClusterSize();
    }
    ApproxGeneTreeProbCluster( const ApproxGeneTreeProbCluster &rhs) : pClusterRoot(rhs.pClusterRoot), clusSize(rhs.clusSize) {}
    
    // the size of this cluster
    //int GetSize() const;
    TreeNode *GetClusterRoot() const { return pClusterRoot; }
    int GetClusterSize() const
    {
        //set<TreeNode *> setDescendents;
        //pClusterRoot->GetAllLeavesUnder(setDescendents);
        //return setDescendents.size();
        return clusSize;
    }
    void InitClusterSize()
    {
        set<TreeNode *> setDescendents;
        pClusterRoot->GetAllLeavesUnder(setDescendents);
        clusSize = setDescendents.size();
    }
    
private:
    TreeNode *pClusterRoot;
    int clusSize;
};

//////////////////////////////////////////////////////////////////////////////////
// Compact Approx gene tree prob helper

class CompactApproxGeneTreeProbHelper
{
public:
    CompactApproxGeneTreeProbHelper( MarginalTree &treeSpeciesIn, PhylogenyTreeBasic &treeGeneIn );
    //void GetGeneAllelesForSpecies( int taxon, set<int> &geneAlleles);
    const set<int> & GetGeneAllelesForSpecies( int taxon);
    //double CalcBranchFactorAtFake(int numDupEvts, int numTotLins, int numActiveLins) const;
    //int GetTaxaAtLeafST(int snid);
    //void GetLeafLinCountsAtSTNode(int nid, vector<int> &listCounts) const;
    //const vector<int> &GetLeafLinCountsAtSTNode(int nid) const;
    const vector<int> &GetCfgBoundsAt(int nodeST) const { return listSTLinClusterSizeList[nodeST]; }
    int GetNumAncClusters(int stNode) const { return listSTLinClusterList[stNode].size(); }
    //MarginalTree &GetSpeciesTree() const { return treeSpecies; }

    
private:
    //void Init();
    void ConsClusterInfo();
    void GetGeneAllelesForSpeciesSet( const set<int> &setTaxa, set<int> &geneAlleles);
    void GetNodesForIds(const set<int> &nids, const map<int, TreeNode*> &mapIdNode, set<TreeNode *> &setNodes) const;
    
    // input ST/GT trees
    MarginalTree &treeSpecies;
    PhylogenyTreeBasic &treeGene;
    
    // for each species tree branch (bottom), the list of clusters; listed in the order of species tree branch
    vector< vector<ApproxGeneTreeProbCluster> > listSTLinClusterList;
    
    // for each species tree branch, max/min number of lineages for each cluster
    vector< vector<int> > listSTLinClusterSizeList;
    //vector< vector<int> > listSTLinClusterSizeListMin;
    //vector< vector<int> > listSTLinClusterSizeListNet;     // = ub-lb
    
    // quick find cluster based on its root
    //vector< map< TreeNode *, int > > listSTLinClusterRootMap;
    
    // quick find absorbed child root
    //vector< map< TreeNode *, int> > listAborbedTreeNodes;
    
    // quicking checking whether a node is ancestral to another node
    //map<TreeNode *, set<TreeNode *> > mapGTNodeDesc;
    
    // useful info
    //map<int, TreeNode*> mapIdNode;                // in GT, which id is for which TreeNode
    
    // various cache to speedup
    vector< set<int> > listTaxonGids;               // for each taxon, the set of gene lineages
    //map<int, set<int> > mapTaxonGids;           // for each taxon, the set of gene lineages
    //map<pair<int,TreeNode *>, int> mapCLusterRoots;  // for each node at a species tree branch, its new cluster root

};



//////////////////////////////////////////////////////////////////////////////////
// Heuristic prob calculation

class ApproxGeneTreeProbHeu2
{
public:
    ApproxGeneTreeProbHeu2(MarginalTree &treeSpeciesIn, const vector<PhylogenyTreeBasic *> &listGeneTreePtrsIn, int threadIndex);
    virtual ~ApproxGeneTreeProbHeu2();
    void CalcLogProb(vector<double> &listLogProbs);
    void UpdateSTBrLen(int nodeST, double brLenNew);
    void OnUpdateSTBrLen();
    string GetSpeciesTreeNW() const;
    //void SnapSpeciesTreeLens();
    //void SnapSpeciesTreeLenAt(int node);
    //void SetThreadsNum(int nt) { numThreads = nt; }
    
    // make public for threading, well...
    void CalcUpperProbAt(int indexGT, int nodeST);
    void ConsLowerProbAt(int indexGT, int nodeST);
    int GetNumDPCellsAt(int indexGT, int nodeST) const;
    int GetNumTaxa() const { return treeSpecies.GetNumLeaves(); }
    int GetNumGTs() const { return listGeneTreePtrs.size(); }
    
private:
    // utility code
    void Init();
    void InitProbs();
    //void InitDPTblMT();
    double CalcLogProbForTree(int indexGT);
    int GetMaxLinGTLinNumAt(int indexGT, int nodeST) const;
    int GetMinLinGTLinNumAt(int indexGT, int nodeST) const;
    int ConvDPCellIndexToLinNum(int indexGT, int nodeST, int nc) const;
    int ConvLinNumToDPCellIndex(int indexGT, int nodeST, int nl) const;
    // prob of u lineages coalescing into v along a branch of length brLen
    double CalcCoalProbFor(double brLen, int u, int v);
    double CalcCoalProbForLenGrid(int grLen, double brLen, int u, int v);
    double CalcCoalCoeff(int u, int v);
    //void UpdateSTBrLenMTImp(int nodeST);
    double GetCoalProbBtwCfgLins2ConvertedLen( int u, int v, double brLen, int grLen );     // grLen: grid of length
    
    MarginalTree &treeSpecies;
    const vector<PhylogenyTreeBasic *> &listGeneTreePtrs;
    int threadIndex;
    vector<CompactApproxGeneTreeProbHelper *> listGTPtHelpers;
    
    // for DP
    vector<vector<vector<double> > > tblDPLower;
    vector<vector<vector<double> > > tblDPUpper;
    
    // for speedup
    vector<vector<int> > listMaxLinGTLinNumAt;
    vector<vector<double> > listLogCoeffs;
    map<pair<pair<int,int>,int>, double> listCachedCoalProb;
    //const int LEN_FAC_APTP2 = 10000;  // leave 4 digits after decimal point for branch length
    //map<int, map<int, map<int, double> > > mapCacheCoalProbs;
    
    // multithreading
    //int numThreads;
};


//////////////////////////////////////////////////////////////////////////////////
// Analyze common subtrees in species trees

class STSubtreeDepot
{
public:
    static STSubtreeDepot& Instance();
    int GetSTIdFor(const string &subtreeNW);
    void DumpStats() const;
    
private:
    STSubtreeDepot();

    // next id to use
    int idNextToUse;
    
    // for each string, assign an id
    map<string, int> mapNWtoId;
};



#endif /* defined(____ApproxGeneTreeProb2__) */
