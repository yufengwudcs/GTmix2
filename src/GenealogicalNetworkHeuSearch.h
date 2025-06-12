//
//  GenealogicalNetworkHeuSearch.h
//  
//
//  Created by Yufeng Wu on 1/14/16.
//
//

#ifndef ____GenealogicalNetworkHeuSearch__
#define ____GenealogicalNetworkHeuSearch__

#include <set>
#include <vector>
#include <map>
using namespace std;

#include "PhylogenyTreeBasic.h"

class GenealogicalNetwork;
class MarginalTree;

//***********************************************************************************
// Helper routine
void GNHeuSearchInitNet( const vector<PhylogenyTreeBasic *> &listGenetreePtrs, GenealogicalNetwork &netInit  );


//***********************************************************************************
// Basic info about gene tree

class GNHeuSearchGTreeInfo
{
public:
    GNHeuSearchGTreeInfo( PhylogenyTreeBasic &gtree );
    int GetNumOfAlleles(int taxon) const;
    int GetNumOfAllelesFor(const set<int> &taxa) const;
    int GetNumOfLeavesUnder(TreeNode *pncurr) const;
    set<int> &GetAllTaxa() { return setTaxa; }
    set<TreeNode *> &GetSetLeavesFor(int taxon);
    int GetNumEdgesBtwLeaves( TreeNode *pLeaf1, TreeNode *pLeaf2 );
    
private:
    void Init( PhylogenyTreeBasic &gtree );
    
    map<int,int> mapTaxonAlleleCounts;
    map<TreeNode *, multiset<int> > mapNodeLeafIds;
    set<int> setTaxa;
    map<int, set<TreeNode *> > mapTaxonToLeaves;
    map< pair<TreeNode *, TreeNode *>, int > mapEdgeNumBetweenPairLeaves;
};



//***********************************************************************************
// How much support a cluster receives

class GNHeuSearchClusterSupport
{
public:
    GNHeuSearchClusterSupport() : outgroup(-1), numTrees(100) {}
    void AddSupport( const set<int> &setTaxa, double val );
    void AddSupport2( const set<int> &setTaxa, double val );
    void AddDistEst(int taxon1, int taxon2, double distEst);
    int GetNumTaxa() const;
    string BuildNJTree() const;
    string ReWeightTree(const string &strTree) const;
    int GetOutgroup() const { return outgroup; }
    void SetOutgroup(int og) { outgroup = og; }
    void SetNumTrees(int nt) { numTrees = nt; }
    void Dump() const;
    
private:
    map<set<int>, double> mapClusterSupport;
    map<set<int>, double> mapClusterSupport2;
    map< pair<int,int>, double > mapTaxaPairDist;
    map< pair<int,int>, int > mapTaxaPairDistNumEst;
    int outgroup;
    int numTrees;
};


//***********************************************************************************
// Evaluate the tree

class GNHeuSearchEvaluator
{
public:
    void Evaluate( PhylogenyTreeBasic &treeGene, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport );
    void Evaluate2( PhylogenyTreeBasic &treeGene, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport );
    
private:
    void EvaluateTaxaDist( int taxon1, int taxon2, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport  );
    void EvaluateNode( TreeNode *pncurr, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport );
    //void EvaluateNode2( TreeNode *pncurr, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport, const map<int,double> &mapFreqs );
};

//***********************************************************************************
// Search network based on subtree cluster scoring
// i.e. prefer networks with marginal tree that provides more support from subnet structures

class GNHeuInfSubtreeCluster
{
public:
    GNHeuInfSubtreeCluster(GenealogicalNetwork &netInitIn, vector<PhylogenyTreeBasic *> &listGenetreePtrsIn, TaxaMapper *pMapperTaxonId, map< set<int>,double > &mapClusterSupportIn);
    ~GNHeuInfSubtreeCluster();
    void Infer();
    GenealogicalNetwork *GetOptNet() { return pnetworkOpt; }
    
private:
    double SearchWithNewMix();
    double ScoreNetwork( GenealogicalNetwork *pNetToEval );
    double ScoreMargTree(MarginalTree *ptree);
    double ScoreSplit(const set<int> &split);
    
    GenealogicalNetwork &netInit;
    vector<PhylogenyTreeBasic *> &listGenetreePtrs;
    TaxaMapper *pMapperTaxonId;
    map<set<int>,double> mapClusterSupport;
    GenealogicalNetwork *pnetworkOpt;
};




#endif /* defined(____GenealogicalNetworkHeuSearch__) */
