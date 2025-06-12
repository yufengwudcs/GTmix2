//
//  GenealogicalNetwork.h
//  
//
//  Created by Yufeng Wu on 12/21/15.
//
//

#ifndef ____GenealogicalNetwork__
#define ____GenealogicalNetwork__

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <string>
using namespace std;

#include "MarginalTree.h"
class TaxaMapper;

//***********************************************************************************
// Assume:
// each node either has two descendents (non-leaf) or none
// each node has either one ancestor (tree nodes), or two ancecstors (mix node) or none (root)
// assume taxaid starts from 0

//***********************************************************************************
// Network node

class GenealogicalNetworkBranch;

class GenealogicalNetworkNode
{
public:
    GenealogicalNetworkNode();
    ~GenealogicalNetworkNode();
    int GetID() const {return nodeId;}
    int GetTaxonId() const { return taxonId; }
    void SetTaxonId(int tid) { taxonId = tid; }
    string GetTaxonStrUser() const { return taxonStrUser; }
    void SetTaxonStrUser(const string &str) { taxonStrUser = str; }
    bool IsMixNode() const;
    bool IsLeaf() const;
    bool IsRoot() const;
    bool HasTaxon() const;
    GenealogicalNetworkBranch *GetLeftDesc() const { return pDescBranch1; }
    GenealogicalNetworkNode *GetLeftDescNode() const;
    GenealogicalNetworkBranch *GetRightDesc() const { return pDescBranch2; }
    GenealogicalNetworkNode *GetRightDescNode() const;
    GenealogicalNetworkNode *GetOtherDescNode(GenealogicalNetworkNode *pDescOne) const;
    int GetNumChildren() const;
    GenealogicalNetworkNode *GetChildSingle() const;
    GenealogicalNetworkBranch *GetAnces1() const { return pSrcBranch1; }
    GenealogicalNetworkBranch *GetAnces2() const { return pSrcBranch2; }
    void SwapAncesBrs();    // swap the left and right branches
    void SetParentBranch( GenealogicalNetworkBranch *pAnc );
    void SetChildBranch( GenealogicalNetworkBranch *pChild);
    void RemoveChild( GenealogicalNetworkNode *pNodeChild );
    void RemoveParent( GenealogicalNetworkNode *pParent );
    int GetNumParents() const;
    GenealogicalNetworkNode *GetParentSingle() const;
    GenealogicalNetworkBranch *GetParentBrSingle() const;
    void GetParents( vector<GenealogicalNetworkNode *> &listPars ) const;
    GenealogicalNetworkBranch *GetMixNodeDesc() const;
    void GetAllDescendantNodesUnder( set<GenealogicalNetworkNode *> &setDescsBelow) const;
    void GetAllLeavesUnder( set<GenealogicalNetworkNode *> &setLeavesBelow) const;
    void GetAllLeafTaxonIdsUnder(set<string> &setLeavesIds ) const;
    void GetAllLeafTaxonIdsUnder(set<int> &setLeavesIds ) const;
    bool IsLeftAncesBranch( GenealogicalNetworkBranch *pBr ) const;
    void SetDescendents( GenealogicalNetworkBranch *pDest1, GenealogicalNetworkBranch *pDest2 );
    double GetMixRatio() const { return ratioMix; }
    void SetMixRatio(double r) { ratioMix = r; }
    GenealogicalNetworkNode *AddToSib( GenealogicalNetworkNode *pSibToBe, double lenThis, double lenSib, GenealogicalNetworkBranch * &pbrParToOld, GenealogicalNetworkBranch * &pbrParToSib );
    bool IsAncestralTo(GenealogicalNetworkNode *pDescendent) const;
    void Dump() const;
    
private:
    void DumpTaxonId() const;
    
    int taxonId;
    string taxonStrUser;        // user assigned taxon string (can be empty; and if empty, not displayed)
    GenealogicalNetworkBranch *pDescBranch1;
    GenealogicalNetworkBranch *pDescBranch2;
    GenealogicalNetworkBranch *pSrcBranch1;
    GenealogicalNetworkBranch *pSrcBranch2;
    double ratioMix;
    int nodeId;
    
    static int idToUseNext;
};

//***********************************************************************************
// Edge

class GenealogicalNetworkBranch
{
public:
    GenealogicalNetworkBranch( GenealogicalNetworkNode *pSrc, GenealogicalNetworkNode *pDest );
    double GetLength() const;
    void SetLength(double len);
    GenealogicalNetworkNode *GetSrcNode() const { return pNodeSrc; }
    GenealogicalNetworkNode *GetDestNode() const { return pNodeDest; }
    void SetSrcNode(GenealogicalNetworkNode *pSrcNew) { pNodeSrc=pSrcNew; }
    void SetDestNode(GenealogicalNetworkNode *pDestNew) { pNodeDest=pDestNew; }
    bool IsAncestralTo( GenealogicalNetworkBranch *pBrDescend ) const;
    bool IsAdjacentToOutgroup(int outgroup) const;
    bool IsMixing() const;
    bool IsSrcMixing() const;
    bool IsLeafBranch() const;
    void DetachFromSrc();
    void DetachFromDest();
    void AttachToSrc(GenealogicalNetworkNode *pSrc, GenealogicalNetworkNode *pBrNodeeUnder);
    //void AttachToDest(GenealogicalNetworkNode *pDest, GenealogicalNetworkNode *pBrNodeeUnder);
    int GetUserData() const { return userData; }
    void SetUserData(int t) { userData = t; }
    void Dump() const;
    
private:
    GenealogicalNetworkNode *pNodeSrc;
    GenealogicalNetworkNode *pNodeDest;
    double lenBranch;
    int userData;
};

//***********************************************************************************
// Marg-tree info

class GenealogicalNetworkMTreeInfo
{
public:
    GenealogicalNetworkMTreeInfo();
    GenealogicalNetworkMTreeInfo( MarginalTree *pMargTreeIn, double freqIn );
    GenealogicalNetworkMTreeInfo(const GenealogicalNetworkMTreeInfo &rhs);
    ~GenealogicalNetworkMTreeInfo();
    MarginalTree *GetMargTree() const { return pMargTree; }
    double GetFreq() const { return freq; }
    void FreeTree() const;
    bool IsTopologicalEquiv(const GenealogicalNetworkMTreeInfo &other) const;
    void CombineWith( const GenealogicalNetworkMTreeInfo &other );
    void Dump() const;
    
private:
    MarginalTree *pMargTree;
    double freq;
};

//***********************************************************************************
// Marg-tree coding: how this tree is obtained from the network

class GenealogicalNetworkMTreeCode
{
public:
    GenealogicalNetworkMTreeCode( const vector<int> &ncode, const vector<GenealogicalNetworkNode *> &listMixNodes);
    GenealogicalNetworkMTreeCode(const GenealogicalNetworkMTreeCode &rhs);
    ~GenealogicalNetworkMTreeCode();
    bool operator<(const GenealogicalNetworkMTreeCode &rhs) const;
    const vector<int> &GetCode() const { return codeMTree; }
    const vector<GenealogicalNetworkNode *> &GetMixNodes() const { return listMixNodes; }
    void GetChosenBrs(set<GenealogicalNetworkBranch *> &setChosenBrs) const;
    bool IsLeftAtMixNode(int indexMN) const;
    void Dump() const;
    static GenealogicalNetworkMTreeCode & GetEmptyCode();
    
private:
    vector<int> codeMTree;
    vector<GenealogicalNetworkNode *> listMixNodes;
};


//***********************************************************************************
// Main interface for genealogical network (with admixture)

class GenealogicalNetwork
{
public:
    GenealogicalNetwork();
    ~GenealogicalNetwork();
    void InitWithTree(string &treeNW, int numTaxa);
    GenealogicalNetwork *Copy(map<GenealogicalNetworkNode *, GenealogicalNetworkNode *> *pMapOldNodeToNew = NULL) const;
    void AddNode( GenealogicalNetworkNode *pNode );
    void AddBranch( GenealogicalNetworkNode *pSrcNode, GenealogicalNetworkNode *pDestNode, double len);
    void RemoveBranch( GenealogicalNetworkNode *pSrcNode, GenealogicalNetworkNode *pDestNode, bool fMem = false );
    void RemoveNode(GenealogicalNetworkNode *pNode, bool fRmBrs = true);
    GenealogicalNetworkBranch *GetBranch( GenealogicalNetworkNode *pSrc, GenealogicalNetworkNode *pDest ) const;
    void GetAllBranches( set<GenealogicalNetworkBranch *> &setAllBranches ) const;
    void GetAllBranchesOrdered( vector<GenealogicalNetworkBranch *> &setAllBranches ) const;
    void GetAllBranchesOrdered2( vector<GenealogicalNetworkBranch *> &setAllBranches ) const;
    void GetAllBranchesBottomUp( vector<GenealogicalNetworkBranch *> &setAllBranches ) const;
    void RetriveAllMarginalTrees( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapMargTreesWithFreq ) const;
    //void RetriveAllSubMarginalTrees(int szSubsetTaxa, map<set<int>, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> > &mapSubMargTreesWithFreq ) const;
    void RetriveAllNetCode( set<GenealogicalNetworkMTreeCode > &setMTCodes ) const;
    void UpdateMargTreesForChangedBrLen(GenealogicalNetworkBranch *pBrToChange, double brLenNew, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapMargTreesWithFreq);
    //void RetriveAffectedMargTreesForChangedBr(GenealogicalNetworkBranch *pBrToChange, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapMargTreesWithFreq);
    void RetriveNgbrNetsOneEvt( vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup=-1 ) const;
    void RetriveNgbrNetsOneEvt2( vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup=-1, vector<set<GenealogicalNetworkBranch *> > *plistNgbrNetsCritcalBrs = NULL  ) const;
    void RetriveNgbrNetsOneNewMix( vector<GenealogicalNetwork *> &listNetsOneNewMix ) const;
    void RetriveNgbrNetsOneNewMixForEdge( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetwork *> &listNetsOneNewMix ) const;
    //void RetriveNgbrNetsChangeOneMix(vector<GenealogicalNetwork *> &listNgbrNetsDist1);
    int GetNumLeaves() const;
    int GetNumMixNodes() const { vector<GenealogicalNetworkNode*> listNodes; GetMixNodes(listNodes); return listNodes.size(); }
    GenealogicalNetworkNode *GetRoot() const;
    GenealogicalNetworkNode *GetLeafWithTaxon(int taxon) const;
    void GetWtsForAllMargTrees(vector<double> &listWts) const;
    void GetMixNodes( vector<GenealogicalNetworkNode *> &listMixNodes ) const;
    void RemoveTrivialCycles();
    void RemoveIsolatedNodes();
    bool CheckTrivialCycles() const;
    bool CheckSingleRoot() const;
    bool CheckNoChangeForReattach(GenealogicalNetworkBranch *pBrSrc, GenealogicalNetworkBranch *pBrDest) const;
    bool CheckOutgroup(int og) const;
    void MapBackUserLabels(TaxaMapper &mapperTaxa);
    void MapBackUserLabelsByUserLabels(TaxaMapper &mapperTaxa);
    void OutputGML ( const char *inFileName );
    void Dump() const;
    void DumpMargTrees(bool fConv = true) const;
    void DumpAdmixNodes() const;
    void DumpAdmixNodeIds() const;
    void ReadFromFileTreeMixFormat( const string &strFileTreemix, TaxaMapper *pMap );
    void CreateMapTaxonIdToUserId( map<string,string> &mapTaxonIdToUserId ) const;
    void GetAllTaxa(set<int> &setTaxa) const;
    void GetListofNodesTopdown( vector<GenealogicalNetworkNode *> &listNodes ) const;
    void GetListofNodesOrdered( vector<GenealogicalNetworkNode *> &listNodes ) const;
    void FindAllPathsBtwAllPairTaxa(map<pair<int,int>, set<set<GenealogicalNetworkBranch *> > > &mapPathsBtwPairTaxa) const;
    void GetLeaves( vector<GenealogicalNetworkNode *> &listLeaves ) const;
    bool IsOutgroup(int taxon) const;
    void AddMixBtwBranches( GenealogicalNetworkBranch *pBrMixSink, GenealogicalNetworkBranch *pBrMixSource );
    void MapNetCodedEdgeToMTBrs(const GenealogicalNetworkMTreeCode &ncode, const MarginalTree &margTree, std::map<GenealogicalNetworkBranch *, int> &mapChosenEdgeToMargBr ) const;
    void GetTopologySignature(std::map<string,int> &margTreeCnts) const;
    static void SetNNIMode(bool f );
    static bool IsNNIMode();
    
private:
    void GetMTCodesForMixNodes( const vector<GenealogicalNetworkNode *> &listMixNodes, set<GenealogicalNetworkMTreeCode > &setMTCodes ) const;
    void RetriveMargTreeForCode( const GenealogicalNetworkMTreeCode &ncode, MarginalTree &margTree ) const;
    double GetMargTreeFreqForCode( const GenealogicalNetworkMTreeCode &ncode ) const;
    void GetChangedMargTreeCodesByBrLen( GenealogicalNetworkBranch *pBranchChanged, set<GenealogicalNetworkMTreeCode> &setChangedMargTreeCodes );
    void RetriveNgbrNetsOneEvtForEdge( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup, vector<set<GenealogicalNetworkBranch *> > *pListCriticalBrs = NULL ) const;
    void RetriveNgbrNetsOneEvtForMixEdge( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetwork *> &listNgbrNetsDist1, vector<set<GenealogicalNetworkBranch *> > *pListCriticalBrs = NULL ) const;
    void RetriveNgbrNetsOneEvtChangeMixEdge( GenealogicalNetworkBranch *pMixBrToChange, vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup, vector<set<GenealogicalNetworkBranch *> > *pListCriticalBrs = NULL ) const;
    bool RecSearchForClusters( const set<GenealogicalNetworkBranch *> &setChosenBrs, GenealogicalNetworkNode *pnCurr, int &posToUseNext, stack< int > &stackParPos, map<int,int> &mapParPos, map<int, double> &mapBrLen ) const;
    // YW: this assumes taxa id starts from 0
    int GetLeafPosById(int taxonId) const { return taxonId; }
    GenealogicalNetworkBranch * GetSelAncBrFrom( GenealogicalNetworkNode *pnCurr, const set<GenealogicalNetworkBranch *> &setChosenBrs ) const;
    void FindReattachBranches( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetworkBranch *> &setBranchesToAttach, bool fUseNNI ) const;
    GenealogicalNetworkNode * ReattachBranchTo( GenealogicalNetworkBranch *pBrSrc, GenealogicalNetworkBranch *pBrDest );
    void ReattachBranchAsRootSib(GenealogicalNetworkBranch *pBrSrc);
    void CleanupAt( GenealogicalNetworkNode *pNodeChanged );
    //void RetriveSubMarginalTreesFromCompleteTrees(const set<int> &taxaSet, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapCompleteMargTreesWithFreq, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapSubMargTreesWithFreq ) const;
    bool IsBranchIn( GenealogicalNetworkBranch *pBr, const set<GenealogicalNetworkBranch *> &setBrs ) const;
    void GetAncesBranchesFrom(GenealogicalNetworkNode *pNode, set<GenealogicalNetworkBranch *> &setAncBrs) const;
    void PruneNonexistBrs(set<GenealogicalNetworkBranch *> &setAncBrs) const;
    
    set<GenealogicalNetworkNode *> setNodes;
    set<GenealogicalNetworkBranch *> setBranches;
    vector<GenealogicalNetworkNode *> listNodesAdded;
    static bool fNNIMode;
};


#endif /* defined(____GenealogicalNetwork__) */
