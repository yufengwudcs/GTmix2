//
//  GenealogicalNetwork.cpp
//  
//
//  Created by Yufeng Wu on 12/21/15.
//
//

#include "GenealogicalNetwork.h"
#include "Utils4.h"
#include "MarginalTree.h"
#include "PhylogenyTreeBasic.h"
#include "UnWeightedGraph.h"
#include <string>
#include <algorithm>
#include <sstream>

//***********************************************************************************
// Network node

int GenealogicalNetworkNode :: idToUseNext = 1;

GenealogicalNetworkNode :: GenealogicalNetworkNode(): taxonId(-1), pDescBranch1(NULL), pDescBranch2(NULL), pSrcBranch1(NULL), pSrcBranch2(NULL), ratioMix(0.0)
{
    nodeId = idToUseNext++;
}
GenealogicalNetworkNode :: ~GenealogicalNetworkNode()
{
}
bool GenealogicalNetworkNode :: IsMixNode() const
{
    return GetNumParents() == 2;
}
bool GenealogicalNetworkNode :: IsLeaf() const
{
    return GetNumChildren() == 0;
    // isolated nodes are not leaves
    //return GetNumChildren() == 0 && GetNumParents() > 0;
}
bool GenealogicalNetworkNode :: IsRoot() const
{
    return GetNumParents() == 0;
}
bool GenealogicalNetworkNode :: HasTaxon() const
{
    return GetTaxonId() >= 0;
}
GenealogicalNetworkNode * GenealogicalNetworkNode :: GetOtherDescNode(GenealogicalNetworkNode *pDescOne) const
{
    // get the descendent other than pDescOne
    if( GetLeftDescNode()==pDescOne )
    {
        return GetRightDescNode();
    }
    else
    {
        return GetLeftDescNode();
    }
}
int GenealogicalNetworkNode :: GetNumChildren() const
{
    int res = 0;
    if( pDescBranch1 != NULL )
    {
        ++res;
    }
    if( pDescBranch2 != NULL)
    {
        ++res;
    }
    return res;
}
GenealogicalNetworkNode * GenealogicalNetworkNode :: GetChildSingle() const
{
    // return one child
    if( pDescBranch1 != NULL)
    {
        return pDescBranch1->GetDestNode();
    }
    else if( pDescBranch2 != NULL )
    {
        return pDescBranch2->GetDestNode();
    }
    else
    {
        return NULL;
    }
}
void GenealogicalNetworkNode :: RemoveChild( GenealogicalNetworkNode *pNodeChild )
{
    pNodeChild->RemoveParent(this);
    // remove the specified child
    if( pDescBranch1 != NULL)
    {
        if( pDescBranch1->GetDestNode() == pNodeChild )
        {
            //
            //delete pDescBranch1;
            pDescBranch1 = NULL;
            return;
        }
    }
    if( pDescBranch2 != NULL )
    {
        if( pDescBranch2->GetDestNode() == pNodeChild )
        {
            //delete pDescBranch2;
            pDescBranch2 = NULL;
            return;
        }
    }
}
void GenealogicalNetworkNode :: RemoveParent( GenealogicalNetworkNode *pParent )
{
    //
    if( pSrcBranch1 != NULL)
    {
        if( pSrcBranch1->GetSrcNode() == pParent)
        {
            pSrcBranch1 = NULL;
            return;
        }
    }
    if( pSrcBranch2 != NULL)
    {
        if( pSrcBranch2->GetSrcNode() == pParent)
        {
            pSrcBranch2 = NULL;
            return;
        }
    }
}
int GenealogicalNetworkNode :: GetNumParents() const
{
    int res = 0;
    if( pSrcBranch1 != NULL )
    {
        ++res;
    }
    if( pSrcBranch2 != NULL )
    {
        ++res;
    }
    return res;
}

GenealogicalNetworkBranch * GenealogicalNetworkNode :: GetParentBrSingle() const
{
    if( pSrcBranch1 != NULL)
    {
        return pSrcBranch1;
    }
    else if( pSrcBranch2 != NULL )
    {
        return pSrcBranch2;
    }
    else
    {
        return NULL;
    }
}

GenealogicalNetworkNode * GenealogicalNetworkNode :: GetParentSingle() const
{
    if( pSrcBranch1 != NULL)
    {
        return pSrcBranch1->GetSrcNode();
    }
    else if( pSrcBranch2 != NULL )
    {
        return pSrcBranch2->GetSrcNode();
    }
    else
    {
        return NULL;
    }
}

void GenealogicalNetworkNode :: GetAllDescendantNodesUnder( set<GenealogicalNetworkNode *> &setDescsBelow) const
{
    //
    //
    if(IsLeaf())
    {
        GenealogicalNetworkNode *pthis = const_cast<GenealogicalNetworkNode *>(this);
        setDescsBelow.insert(pthis);
        return;
    }
    if( GetLeftDesc() != NULL )
    {
        GenealogicalNetworkNode *pLeft = GetLeftDescNode();
        if( pLeft != NULL )
        {
            setDescsBelow.insert(pLeft);
            pLeft->GetAllDescendantNodesUnder(setDescsBelow);
        }
    }
    if( GetRightDesc() != NULL )
    {
        GenealogicalNetworkNode *pRight = GetRightDescNode();
        if( pRight != NULL )
        {
            setDescsBelow.insert(pRight);
            pRight->GetAllDescendantNodesUnder(setDescsBelow);
        }
    }
}

void GenealogicalNetworkNode :: GetAllLeavesUnder( set<GenealogicalNetworkNode *> &setLeavesBelow) const
{
    //
    if(IsLeaf())
    {
        GenealogicalNetworkNode *pthis = const_cast<GenealogicalNetworkNode *>(this);
        setLeavesBelow.insert(pthis);
        return;
    }
    if( GetLeftDesc() != NULL )
    {
        GenealogicalNetworkNode *pLeft = GetLeftDescNode();
        if( pLeft != NULL )
        {
            pLeft->GetAllLeavesUnder(setLeavesBelow);
        }
    }
    if( GetRightDesc() != NULL )
    {
        GenealogicalNetworkNode *pRight = GetRightDescNode();
        if( pRight != NULL )
        {
            pRight->GetAllLeavesUnder(setLeavesBelow);
        }
    }
}

void GenealogicalNetworkNode :: GetAllLeafTaxonIdsUnder(set<string> &setLeavesIds ) const
{
    //
    setLeavesIds.clear();
    set<GenealogicalNetworkNode *> setLeavesBelow;
    GetAllLeavesUnder( setLeavesBelow );
    for( set<GenealogicalNetworkNode *> :: iterator it = setLeavesBelow.begin(); it != setLeavesBelow.end(); ++it )
    {
        string strLv = (*it)->GetTaxonStrUser();
        setLeavesIds.insert(strLv);
    }
}
void GenealogicalNetworkNode :: GetAllLeafTaxonIdsUnder(set<int> &setLeavesIds ) const
{
    //
    setLeavesIds.clear();
    set<GenealogicalNetworkNode *> setLeavesBelow;
    GetAllLeavesUnder( setLeavesBelow );
    for( set<GenealogicalNetworkNode *> :: iterator it = setLeavesBelow.begin(); it != setLeavesBelow.end(); ++it )
    {
        setLeavesIds.insert((*it)->GetTaxonId());
    }
}

bool GenealogicalNetworkNode :: IsLeftAncesBranch( GenealogicalNetworkBranch *pBr ) const
{
    if( pBr == pSrcBranch1 )
    {
        return true;
    }
    else
    {
        return false;
    }
}
void GenealogicalNetworkNode :: SetDescendents( GenealogicalNetworkBranch *pDest1, GenealogicalNetworkBranch *pDest2 )
{
    pDescBranch1 = pDest1;
    pDescBranch2 = pDest2;
}
GenealogicalNetworkNode * GenealogicalNetworkNode :: AddToSib( GenealogicalNetworkNode *pSibToBe, double lenThis, double lenSib, GenealogicalNetworkBranch * &pbrParToOld, GenealogicalNetworkBranch * &pbrParToSib )
{
    // create a new node and add this node as sibling
#if 0
cout << "AddToSib: lenThis: " << lenThis << ", lenSib: " << lenSib << ": pSibToBe: ";
pSibToBe->Dump();
cout << endl;
#endif
    GenealogicalNetworkNode *pNodeParNew = new GenealogicalNetworkNode();
    pbrParToOld = new GenealogicalNetworkBranch( pNodeParNew, this );
    pbrParToOld->SetLength(lenThis);
    pbrParToSib = new GenealogicalNetworkBranch( pNodeParNew, pSibToBe);
    pbrParToSib->SetLength(lenSib);
    pNodeParNew->SetChildBranch( pbrParToOld );
    pNodeParNew->SetChildBranch( pbrParToSib);
    this->SetParentBranch(pbrParToOld);
    pSibToBe->SetParentBranch(pbrParToSib);
    return pNodeParNew;
}
void GenealogicalNetworkNode :: SetParentBranch( GenealogicalNetworkBranch *pAnc )
{
    //
    pAnc->SetDestNode(this);
    if( pSrcBranch1 == NULL )
    {
        //
        pSrcBranch1 = pAnc;
    }
    else if( pSrcBranch2 == NULL)
    {
        pSrcBranch2 = pAnc;
    }
    else
    {
        YW_ASSERT_INFO(false, "SetParentBranch");
    }
}
void GenealogicalNetworkNode :: SetChildBranch( GenealogicalNetworkBranch *pChild)
{
#if 0
this->Dump();
cout << endl;
#endif
    
    pChild->SetSrcNode(this);
    
    //
    if( pDescBranch1 == NULL )
    {
        //
        pDescBranch1 = pChild;
    }
    else if( pDescBranch2 == NULL)
    {
        pDescBranch2 = pChild;
    }
    else
    {
        this->Dump();
        cout << endl;
        YW_ASSERT_INFO(false, "SetChildBranch");
    }
}
bool GenealogicalNetworkNode :: IsAncestralTo(GenealogicalNetworkNode *pDescendent) const
{
    // allow multiple hops
    if( pDescendent == NULL )
    {
        return false;
    }
    GenealogicalNetworkNode *pthis=const_cast<GenealogicalNetworkNode*>(this);
    if( pDescendent == pthis)
    {
        return true;
    }
    bool fLeft = false;
    bool fRight = false;
    if( pDescBranch1 != NULL )
    {
        fLeft = pDescBranch1->GetDestNode()->IsAncestralTo(pDescendent);
    }
    if( fLeft == true )
    {
        return true;
    }
    if(pDescBranch2 != NULL)
    {
        fRight = pDescBranch2->GetDestNode()->IsAncestralTo(pDescendent);
    }
    return fRight;
}
void GenealogicalNetworkNode :: Dump() const
{
    cout << "NodeID[" << nodeId << "]:tid[";
    DumpTaxonId();
    cout << "]:mr[" << ratioMix << "]:";
    if( pDescBranch1!=NULL)
    {
        YW_ASSERT_INFO( pDescBranch1->GetDestNode()!=NULL, "NULL1" );
        cout << "Destination node(left): " << pDescBranch1->GetDestNode()->GetID()<<"  ";
        cout << "length: " << pDescBranch1->GetLength() << " ";
    }
    if( pDescBranch2!=NULL)
    {
        YW_ASSERT_INFO( pDescBranch2->GetDestNode()!=NULL, "NULL2" );
        cout << "Destination node(right): " << pDescBranch2->GetDestNode()->GetID()<<"  ";
        cout << "length: " << pDescBranch2->GetLength() << " ";
    }
    if( pSrcBranch1!=NULL)
    {
        YW_ASSERT_INFO( pSrcBranch1->GetSrcNode()!=NULL, "NULL3" );
        cout << "Source node(left): " << pSrcBranch1->GetSrcNode()->GetID()<<"  ";
        cout << "length: " << pSrcBranch1->GetLength() << " ";
    }
    if( pSrcBranch2!=NULL)
    {
        YW_ASSERT_INFO( pSrcBranch2->GetSrcNode()!=NULL, "NULL4" );
        cout << "Source node(right): " << pSrcBranch2->GetSrcNode()->GetID()<<"  ";
        cout << "length: " << pSrcBranch2->GetLength() << " ";
    }
    cout << endl;
//cout << "Done with node\n";
}
void GenealogicalNetworkNode :: DumpTaxonId() const
{
    // if user id string is set, use it; otherwise use the taxon id
    if( taxonStrUser.size() > 0 )
    {
        cout << taxonStrUser;
    }
    else
    {
        cout << taxonId;
    }
}
GenealogicalNetworkNode * GenealogicalNetworkNode :: GetLeftDescNode() const
{
    if( GetLeftDesc() == NULL )
    {
        return NULL;
    }
    return GetLeftDesc()->GetDestNode();
}

GenealogicalNetworkNode * GenealogicalNetworkNode :: GetRightDescNode() const
{
    if( GetRightDesc() == NULL )
    {
        return NULL;
    }
    return GetRightDesc()->GetDestNode();
}

void GenealogicalNetworkNode :: GetParents( vector<GenealogicalNetworkNode *> &listPars ) const
{
    listPars.clear();
    GenealogicalNetworkBranch *pparBr1 = GetAnces1();
    GenealogicalNetworkNode *ppar1 = NULL;
    if( pparBr1 != NULL )
    {
        GenealogicalNetworkNode *pnca = pparBr1->GetSrcNode();
        listPars.push_back(pnca);
    }
    GenealogicalNetworkBranch *pparBr2 = GetAnces2();
    GenealogicalNetworkNode *ppar2 = NULL;
    if( pparBr2 != NULL )
    {
        GenealogicalNetworkNode *pnca = pparBr2->GetSrcNode();
        listPars.push_back(pnca);
    }
}

void GenealogicalNetworkNode :: SwapAncesBrs()
{
    // swap the left and right branches
    // do nothing if the right branch is NULL
    if( GetAnces2() != NULL )
    {
        GenealogicalNetworkBranch *pbrt = GetAnces2();
        pSrcBranch2 = pSrcBranch1;
        pSrcBranch1 = pbrt;
    }
}

//***********************************************************************************
// Edge

GenealogicalNetworkBranch :: GenealogicalNetworkBranch( GenealogicalNetworkNode *pSrc, GenealogicalNetworkNode *pDest ): pNodeSrc(pSrc), pNodeDest(pDest), lenBranch(0.0), userData(0)
{
}

double GenealogicalNetworkBranch :: GetLength() const
{
    return lenBranch;
}
void GenealogicalNetworkBranch :: SetLength(double len)
{
    lenBranch = len;
}

bool GenealogicalNetworkBranch :: IsMixing() const
{
    //  mixing if the dest node is mixing node
    return pNodeDest->IsMixNode();
}
bool GenealogicalNetworkBranch :: IsSrcMixing() const
{
    return pNodeSrc->IsMixNode();
}
bool GenealogicalNetworkBranch :: IsLeafBranch() const
{
    return pNodeDest->IsLeaf();
}
void GenealogicalNetworkBranch :: DetachFromSrc()
{
    GenealogicalNetworkNode *pSrc = GetSrcNode();
    GenealogicalNetworkNode *pDest = GetDestNode();
    pSrc->RemoveChild(pDest);
}
void GenealogicalNetworkBranch :: DetachFromDest()
{
    // exactly the same as detach from src???
    GenealogicalNetworkNode *pSrc = GetSrcNode();
    GenealogicalNetworkNode *pDest = GetDestNode();
    pSrc->RemoveChild(pDest);
}
void GenealogicalNetworkBranch :: AttachToSrc(GenealogicalNetworkNode *pSrc, GenealogicalNetworkNode *pBrNodeeUnder)
{
    SetSrcNode(pSrc);
    SetDestNode(pBrNodeeUnder);
    pSrc->SetChildBranch(this);
    pBrNodeeUnder->SetParentBranch(this);
}
//void GenealogicalNetworkBranch :: AttachToDest(GenealogicalNetworkNode *pDest)
//{
//    SetDestNode(pDest);
//    //SetDestNode(pBrNodeeUnder);
//    //pSrc->SetChildBranch(this);
//    pDest->SetParentBranch(this);
//}
void GenealogicalNetworkBranch :: Dump() const
{
    //
    cout << "Br[" << GetLength() << "]: from node:";
    cout << endl;
    YW_ASSERT_INFO(pNodeSrc != NULL, "src: null");
    YW_ASSERT_INFO(pNodeDest!=NULL, "dest: null");
    pNodeSrc->Dump();
    cout << " to node: ";
    pNodeDest->Dump();
    cout << "    ";
}

bool GenealogicalNetworkBranch :: IsAncestralTo( GenealogicalNetworkBranch *pBrDescend ) const
{
    //
    GenealogicalNetworkNode *pNodeAnc = this->GetDestNode();
    GenealogicalNetworkNode *pNodeDesc = pBrDescend->GetSrcNode();
    return pNodeAnc->IsAncestralTo(pNodeDesc);
}

bool GenealogicalNetworkBranch :: IsAdjacentToOutgroup(int outgroup) const
{
//cout << "IsAdjacentToOutgroup: branch: ";
//Dump();
//cout << endl;
    // is this branch next to outgroup (either descendent is or its )
    if(outgroup < 0)
    {
        return false;
    }
    //int outgroupUse = outgroup-1;
    int outgroupUse = outgroup;
    YW_ASSERT_INFO(outgroupUse>=0, "Fail");
    if( GetDestNode()->IsLeaf() && GetDestNode()->GetTaxonId()==outgroupUse )
    {
//cout << "adj1\n";
        return true;
    }
#if 0
    GenealogicalNetworkNode *pOther = GetSrcNode()->GetOtherDescNode( GetDestNode() );
    if( pOther != NULL && pOther->IsLeaf() && pOther->GetTaxonId() ==outgroupUse)
    {
cout << "adj2\n";
        return true;
    }
#endif
//cout << "no-adj\n";
    return false;
}

//***********************************************************************************
// Marg-tree info

GenealogicalNetworkMTreeInfo :: GenealogicalNetworkMTreeInfo() : pMargTree(NULL), freq(0.0)
{
}

GenealogicalNetworkMTreeInfo :: GenealogicalNetworkMTreeInfo( MarginalTree *pMargTreeIn, double freqIn ) : pMargTree(pMargTreeIn), freq(freqIn)
{
}

GenealogicalNetworkMTreeInfo :: GenealogicalNetworkMTreeInfo(const GenealogicalNetworkMTreeInfo &rhs) : pMargTree(rhs.pMargTree), freq(rhs.freq)
{
}
GenealogicalNetworkMTreeInfo :: ~GenealogicalNetworkMTreeInfo()
{
}

bool GenealogicalNetworkMTreeInfo :: IsTopologicalEquiv(const GenealogicalNetworkMTreeInfo &other) const
{
    //
    return pMargTree->IsToplogicSame(*other.pMargTree);
}
void GenealogicalNetworkMTreeInfo :: CombineWith( const GenealogicalNetworkMTreeInfo &other )
{
    // assume the two are identical
    this->freq += other.freq;
}
void GenealogicalNetworkMTreeInfo :: FreeTree() const
{
    GenealogicalNetworkMTreeInfo *pthis=const_cast<GenealogicalNetworkMTreeInfo *>(this);
    if( pMargTree != NULL )
    {
        delete pthis->pMargTree;
        pthis->pMargTree = NULL;
    }
}

void GenealogicalNetworkMTreeInfo :: Dump() const
{
    // format: [freq] newick-format-tree
    cout << "[" << GetFreq() << "] ";
    string strNWTree = GetMargTree()->GetNewick() ;
    cout << strNWTree << endl;
}

//***********************************************************************************
// Marg-tree coding: how this tree is obtained from the network


GenealogicalNetworkMTreeCode :: GenealogicalNetworkMTreeCode( const vector<int> &ncode, const vector<GenealogicalNetworkNode *> &listMixNodesIn) : codeMTree(ncode), listMixNodes(listMixNodesIn)
{
}

GenealogicalNetworkMTreeCode:: GenealogicalNetworkMTreeCode(const GenealogicalNetworkMTreeCode &rhs) : codeMTree(rhs.codeMTree), listMixNodes(rhs.listMixNodes)
{
}

GenealogicalNetworkMTreeCode :: ~GenealogicalNetworkMTreeCode()
{
}

bool GenealogicalNetworkMTreeCode :: operator<(const GenealogicalNetworkMTreeCode &rhs) const
{
    return codeMTree < rhs.codeMTree;
}

void GenealogicalNetworkMTreeCode :: GetChosenBrs(set<GenealogicalNetworkBranch *> &setChosenBrs) const
{
    // given the choice, find the the list of branches chosen
    for(int i=0; i<(int)codeMTree.size(); ++i)
    {
        bool fLeft = IsLeftAtMixNode(i);
        GenealogicalNetworkBranch *pBrChosen;
        if( fLeft)
        {
            pBrChosen = listMixNodes[i]->GetAnces1();
        }
        else
        {
            pBrChosen = listMixNodes[i]->GetAnces2();
        }
        setChosenBrs.insert(pBrChosen);
    }
}

bool GenealogicalNetworkMTreeCode :: IsLeftAtMixNode(int indexMN) const
{
    YW_ASSERT_INFO(indexMN <(int)codeMTree.size(), "Overflow");
    int bpc = codeMTree[indexMN];
    if( bpc == 0 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void GenealogicalNetworkMTreeCode :: Dump() const
{
    cout << "Code: ";
    DumpIntVec( codeMTree );
}

GenealogicalNetworkMTreeCode & GenealogicalNetworkMTreeCode :: GetEmptyCode()
{
    //
    static vector<int> codeMTreeEmpty;
    static vector<GenealogicalNetworkNode *> listMixNodesEmpty;
    static GenealogicalNetworkMTreeCode instEmpty( codeMTreeEmpty, listMixNodesEmpty );
    return instEmpty;
}

//***************************************************************************
// Main interface for genealogical network (with admixture)

bool GenealogicalNetwork :: fNNIMode = false;


GenealogicalNetwork:: GenealogicalNetwork()
{
}

GenealogicalNetwork :: ~GenealogicalNetwork()
{
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it)
    {
        delete *it;
    }
    for( set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
    {
        delete *it;
    }
}

void GenealogicalNetwork :: SetNNIMode(bool f )
{
    fNNIMode = f;
}
bool GenealogicalNetwork :: IsNNIMode()
{
    return fNNIMode;
}

void GenealogicalNetwork :: InitWithTree(string &treeNW, int numTaxa)
{
    // initialize network with a single tree
    MarginalTree margTree;
    bool fTree = ReadinMarginalTreesNewickWLenString(treeNW, numTaxa, margTree);
    if(fTree == false)
    {
        cout << "Fail to construct the tree" << endl;
    }
    // add nodes first
    vector<GenealogicalNetworkNode *> listGNNodes;
    for(int i=0; i<margTree.GetTotNodesNum(); ++i)
    {
        GenealogicalNetworkNode *pnode = new GenealogicalNetworkNode;
        if( margTree.IsLeaf( i ) == true )
        {
            pnode->SetTaxonId( margTree.GetLabel(i) );
        }
        listGNNodes.push_back(pnode);
        this->AddNode(pnode);
    }
    for(int i=0; i<margTree.GetTotNodesNum(); ++i)
    {
        int ppos = margTree.GetParent(i);
        if( ppos < 0)
        {
            continue;
        }
        // add the edge
        double len = margTree.GetEdgeLen(i);
        this->AddBranch( listGNNodes[ppos], listGNNodes[i], len );
    }
}

void GenealogicalNetwork :: AddNode( GenealogicalNetworkNode *pNode )
{
    //
    setNodes.insert(pNode);
    listNodesAdded.push_back(pNode);
}

void GenealogicalNetwork :: AddBranch( GenealogicalNetworkNode *pSrcNode, GenealogicalNetworkNode *pDestNode, double len)
{
    YW_ASSERT_INFO(pSrcNode!=NULL, "AddBranch: src null");
    YW_ASSERT_INFO(pDestNode!=NULL, "AddBranch: dest null");
    //
    GenealogicalNetworkBranch *pBr = new GenealogicalNetworkBranch( pSrcNode, pDestNode );
    pBr->SetLength(len);
    pSrcNode->SetChildBranch(pBr);
    pDestNode->SetParentBranch(pBr);
    setBranches.insert(pBr);
}

void GenealogicalNetwork :: GetAllBranches( set<GenealogicalNetworkBranch *> &setAllBranches ) const
{
    //
    setAllBranches = this->setBranches;
}

void GenealogicalNetwork :: GetAllBranchesOrdered( vector<GenealogicalNetworkBranch *> &setAllBranches ) const
{
    // first order the nodes
    vector<GenealogicalNetworkNode *> listNodes ;
    map<GenealogicalNetworkNode *, int> mapIndex;
    GetListofNodesOrdered( listNodes );
    for(unsigned int i=0; i<listNodes.size(); ++i)
    {
        mapIndex[listNodes[i]] = i;
    }
    //
    vector<pair<pair<int,int>, GenealogicalNetworkBranch *> > listCmps;
    for(auto x: setBranches)
    {
        pair<int,int> pp( mapIndex[x->GetSrcNode()], mapIndex[x->GetDestNode()] );
        listCmps.push_back( std::make_pair(pp, x) );
    }
    std::sort(listCmps.begin(), listCmps.end());
    for(unsigned int i=0; i<listCmps.size(); ++i)
    {
        setAllBranches.push_back(listCmps[i].second);
    }
}

void GenealogicalNetwork :: GetAllBranchesOrdered2( vector<GenealogicalNetworkBranch *> &setAllBranches ) const
{
    // order by node ids of the two end nodes of the edge
    setAllBranches.clear();
    vector<pair<pair<int,int>, GenealogicalNetworkBranch *> > listCmps;
    for(auto x: setBranches)
    {
        pair<int,int> pp( x->GetSrcNode()->GetID(), x->GetDestNode()->GetID() );
        listCmps.push_back( std::make_pair(pp, x) );
    }
    std::sort(listCmps.begin(), listCmps.end());
    for(unsigned int i=0; i<listCmps.size(); ++i)
    {
        setAllBranches.push_back(listCmps[i].second);
    }
}

void GenealogicalNetwork :: GetAllBranchesBottomUp( vector<GenealogicalNetworkBranch *> &setAllBranches ) const
{
    // first order the nodes
    vector<GenealogicalNetworkNode *> listNodes ;
    GetListofNodesTopdown( listNodes );
    map<GenealogicalNetworkNode *, int> mapIndex;
    for(unsigned int i=0; i<listNodes.size(); ++i)
    {
        mapIndex[listNodes[i]] = -i;
    }
    //
    vector<pair<pair<int,int>, GenealogicalNetworkBranch *> > listCmps;
    for(auto x: setBranches)
    {
        pair<int,int> pp( mapIndex[x->GetSrcNode()], mapIndex[x->GetDestNode()] );
        listCmps.push_back( std::make_pair(pp, x) );
    }
    std::sort(listCmps.begin(), listCmps.end());
    for(unsigned int i=0; i<listCmps.size(); ++i)
    {
        setAllBranches.push_back(listCmps[i].second);
    }
}

void GenealogicalNetwork :: RetriveAllMarginalTrees( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapMargTreesWithFreq ) const
{
    // retrive all maginal trees (w/ their proper coding); note: marg trees are dynamic allocated
    // and need to be freed up
    vector<GenealogicalNetworkNode *> listMixNodes ;
    GetMixNodes( listMixNodes );
    
    // all code
    set<GenealogicalNetworkMTreeCode > setMTCodes;
    GetMTCodesForMixNodes( listMixNodes, setMTCodes );
#if 0
cout << "RetriveAllMarginalTrees: number of codes: " << setMTCodes.size() << "\n";
for(set<GenealogicalNetworkMTreeCode> :: iterator it = setMTCodes.begin(); it != setMTCodes.end(); ++it)
{
it->Dump();
}
#endif
    
    // ret. marg tree for each code
    if( setMTCodes.size() > 0 )
    {
        for( set<GenealogicalNetworkMTreeCode > :: iterator it = setMTCodes.begin(); it != setMTCodes.end(); ++it )
        {
            //
            MarginalTree *pMTCur = new MarginalTree;
            RetriveMargTreeForCode( *it, *pMTCur );
            double freq = GetMargTreeFreqForCode(*it);
            GenealogicalNetworkMTreeInfo mtInfo( pMTCur, freq );
            mapMargTreesWithFreq.insert( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: value_type( *it, mtInfo ) );
        }
    }
    else
    {
        // if there are no network node, then just return the single tree
        GenealogicalNetworkMTreeCode codeEmpty( GenealogicalNetworkMTreeCode :: GetEmptyCode() );
        MarginalTree *pMTCur = new MarginalTree;
        RetriveMargTreeForCode( codeEmpty, *pMTCur );
        double freq = 1.0;
        GenealogicalNetworkMTreeInfo mtInfo( pMTCur, freq );
        mapMargTreesWithFreq.insert( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: value_type( codeEmpty, mtInfo ) );
    }
}
/*
void GenealogicalNetwork :: RetriveAllSubMarginalTrees(int szSubsetTaxa, map<set<int>, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> > &mapSubMargTreesWithFreq ) const
{
//cout << "RetriveAllSubMarginalTrees: subset size: " << szSubsetTaxa << endl;
    // enumerate all subset of taxa
    set<int> setTaxaAll;
    GetAllTaxa(setTaxaAll);
    vector<int> posvec;
    YW_ASSERT_INFO( szSubsetTaxa <= (int)setTaxaAll.size(), "Subset size: too large" );
    
    // get all complete marg trees
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapCompleteMargTreesWithFreq;
    RetriveAllMarginalTrees( mapCompleteMargTreesWithFreq );
    
    GetFirstCombo( szSubsetTaxa, (int)setTaxaAll.size(), posvec );
    while (true)
    {
        set<int> posvecss;
        PopulateSetByVec(posvecss, posvec);
        map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapSubMargTreesWithFreqStep;
        RetriveSubMarginalTreesFromCompleteTrees(posvecss, mapCompleteMargTreesWithFreq, mapSubMargTreesWithFreqStep);
        mapSubMargTreesWithFreq.insert( map<set<int>, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> > :: value_type( posvecss, mapSubMargTreesWithFreqStep ) );
#if 0
cout << "** Subtrees for subset: ";
DumpIntSet(posvecss);
for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: iterator it = mapSubMargTreesWithFreqStep.begin(); it != mapSubMargTreesWithFreqStep.end(); ++it )
{
it->second.Dump();
}
#endif
        if( GetNextCombo( szSubsetTaxa, (int)setTaxaAll.size(), posvec ) ==false)
        {
            break;
        }
    }
    // free
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: iterator it = mapCompleteMargTreesWithFreq.begin(); it != mapCompleteMargTreesWithFreq.end(); ++it )
    {
        it->second.FreeTree();
    }
}*/
/*
void GenealogicalNetwork :: RetriveSubMarginalTreesFromCompleteTrees(const set<int> &taxaSet, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapCompleteMargTreesWithFreq, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapSubMargTreesWithFreq ) const
{
    // retrieve subtrees from the given complete trees; will combine trees if there are equivalent trees
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: iterator it = mapCompleteMargTreesWithFreq.begin(); it != mapCompleteMargTreesWithFreq.end(); ++it )
    {
        // obtain the subtree
        MarginalTree *pSubtree = new MarginalTree;
        map<int,int> mapTaxaTemp;
        CreateSubtreeFromLeaves(*(it->second.GetMargTree()), taxaSet, *pSubtree, mapTaxaTemp);
        GenealogicalNetworkMTreeInfo subtreeInfo( pSubtree, it->second.GetFreq() );
        
        // is this new?
        bool fNew = true;
        for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo>  :: iterator it2= mapSubMargTreesWithFreq.begin(); it2 != mapSubMargTreesWithFreq.end(); ++it2 )
        {
            if( it2->second.IsTopologicalEquiv(subtreeInfo) == true )
            {
                it2->second.CombineWith( subtreeInfo );
                fNew = false;
                delete pSubtree;
                break;;
            }
        }
        if( fNew == true )
        {
            mapSubMargTreesWithFreq.insert( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: value_type( it->first, subtreeInfo ) );
        }
    }
}
*/
void GenealogicalNetwork :: GetWtsForAllMargTrees(vector<double> &listWts) const
{
    //
    listWts.clear();
    // retrive all maginal trees (w/ their proper coding); note: marg trees are dynamic allocated
    // and need to be freed up
    vector<GenealogicalNetworkNode *> listMixNodes ;
    GetMixNodes( listMixNodes );
    
    // all code
    set<GenealogicalNetworkMTreeCode > setMTCodes;
    GetMTCodesForMixNodes( listMixNodes, setMTCodes );
    
    // ret. marg tree for each code
    for( set<GenealogicalNetworkMTreeCode > :: iterator it = setMTCodes.begin(); it != setMTCodes.end(); ++it )
    {
        //
        double freq = GetMargTreeFreqForCode(*it);
        listWts.push_back(freq);
    }

}

void GenealogicalNetwork :: RetriveAllNetCode( set<GenealogicalNetworkMTreeCode > &setMTCodes ) const
{
    // retrive all maginal trees (w/ their proper coding); note: marg trees are dynamic allocated
    // and need to be freed up
    vector<GenealogicalNetworkNode *> listMixNodes ;
    GetMixNodes( listMixNodes );
    
    // all code
    GetMTCodesForMixNodes( listMixNodes, setMTCodes );
}

void GenealogicalNetwork :: RetriveNgbrNetsOneEvt( vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup ) const
{
#if 0
cout << "RetriveNgbrNetsOneEvt: net:";
this->Dump();
this->DumpMargTrees();
#endif
    //
    vector<GenealogicalNetworkBranch *> setAllBranches;
    //GetAllBranches( setAllBranches);
    // use ordered branches to remove randomness
    //GetAllBranchesBottomUp( setAllBranches);
    // use another version of get branches
    GetAllBranchesOrdered2( setAllBranches);
    for( vector<GenealogicalNetworkBranch *>:: iterator it = setAllBranches.begin(); it != setAllBranches.end(); ++it )
    {
        if( (*it)!=NULL && (*it)->IsAdjacentToOutgroup(outgroup) == false )
        {
#if 0
cout << "Processing edge: ";
(*it)->Dump();
cout << endl;
#endif
            RetriveNgbrNetsOneEvtForEdge(*it, listNgbrNetsDist1, outgroup);
        }
    }
//cout << "TOTAL NUMBER OF NEIGHBORING NETWORKS: " << listNgbrNetsDist1.size() << endl;
}

void GenealogicalNetwork :: RetriveNgbrNetsOneEvt2( vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup, vector<set<GenealogicalNetworkBranch *> > *plistNgbrNetsCritcalBrs ) const
{
#if 0
cout << "RetriveNgbrNetsOneEvt2: net:";
this->Dump();
this->DumpMargTrees();
#endif
    //
    vector<GenealogicalNetworkBranch *> setAllBranches;
    //GetAllBranches( setAllBranches);
    // use ordered branches to remove randomness
    //GetAllBranchesBottomUp( setAllBranches);
    // use another version of get branches
    GetAllBranchesOrdered2( setAllBranches);
    for( vector<GenealogicalNetworkBranch *>:: iterator it = setAllBranches.begin(); it != setAllBranches.end(); ++it )
    {
        if( (*it) != NULL && (*it)->IsAdjacentToOutgroup(outgroup) == false )
        {
#if 0
cout << "Processing edge: ";
(*it)->Dump();
cout << endl;
#endif
            RetriveNgbrNetsOneEvtForEdge(*it, listNgbrNetsDist1, outgroup, plistNgbrNetsCritcalBrs);
        }
    }
    
    if( plistNgbrNetsCritcalBrs != NULL )
    {
        YW_ASSERT_INFO( plistNgbrNetsCritcalBrs->size() == listNgbrNetsDist1.size(), "Size: mismatch" );
    }
    
//cout << "TOTAL NUMBER OF NEIGHBORING NETWORKS: " << listNgbrNetsDist1.size() << endl;
}

void GenealogicalNetwork :: RetriveNgbrNetsOneNewMix( vector<GenealogicalNetwork *> &listNetsOneNewMix ) const
{
    // first the destination branch
    vector<GenealogicalNetworkBranch *> setAllBranches;
    //GetAllBranches( setAllBranches);
    // use ordered branch set to remove randomness
    GetAllBranchesBottomUp(setAllBranches);
    
    for( vector<GenealogicalNetworkBranch *>:: iterator it = setAllBranches.begin(); it != setAllBranches.end(); ++it )
    {
        RetriveNgbrNetsOneNewMixForEdge(*it, listNetsOneNewMix);
    }
}

void GenealogicalNetwork :: RemoveTrivialCycles()
{
// hack
//return;
    // remove simple cycles of a network that has a single node for both sources at a network node
    bool fChanged = false;
    while(true)
    {
        bool fCont = false;
        //
        vector<GenealogicalNetworkNode *> listMixNodes;
        GetMixNodes( listMixNodes );
        for( int i=0; i<(int)listMixNodes.size(); ++i  )
        {
            //
            if( listMixNodes[i]->GetAnces1()->GetSrcNode() == listMixNodes[i]->GetAnces2()->GetSrcNode() )
            {
#if 0
cout << "--Trivial cycle detected in network: ";
this->Dump();
#endif
                fChanged = true;
                // this is a node to remove since it has duplicate parents
                GenealogicalNetworkNode *pAncOfTwo = listMixNodes[i]->GetAnces1()->GetSrcNode();
                YW_ASSERT_INFO( pAncOfTwo != NULL, "Cannot be null1");
                YW_ASSERT_INFO(pAncOfTwo->IsMixNode() == false, "Fatal error: this node cannot be mixing node");
                GenealogicalNetworkNode *pAncOfTwoPar = pAncOfTwo->GetParentSingle();
                
#if 0
                if( pAncOfTwoPar == NULL )
                {
                    cout << "In RemoveTrivialCycles: current network is: \n";
                    Dump();
                    OutputGML("tmp.debug.gml");
                    cout << "Processing mixing node: ";
                    listMixNodes[i]->Dump();
                    cout << endl;
                }
#endif
                
                
//#if 0
                //YW_ASSERT_INFO(pAncOfTwoPar != NULL, "Cannot be null");
                // pAncofTwo must be root; then take this into consideration
                GenealogicalNetworkNode *pDescOfTwo = listMixNodes[i]->GetChildSingle();
                double len1 = GetBranch(pAncOfTwo, listMixNodes[i])->GetLength();
                RemoveBranch( pAncOfTwo, listMixNodes[i], true );
                double len2 = GetBranch( pAncOfTwo, listMixNodes[i] )->GetLength();
                RemoveBranch( pAncOfTwo, listMixNodes[i], true );
                double len4 = GetBranch(listMixNodes[i], pDescOfTwo)->GetLength();
                double lenNew = 0.5*(len1+len2)+len4;
                pDescOfTwo->RemoveParent( listMixNodes[i] );
                RemoveBranch( listMixNodes[i], pDescOfTwo, true );
                RemoveNode( listMixNodes[i], false );
                //GenealogicalNetworkNode *pNodeParNew = pAncOfTwo;
                if( pAncOfTwoPar != NULL )
                {
                    pAncOfTwoPar->RemoveChild(pAncOfTwo);
                    double len3 = GetBranch(pAncOfTwoPar, pAncOfTwo)->GetLength();
                    lenNew += len3;
                    RemoveBranch( pAncOfTwoPar, pAncOfTwo, true );
                    //pNodeParNew = pAncOfTwoPar;
                    // finally add the remaining link
                    AddBranch( pAncOfTwoPar, pDescOfTwo, lenNew );
                }
                // 10/8/24: remove after using it.
                RemoveNode( pAncOfTwo, false );

//#endif
                
#if 0
                if( pAncOfTwoPar != NULL )
                {
                    GenealogicalNetworkNode *pDescOfTwo = listMixNodes[i]->GetChildSingle();
                    double len1 = GetBranch(pAncOfTwo, listMixNodes[i])->GetLength();
                    YW_ASSERT_INFO(pAncOfTwoPar != NULL, "Cannot be null");
                    pAncOfTwoPar->RemoveChild(pAncOfTwo);
                    RemoveBranch( pAncOfTwo, listMixNodes[i], true );
                    double len2 = GetBranch( pAncOfTwo, listMixNodes[i] )->GetLength();
                    RemoveBranch( pAncOfTwo, listMixNodes[i], true );
                    double len3 = GetBranch(pAncOfTwoPar, pAncOfTwo)->GetLength();
                    RemoveBranch( pAncOfTwoPar, pAncOfTwo, true );
                    double len4 = GetBranch(listMixNodes[i], pDescOfTwo)->GetLength();
                    double lenNew = 0.5*(len1+len2)+len3+len4;
                    pDescOfTwo->RemoveParent( listMixNodes[i] );
                    RemoveBranch( listMixNodes[i], pDescOfTwo, true );
                    RemoveNode( listMixNodes[i], false );
                    RemoveNode( pAncOfTwo, false );
                    AddBranch( pAncOfTwoPar, pDescOfTwo, lenNew );
                }
//#if 0
                else
                {
                    // pAncofTwo must be root; then take this into consideration
                    // here we simply remove the root as well
#if 0
cout << "before start, pAncOfTwo: ";
pAncOfTwo->Dump();
cout << endl;
#endif
                    pAncOfTwo->RemoveChild(listMixNodes[i]);
                    pAncOfTwo->RemoveChild(listMixNodes[i]);
                    GenealogicalNetworkNode *pDescOfTwo = listMixNodes[i]->GetChildSingle();
                    listMixNodes[i]->RemoveChild(pDescOfTwo);
                    listMixNodes[i]->RemoveParent(pAncOfTwo);
                    double len1 = GetBranch(pAncOfTwo, listMixNodes[i])->GetLength();
//cout << "len1: " << len1 << endl;
                    RemoveBranch( pAncOfTwo, listMixNodes[i], true );
#if 0
cout << "(1), pAncOfTwo: ";
pAncOfTwo->Dump();
cout << endl;
cout << "remove branch1\n";
#endif
                    double len2 = GetBranch( pAncOfTwo, listMixNodes[i] )->GetLength();
//cout << "len2: " << len2 << endl;
                    RemoveBranch( pAncOfTwo, listMixNodes[i], true );
                    //pAncOfTwo->RemoveChild(listMixNodes[i]);
                    
#if 0
cout << "(2) pAncOfTwo: ";
pAncOfTwo->Dump();
cout << endl;
cout << "remove branch2\n";
#endif
                    double len4 = GetBranch(listMixNodes[i], pDescOfTwo)->GetLength();
//cout << "len4:" << len4 << endl;
                    double lenNew = 0.5*(len1+len2)+len4;
                    pDescOfTwo->RemoveParent( listMixNodes[i] );
//cout << "remove parent for descoftwo\n";
                    RemoveBranch( listMixNodes[i], pDescOfTwo, true );
//cout << "remove branch between mix and child\n";
                    RemoveNode( listMixNodes[i], false );
//cout << "remove mix node\n";
                    YW_ASSERT_INFO(pAncOfTwo!=NULL, "pAncofTwo: NULL");
                    //YW_ASSERT_INFO(pDescOfTwo!=NULL, "pDescOfTwo: NULL");
                    RemoveNode( pAncOfTwo, false );
#if 0
cout << "pAncOfTwo: ";
pAncOfTwo->Dump();
cout << endl;
cout << "pDescOfTwo:";
pDescOfTwo->Dump();
cout << endl;
#endif
//cout << "Before adding new branch: network \n";
//OutputGML("tmp.endofcycleremoval.debug2.gml");
                    //AddBranch( pAncOfTwo, pDescOfTwo, lenNew );
cout << "After cycle removing: network \n";
OutputGML("tmp.endofcycleremoval.debug.gml");
                }
//#endif
#endif
                
                fCont = true;
                break;
            }
        }
        
        if(fCont == false)
        {
            break;
        }
    }
    if( fChanged )
    {
        //cout << "After trivial cycle removal: ";
        //this->Dump();
    }
//cout << "After cycle removing: network \n";
//OutputGML("tmp.endofcycleremoval.debug.gml");
}

void GenealogicalNetwork :: RemoveIsolatedNodes()
{
    // remove node that is isolated
    for (auto it = setNodes.begin(); it != setNodes.end(); )
    {
        if((*it)->IsLeaf() && (*it)->IsRoot())
        {
#if 0
cout << "WARNING: this node is isolated: ";
(*it)->Dump();
cout << endl;
cout << "network is: ";
Dump();
#endif
            setNodes.erase(it++);
        }
        else
        {
            ++it;
        }
    }
}

bool GenealogicalNetwork :: CheckTrivialCycles() const
{
    //
    vector<GenealogicalNetworkNode *> listMixNodes;
    GetMixNodes( listMixNodes );
    for( int i=0; i<(int)listMixNodes.size(); ++i  )
    {
        //
        if( listMixNodes[i]->GetAnces1()->GetSrcNode() == listMixNodes[i]->GetAnces2()->GetSrcNode() )
        {
            return true;
        }
    }
    return false;
}

bool GenealogicalNetwork :: CheckSingleRoot() const
{
    //
    int numRoots = 0;
    for(auto x: setNodes)
    {
        if( x->IsRoot() )
        {
            ++numRoots;
        }
    }
    return numRoots == 1;
}

bool GenealogicalNetwork :: CheckNoChangeForReattach(GenealogicalNetworkBranch *pBrSrc, GenealogicalNetworkBranch *pBrDest) const
{
    // return true if removing pBrDest causes a trivial cycle
    if( pBrSrc == NULL )
    {
        return false;
    }
    GenealogicalNetworkNode *pSrcSrc = pBrSrc->GetSrcNode();
    GenealogicalNetworkNode *pSrcDest = pBrSrc->GetDestNode();
    if( pSrcSrc->IsMixNode() )
    {
        return false;
    }
    GenealogicalNetworkNode *pSrcSrcPar = pSrcSrc->GetParentSingle();
    if( pSrcSrcPar == NULL )
    {
        return false;
    }
    GenealogicalNetworkNode *pSrcSrcSib = pSrcSrc->GetOtherDescNode(pSrcDest);
    if( pSrcSrcSib->IsMixNode() == false)
    {
        return false;
    }
    vector<GenealogicalNetworkNode *> listPars;
    pSrcSrcSib->GetParents( listPars );
    for(auto x: listPars)
    {
        if( x == pSrcSrcPar )
        {
            return true;
        }
    }
    return false;
}

bool GenealogicalNetwork :: CheckOutgroup(int og) const
{
    // check if og is valid in this network
    if( og < 0 )
    {
        return true;
    }
    // make sure there is a single child of root that is the og
    GenealogicalNetworkNode *pr = GetRoot();
    YW_ASSERT_INFO(pr != NULL, "Wrong: networ has no root");
    GenealogicalNetworkNode *prleft = pr->GetLeftDescNode();
    GenealogicalNetworkNode *prright = pr->GetRightDescNode();
    YW_ASSERT_INFO( prleft != NULL && prright != NULL, "Network root cannot only have a single child" );
    if( prleft->IsLeaf() && prleft->GetTaxonId() == og )
    {
        return true;
    }
    if( prright->IsLeaf() && prright->GetTaxonId() == og )
    {
        return true;
    }
    return false;
}

void GenealogicalNetwork :: Dump() const
{
    //
    cout << "[#nodes; " << setNodes.size() << "], [#edges:" << setBranches.size() << "]\n";
    cout << "--list of nodes: ";
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it)
    {
        (*it)->Dump();
        cout << endl;
    }
    cout << "--list of edges: ";
    for(set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
    {
        (*it)->Dump();
        cout << endl;
    }
}

void GenealogicalNetwork :: GetMixNodes( vector<GenealogicalNetworkNode *> &listMixNodes ) const
{
    listMixNodes.clear();
#if 0
    for( vector<GenealogicalNetworkNode *> :: const_iterator it = listNodesAdded.begin(); it != listNodesAdded.end(); ++it)
    {
        if( (*it)->IsMixNode() && setNodes.find(*it) != setNodes.end() )
        {
            listMixNodes.push_back(*it);
        }
    }
#endif
//#if 0
    //
    vector<GenealogicalNetworkNode *> listNodesOrdered;
    GetListofNodesOrdered( listNodesOrdered );
    for( vector<GenealogicalNetworkNode *> :: iterator it = listNodesOrdered.begin(); it != listNodesOrdered.end(); ++it)
    {
        if( (*it)->IsMixNode() )
        {
            listMixNodes.push_back(*it);
        }
    }
//#endif
}

void GenealogicalNetwork :: GetLeaves( vector<GenealogicalNetworkNode *> &listLeaves ) const
{
    //
    listLeaves.clear();
    //vector<GenealogicalNetworkNode *> listNodesOrdered;
    //GetListofNodesOrdered( listNodesOrdered );
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it)
    {
        if( (*it)->IsLeaf() )
        {
            listLeaves.push_back(*it);
        }
    }
}

void GenealogicalNetwork :: GetMTCodesForMixNodes( const vector<GenealogicalNetworkNode *> &listMixNodes, set<GenealogicalNetworkMTreeCode > &setMTCodes ) const
{
#if 0
cout << "GetMTCodesForMixNodes: list of mix nodes: \n";
for(int i=0; i<(int)listMixNodes.size(); ++i)
{
listMixNodes[i]->Dump();
}
#endif
    
    //
    int numMN = listMixNodes.size();
    vector<int> listChoices;
    GetFirstMutliChoice( numMN, 2, listChoices );
    while(true)
    {
#if 0
cout << "listChoices: ";
DumpIntVec(listChoices);
#endif
        GenealogicalNetworkMTreeCode ncode( listChoices, listMixNodes );
        setMTCodes.insert(ncode);
#if 0
cout << "Code constructed: ";
ncode.Dump();
#endif
        
        if( GetNextMutliChoice( numMN, 2, listChoices ) == false )
        {
            break;
        }
    }
}

void GenealogicalNetwork :: RetriveMargTreeForCode( const GenealogicalNetworkMTreeCode &ncode, MarginalTree &margTree ) const
{
#if 0
cout << "RetriveMargTreeForCode: nocdoe: ";
ncode.Dump();
#endif
    //
    set<GenealogicalNetworkBranch *> setChosenBrs;
    ncode.GetChosenBrs(setChosenBrs);
    stack<int> stackParPos;
    map<int,int> mapParPos;
    map<int, double> mapBrLen;
    // init the
    int numTaxa = GetNumLeaves();
    int posToUseNext = numTaxa;
    YW_ASSERT_INFO(GetRoot() != NULL, "Root: null");
    RecSearchForClusters( setChosenBrs, GetRoot(), posToUseNext, stackParPos, mapParPos, mapBrLen );
#if 0
cout << "numTaxa: " << numTaxa << ", posToUseNext: " << posToUseNext << ", size of stack: " << stackParPos.size() << ", mapParPos: ";
for(map<int,int> :: iterator it = mapParPos.begin(); it != mapParPos.end(); ++it)
{
cout << "[" << it->first << "," << it->second << "]  ";
}
cout << "mapBrLen: ";
for(map<int,double> :: iterator it = mapBrLen.begin(); it != mapBrLen.end(); ++it)
{
cout << "[" << it->first << "," << it->second << "]  ";
}
#endif
    // create marg tree from the map
    vector<int> vecLabels;
    vector<int> vecpp;
    vector<double> vecBrLen;
    for( int i=0; i<2*numTaxa-1; ++i )
    {
        if( i<numTaxa)
        {
            vecLabels.push_back(i);
        }
        else
        {
            vecLabels.push_back(-1);
        }
        
        int ppos = -1;
        if( mapParPos.find(i) != mapParPos.end() )
        {
            ppos = mapParPos[i];
        }
        if( i >=numTaxa && i <2*numTaxa-2)
        {
//GenealogicalNetwork *pthis=const_cast<GenealogicalNetwork*>(this);
//            pthis->OutputGML("debug.gml");

if( ppos < 0 )
{
cout << "ppos: " << ppos << ", i=" << i << endl;
this->Dump();
GenealogicalNetwork *pthis=const_cast<GenealogicalNetwork*>(this);
pthis->OutputGML("debug.gml");
cout << "Chosen edges: \n";
for(auto xx : setChosenBrs)
{
xx->Dump();
cout << endl;
}
cout << "numTaxa: " << numTaxa << ", posToUseNext: " << posToUseNext << ", size of stack: " << stackParPos.size() << ", mapParPos: ";
for(map<int,int> :: iterator it = mapParPos.begin(); it != mapParPos.end(); ++it)
{
cout << "[" << it->first << "," << it->second << "]  ";
}
}
            YW_ASSERT_INFO(ppos >=0, "Wrong2" );
        }
        vecpp.push_back(ppos);
        double brLen = 0.0;
        if( mapBrLen.find(i) != mapBrLen.end() )
        {
            brLen = mapBrLen[i];
        }
        vecBrLen.push_back(brLen);
    }
#if 0
cout << "vecLabels: ";
DumpIntVec(vecLabels);
cout << "vecpp: ";
DumpIntVec(vecpp);
cout << "vecBrLen: ";
DumpDoubleVec(vecBrLen);
#endif
    margTree.SetNumLeaves(numTaxa);
    margTree.SetLabelList(vecLabels);
    margTree.SetParList(vecpp);
    margTree.SetBranchLenList( vecBrLen );
    margTree.BuildDescendantInfo();
    margTree.SortByLeafId();
    
    // TBD: YW debugging
    // build a standard form
    string strNW = margTree.GetNewickSorted(true);
    ReadinMarginalTreesNewickWLenString(strNW, numTaxa, margTree);
}

double GenealogicalNetwork :: GetMargTreeFreqForCode( const GenealogicalNetworkMTreeCode &ncode ) const
{
    //
    double res = 1.0;
    vector<int> codeChoices = ncode.GetCode();
    vector<GenealogicalNetworkNode *> listMixNodes = ncode.GetMixNodes();
    for(int i=0; i<(int)codeChoices.size(); ++i)
    {
        int bitstep = codeChoices[i];
        double freqNode = listMixNodes[i]->GetMixRatio();
        double resstep = freqNode;
        if( bitstep == 1 )
        {
            resstep = 1.0-freqNode;
        }
//cout << "i = " << i << ": bitstep" << bitstep << ", freqNode: " << freqNode << ", resstep: " << resstep << endl;
        res *= resstep;
    }
//cout << "---Freq for code: " << res << endl;
    return res;
}

bool GenealogicalNetwork :: IsBranchIn( GenealogicalNetworkBranch *pBr, const set<GenealogicalNetworkBranch *> &setBrs ) const
{
    //
    bool res = false;
    for(set<GenealogicalNetworkBranch *> :: const_iterator it = setBrs.begin(); it != setBrs.end(); ++it)
    {
        const GenealogicalNetworkBranch *pBrStep = *it;
        if( pBr->GetSrcNode()==pBrStep->GetSrcNode()  && pBr->GetDestNode()==pBrStep->GetDestNode() )
        {
            res = true;
            break;
        }
    }
    return res;
}

bool GenealogicalNetwork :: RecSearchForClusters( const set<GenealogicalNetworkBranch *> &setChosenBrs, GenealogicalNetworkNode *pnCurr, int &posToUseNext, stack< int > &stackParPos, map<int,int> &mapParPos, map<int, double> &mapBrLen ) const
{
    if( pnCurr == NULL)
    {
        // nothing to pursue
        return false;
    }
#if 0
cout << "RecSearchForClusters: posToUseNext: " << posToUseNext << ", pnCurr: ";
pnCurr->Dump();
cout << "Number of chosen edges: " << setChosenBrs.size() << endl;
#endif
    // for given choice at network nodes, find the list of parent/child position info (as
    // required by MarginalTree class) and also branch length of each cluster
    // pnCurr: where in the network to look for clusters
    double brLen = 0.0;
    GenealogicalNetworkBranch *pAncBr = GetSelAncBrFrom( pnCurr, setChosenBrs );
    if( pAncBr != NULL )
    {
        brLen = pAncBr->GetLength();
    }
    // top-down from pnCurr node by following the branches as given (when there is choice)
    // remember the parent position and branch length traversed
    // ASSUME: leaf already recorded (taxon: id from 1 to n)
    // return true if there are clades found below; false otherwise (no calde below)
    if( pnCurr->IsLeaf() == true )
    {
        int tid = pnCurr->GetTaxonId();
        // position is tid-1 for leaf by default
        int tpos = GetLeafPosById( tid );
//cout << "Leaf: tid = " << tid << ", tpos: " << tpos << endl;
        //mapParPos.insert( map<int,int> :: value_type(tpos, posToUseNext) );
        stackParPos.push( tpos );
        mapBrLen.insert( map<int,double> :: value_type( tpos, brLen )  );
        return true;
    }
    
    // now deal with two descendents; recursively search
    if( pnCurr->IsMixNode() == false )
    {
//cout << "Is not a mixing node: \n";
        // make sure this is selected
        GenealogicalNetworkNode *pnChildLeft = pnCurr->GetLeftDescNode();
        GenealogicalNetworkBranch *pBrChildLeft = GetBranch( pnCurr, pnChildLeft );
        //if( pBrChildLeft->IsMixing() && setChosenBrs.find(pBrChildLeft) == setChosenBrs.end() )
        if( pBrChildLeft->IsMixing() && IsBranchIn( pBrChildLeft, setChosenBrs ) ==false )
        {
            pnChildLeft = NULL;
        }
        GenealogicalNetworkNode *pnChildRight = pnCurr->GetRightDescNode();
        GenealogicalNetworkBranch *pBrChildRight = GetBranch( pnCurr, pnChildRight );
        if( pBrChildRight->IsMixing() && IsBranchIn(pBrChildRight, setChosenBrs) == false )
        //if( pBrChildRight->IsMixing() && setChosenBrs.find(pBrChildRight) == setChosenBrs.end() )
        {
            pnChildRight = NULL;
        }
        
        bool fLeft = RecSearchForClusters( setChosenBrs, pnChildLeft, posToUseNext, stackParPos, mapParPos, mapBrLen );
        bool fRight = RecSearchForClusters( setChosenBrs, pnChildRight, posToUseNext, stackParPos, mapParPos, mapBrLen );
        //GenealogicalNetworkBranch *pAncBr = GetSelAncBrFrom( pnCurr, setChosenBrs );

        // if both true (find a cluster), create a new intrnal node
        if( fLeft && fRight)
        {
            int posIntNode = posToUseNext++;
            int posChildLeft = stackParPos.top();
            stackParPos.pop();
            int posChildRight = stackParPos.top();
            stackParPos.pop();
            mapParPos.insert( map<int,int> :: value_type(posChildLeft, posIntNode) );
            mapParPos.insert( map<int,int> :: value_type(posChildRight, posIntNode) );
            stackParPos.push( posIntNode );
            mapBrLen.insert( map<int,double> :: value_type( posIntNode, brLen )  );
        }
        else if(fLeft || fRight)
        {
            int posDescSingle = stackParPos.top();
            YW_ASSERT_INFO(mapBrLen.find(posDescSingle) != mapBrLen.end(), "Fail to find");
            mapBrLen[posDescSingle] += brLen;
        }
        
        return fLeft || fRight;
    }
    else
    {
//cout << "Is mixing node.\n";
        // here assume mixing node has only one child
        // before continuing: need to check to ensure this branch is actually chosen
        GenealogicalNetworkNode *pnChildDesc = pnCurr->GetChildSingle();
        GenealogicalNetworkBranch *pBrChildDesc = GetBranch( pnCurr, pnChildDesc );
        YW_ASSERT_INFO(pBrChildDesc != NULL, "Branch is not found");
        if( pBrChildDesc->IsMixing() && setChosenBrs.find(pBrChildDesc) == setChosenBrs.end() )
        {
            pnChildDesc = NULL;
        }
        
        bool fDesc = RecSearchForClusters( setChosenBrs, pnChildDesc, posToUseNext, stackParPos, mapParPos, mapBrLen );
        if( fDesc )
        {
            int posDescSingle = stackParPos.top();
            YW_ASSERT_INFO(mapBrLen.find(posDescSingle) != mapBrLen.end(), "Fail to find");
            mapBrLen[posDescSingle] += brLen;
        }
        return fDesc;
    }
}

GenealogicalNetworkBranch * GenealogicalNetwork :: GetSelAncBrFrom( GenealogicalNetworkNode *pnCurr, const set<GenealogicalNetworkBranch *> &setChosenBrs ) const
{
    //
    if(pnCurr->IsMixNode() == false )
    {
        GenealogicalNetworkBranch *pBrAncSel = pnCurr->GetParentBrSingle();
        return pBrAncSel;
    }
    else
    {
        //
        GenealogicalNetworkBranch* pBrAnc1 = pnCurr->GetAnces1();
        GenealogicalNetworkBranch *pBrAnc2 = pnCurr->GetAnces2();
        if( setChosenBrs.find(pBrAnc1) != setChosenBrs.end() )
        {
            YW_ASSERT_INFO( setChosenBrs.find(pBrAnc2) == setChosenBrs.end(), "WRONG: both branches are chosen!"  );
            return pBrAnc1;
        }
        else if( setChosenBrs.find(pBrAnc2) != setChosenBrs.end() )
        {
            //
            return pBrAnc2;
        }
        else
        {
            YW_ASSERT_INFO(false, "Wrong: no choices can be made");
            return NULL;
        }
    }
}

void GenealogicalNetwork :: GetChangedMargTreeCodesByBrLen( GenealogicalNetworkBranch *pBranchChanged, set<GenealogicalNetworkMTreeCode> &setChangedMargTreeCodes )
{
    // retrive all maginal trees (w/ their proper coding)
    vector<GenealogicalNetworkNode *> listMixNodes ;
    GetMixNodes( listMixNodes );
    GetMTCodesForMixNodes( listMixNodes, setChangedMargTreeCodes );
    
    // if this branch is not mixing, then need to get all the marg trees
    if( pBranchChanged->IsMixing() == true)
    {
        //
        GenealogicalNetworkNode *pnMix = pBranchChanged->GetDestNode();
        int posNode = GetItemIndexInVecGen( listMixNodes, pnMix );
        YW_ASSERT_INFO(posNode >= 0, "Wrong");
        bool fLeft = pnMix->IsLeftAncesBranch(pBranchChanged);
        
        set<GenealogicalNetworkMTreeCode> setCodesNew;
        for( set<GenealogicalNetworkMTreeCode> :: iterator it = setChangedMargTreeCodes.begin(); it != setChangedMargTreeCodes.end(); ++it)
        {
            bool fNodeLeft = it->IsLeftAtMixNode( posNode );
            if( (fLeft^fNodeLeft) == false )
            {
                setCodesNew.insert(*it);
            }
        }
        setChangedMargTreeCodes = setCodesNew;
    }
}


int GenealogicalNetwork :: GetNumLeaves() const
{
    //
    vector<GenealogicalNetworkNode *> listLeaves;
    GetLeaves(listLeaves);
    return listLeaves.size();
}

GenealogicalNetworkNode * GenealogicalNetwork :: GetRoot() const
{
    //
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it)
    {
        if( (*it)->GetNumParents() == 0 )
        {
            return *it;
        }
    }
    YW_ASSERT_INFO(false, "Fail to find the root");
    return NULL;
}

GenealogicalNetworkNode * GenealogicalNetwork :: GetLeafWithTaxon(int taxon) const
{
    //
    vector<GenealogicalNetworkNode *> listLeaves;
    GetLeaves(listLeaves);
    for(int i=0; i<(int)listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == taxon )
        {
            return listLeaves[i];
        }
    }
    return NULL;
}


void GenealogicalNetwork :: UpdateMargTreesForChangedBrLen(GenealogicalNetworkBranch *pBrToChange, double brLenNew, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapMargTreesWithFreq)
{
    // create the list of marginal trees after updating one branch length
    pBrToChange->SetLength(brLenNew);
    set<GenealogicalNetworkMTreeCode> setChangedMTCodes;
    GetChangedMargTreeCodesByBrLen(pBrToChange, setChangedMTCodes);
    for( set<GenealogicalNetworkMTreeCode>::iterator it= setChangedMTCodes.begin(); it != setChangedMTCodes.end(); ++it )
    {
        //
        YW_ASSERT_INFO( mapMargTreesWithFreq.find(*it)!=mapMargTreesWithFreq.end(), "Fail to find" );
        RetriveMargTreeForCode( *it, *mapMargTreesWithFreq[*it].GetMargTree() );
    }
}
/*
void GenealogicalNetwork :: RetriveAffectedMargTreesForChangedBr(GenealogicalNetworkBranch *pBrToChange, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> &mapMargTreesWithFreq)
{
    //
    set<GenealogicalNetworkMTreeCode> setChangedMTCodes;
    GetChangedMargTreeCodesByBrLen(pBrToChange, setChangedMTCodes);
    for( set<GenealogicalNetworkMTreeCode>::iterator it= setChangedMTCodes.begin(); it != setChangedMTCodes.end(); ++it )
    {
        //
        MarginalTree *pMTCur = new MarginalTree;
        RetriveMargTreeForCode( *it, *pMTCur );
        double freq = GetMargTreeFreqForCode(*it);
        GenealogicalNetworkMTreeInfo mtInfo( pMTCur, freq );
        mapMargTreesWithFreq.insert( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: value_type( *it, mtInfo ) );
    }

}*/

GenealogicalNetworkBranch * GenealogicalNetwork :: GetBranch( GenealogicalNetworkNode *pSrc, GenealogicalNetworkNode *pDest ) const
{
    //
    GenealogicalNetworkBranch *pres = NULL;
    for( set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
    {
        if( (*it)->GetSrcNode() == pSrc && (*it)->GetDestNode() == pDest )
        {
            pres = *it;
            break;
        }
    }
    return pres;
}

void GenealogicalNetwork :: RetriveNgbrNetsOneEvtForEdge( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup, vector<set<GenealogicalNetworkBranch *> > *pListCriticalBrs ) const
{
//cout << "Entering RetriveNgbrNetsOneEvtForEdge\n";
    YW_ASSERT_INFO(pBranchToChange != NULL, "CANNOT BE NULL");
    // critial branches are those ancestral to affected edges
    if(pBranchToChange->IsSrcMixing() )
    {
//cout << "Mixing edge...\n";
        RetriveNgbrNetsOneEvtForMixEdge(pBranchToChange, listNgbrNetsDist1, pListCriticalBrs);
//cout << "    done mixing edge..\n";
        return;
    }
    
    if( pBranchToChange->IsMixing() )
    {
        //cout << "Now changing this mixing branch:\n";
        RetriveNgbrNetsOneEvtChangeMixEdge(pBranchToChange, listNgbrNetsDist1, outgroup, pListCriticalBrs);
    }
    
//cout << "Non-mixing edge...\n";
    
    // change where the branch attaches to and return all the possible networks
    vector<GenealogicalNetworkBranch *> setBrsToReattach;
    FindReattachBranches( pBranchToChange, setBrsToReattach, IsNNIMode() );
#if 0
cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvRetriveNgbrNetsOneEvtForEdge: pBranchToCHange: ";
pBranchToChange->Dump();
//cout << "  Branches to attach: ";
//for(vector<GenealogicalNetworkBranch *>::iterator it= setBrsToReattach.begin(); it != setBrsToReattach.end(); ++it)
//{
//(*it)->Dump();
//cout << endl;
//}
cout << "^^^^Number of attaching edges: " << setBrsToReattach.size() << endl;
//cout << "^^^ current net: ";
//this->Dump();
#endif
    
    // YW: if no where to re-attach, stop. TBD June 3, 2024
    //if( setBrsToReattach.size() == 0 )
    //{
    //    return;
    //}
    
    for( vector<GenealogicalNetworkBranch *> :: iterator it = setBrsToReattach.begin(); it != setBrsToReattach.end(); ++it )
    {
#if 0
cout << "Now reattach to this branch: ";
(*it)->Dump();
#endif
        
        // YW: April 12, 2024, fix an error in building network
        if( CheckNoChangeForReattach( pBranchToChange, *it ) )
        {
            continue;
        }
        
        map<GenealogicalNetworkNode *, GenealogicalNetworkNode*> mapOldNodeToNew;
        GenealogicalNetwork *pNetCopy = Copy( &mapOldNodeToNew );
#if 0
cout << "Network copy: ";
pNetCopy->Dump();
#endif
        //
        GenealogicalNetworkNode *pnodeDest = pBranchToChange->GetDestNode();
        GenealogicalNetworkNode *pnodeSrc = pBranchToChange->GetSrcNode();
        YW_ASSERT_INFO(mapOldNodeToNew.find(pnodeDest) != mapOldNodeToNew.end(), "Fail to find1");
        YW_ASSERT_INFO(mapOldNodeToNew.find(pnodeSrc) != mapOldNodeToNew.end(), "Fail to find2");
        GenealogicalNetworkNode *pnodeSrcCopy = mapOldNodeToNew[pnodeSrc];
        GenealogicalNetworkNode *pnodeDestCopy = mapOldNodeToNew[pnodeDest];
        
        set<GenealogicalNetworkBranch *> setAncBrs;
        YW_ASSERT_INFO( pnodeSrcCopy != NULL, "Wrong: null for srcopy branch" );
        pNetCopy->GetAncesBranchesFrom(pnodeSrcCopy, setAncBrs);
        
        // now re-attach
        GenealogicalNetworkBranch *pBrSrc = pNetCopy->GetBranch( pnodeSrcCopy,  pnodeDestCopy);
        YW_ASSERT_INFO(pBrSrc != NULL, "Fail to find pBrSrc here");
        YW_ASSERT_INFO(mapOldNodeToNew.find((*it)->GetSrcNode() ) != mapOldNodeToNew.end(), "Fail to find3");
        YW_ASSERT_INFO(mapOldNodeToNew.find((*it)->GetDestNode() ) != mapOldNodeToNew.end(), "Fail to find4");
        GenealogicalNetworkBranch *pBrDestCopy = pNetCopy->GetBranch( mapOldNodeToNew[(*it)->GetSrcNode()],  mapOldNodeToNew[(*it)->GetDestNode()]);
        YW_ASSERT_INFO(pBrDestCopy != NULL, "Fail to find5");
        GenealogicalNetworkNode *pnParAttached = pNetCopy->ReattachBranchTo( pBrSrc, pBrDestCopy );
        
        YW_ASSERT_INFO( pnParAttached != NULL, "Wrong: null for pnParAttached" );
        pNetCopy->GetAncesBranchesFrom(pnParAttached, setAncBrs);
        
        
        // ensure the network is OK
        pNetCopy->RemoveTrivialCycles();
        
        pNetCopy->RemoveIsolatedNodes();
        
        YW_ASSERT_INFO( pNetCopy->CheckSingleRoot(), "Bad net: more than 1 roots1" );
        
#if 0
cout << "*****************Constructed one neighbor network: ";
pNetCopy->Dump();
#endif
        listNgbrNetsDist1.push_back(pNetCopy);
        
        if( pListCriticalBrs != NULL )
        {
//cout << "Now prune non-exist brs...\n";
            // in case some branches are not there...
            pNetCopy->PruneNonexistBrs(setAncBrs);
            pListCriticalBrs->push_back(setAncBrs);
//cout <<"prune: done\n";
#if 0
cout << "**** Prune non-exist bras:neighbor network: ";
pNetCopy->Dump();
#endif
        }
    }
    // YW: this does not seem to work well (01/29/16); SkIP for now
//#if 0
    // finally also allow the branch to attach as sibling of the root (unless it is already incident to the root)
    // as before, don't allow mising edge to attach to root sib
    if( pBranchToChange->GetSrcNode()->IsRoot() == false  && pBranchToChange->IsSrcMixing() == false  )
    {
//cout << "Now consider root sib..\n";
        if( CheckSingleRoot() == false  )
        {
            cout << "*****************BAD SELF  ";
            Dump();
            YW_ASSERT_INFO(false, "Stop000");
        }
        
        //
        map<GenealogicalNetworkNode *, GenealogicalNetworkNode*> mapOldNodeToNew;
        GenealogicalNetwork *pNetCopy = Copy( &mapOldNodeToNew );
#if 0
        cout << "Network copy: ";
        pNetCopy->Dump();
#endif
        if( pNetCopy->CheckSingleRoot() == false  )
        {
            cout << "*****************BAD copy  ";
            pNetCopy->Dump();
            cout << "---------SELF  ";
            Dump();
            YW_ASSERT_INFO(false, "Stop00");
        }
        
        //
        GenealogicalNetworkNode *pnodeDest = pBranchToChange->GetDestNode();
        GenealogicalNetworkNode *pnodeSrc = pBranchToChange->GetSrcNode();
        YW_ASSERT_INFO(mapOldNodeToNew.find(pnodeDest) != mapOldNodeToNew.end(), "Fail to find1");
        YW_ASSERT_INFO(mapOldNodeToNew.find(pnodeSrc) != mapOldNodeToNew.end(), "Fail to find2");
        GenealogicalNetworkNode *pnodeSrcCopy = mapOldNodeToNew[pnodeSrc];
        GenealogicalNetworkNode *pnodeDestCopy = mapOldNodeToNew[pnodeDest];
        GenealogicalNetworkBranch *pBrSrc = pNetCopy->GetBranch( pnodeSrcCopy,  pnodeDestCopy);
        YW_ASSERT_INFO( pBrSrc!=NULL, "Fail to find branch11" );
        pNetCopy->ReattachBranchAsRootSib( pBrSrc );
        
        if( pNetCopy->CheckSingleRoot() == false  )
        {
            cout << "*****************BEFORE removing cycle  ";
            pNetCopy->Dump();
            YW_ASSERT_INFO(false, "Stop0");
        }
        
        // ensure the network is OK
        pNetCopy->RemoveTrivialCycles();
        
        // clean up
        pNetCopy->RemoveIsolatedNodes();
        if( pNetCopy->CheckSingleRoot() == false  )
        {
            cout << "*****************Constructed one neighbor network by attaching to the root sib: ";
            pNetCopy->Dump();
            YW_ASSERT_INFO(false, "Stop");
        }
        
        listNgbrNetsDist1.push_back(pNetCopy);
        if( pListCriticalBrs != NULL )
        {
            // for now, don't create critical edges for these network
            set<GenealogicalNetworkBranch *> setAncBrs;
            pListCriticalBrs->push_back(setAncBrs);
        }
#if 0
cout << "*****************Constructed one neighbor network by attaching to the root sib: ";
pNetCopy->Dump();
#endif
    }
//cout << "Done for this branch: num of networks so far: " << listNgbrNetsDist1.size() << endl;
//#endif
}

void GenealogicalNetwork :: RetriveNgbrNetsOneEvtForMixEdge( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetwork *> &listNgbrNetsDist1, vector<set<GenealogicalNetworkBranch *> > *pListCriticalBrs ) const
{
// hack: don't do this
//return;
//cout << "*** RetriveNgbrNetsOneEvtForMixEdge: network is: \n";
//this->DumpMargTrees();
//cout << "Edge: ";
//pBranchToChange->Dump();
    
    // find neighbors based on an admixed edge
    GenealogicalNetworkNode *pnodeMix=pBranchToChange->GetSrcNode();
    GenealogicalNetworkNode *pnodeMixDesc=pBranchToChange->GetDestNode();
    YW_ASSERT_INFO(pnodeMix->IsMixNode(), "Fatal error");
    GenealogicalNetworkBranch *pmixBr1=pnodeMix->GetAnces1();
    GenealogicalNetworkBranch *pmixBr2=pnodeMix->GetAnces2();
    YW_ASSERT_INFO(pmixBr1!=NULL && pmixBr2!=NULL, "Cannot be null ancestors at mix");
    GenealogicalNetworkNode *pMixAncNode1=pmixBr1->GetSrcNode();
    GenealogicalNetworkNode *pMixAncNode2=pmixBr2->GetSrcNode();
    YW_ASSERT_INFO( pMixAncNode1!=NULL && pMixAncNode2!=NULL, "Cannot be null" );
    //
    vector<GenealogicalNetworkNode *> listAncNodestoSwap;
    bool fAnc1=pMixAncNode1->IsAncestralTo(pMixAncNode2);
    bool fAnc2=pMixAncNode2->IsAncestralTo(pMixAncNode1);
    if(fAnc1==false && fAnc2==false)
    {
        // try both way of changing
        listAncNodestoSwap.push_back( pMixAncNode1 );
        listAncNodestoSwap.push_back( pMixAncNode2 );
    }
    else
    {
#if 0  // HACK 07/11/24
        GenealogicalNetworkNode *pMixAncNode=pMixAncNode1;
        if( fAnc1 == true )
        {
            YW_ASSERT_INFO(fAnc2==false, "Wrong");
            pMixAncNode=pMixAncNode2;
        }
        if(fAnc2 == true )
        {
            YW_ASSERT_INFO(fAnc1==false, "Wrong");
        }
        listAncNodestoSwap.push_back(pMixAncNode);
#endif
    }
    
    for(int i=0; i<(int)listAncNodestoSwap.size(); ++i)
    {
        GenealogicalNetworkNode *pMixAncNode = listAncNodestoSwap[i];

        //
        GenealogicalNetworkNode *pnodeOther = pMixAncNode->GetOtherDescNode(pnodeMix);
        if( pnodeOther==NULL)
        {
            //return;
            // YW: is this OK??? 06/03/24
            continue;
            //return;
        }
        
        // July 11,24: forbid change if there is ancestral relation between pnodeOther and desc
        set<GenealogicalNetworkNode *> setDescsBelow1, setDescsBelow2;
        pnodeMixDesc->GetAllDescendantNodesUnder( setDescsBelow1);
        pnodeOther->GetAllDescendantNodesUnder( setDescsBelow2);
        //if( pnodeOther->IsAncestralTo(pnodeMixDesc) || pnodeMixDesc->IsAncestralTo(pnodeOther) )
        if(AreSetsIntersecting(setDescsBelow1, setDescsBelow2) )
        {
            continue;
        }
        
        map<GenealogicalNetworkNode *, GenealogicalNetworkNode*> mapOldNodeToNew;
        GenealogicalNetwork *pNetCopy = Copy( &mapOldNodeToNew );
        YW_ASSERT_INFO(mapOldNodeToNew.find( pnodeMix ) != mapOldNodeToNew.end(), "Fail to find3");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pnodeMixDesc ) != mapOldNodeToNew.end(), "Fail to find4");
        
        GenealogicalNetworkBranch *pBrMixCopy = pNetCopy->GetBranch( mapOldNodeToNew[ pnodeMix ],  mapOldNodeToNew[ pnodeMixDesc ]);
        YW_ASSERT_INFO(pBrMixCopy != NULL, "NULL001");
        pBrMixCopy->DetachFromSrc();
        YW_ASSERT_INFO(mapOldNodeToNew.find( pMixAncNode ) != mapOldNodeToNew.end(), "Fail to find3");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pnodeOther ) != mapOldNodeToNew.end(), "Fail to find4");
        GenealogicalNetworkBranch *pBrOtherCopy = pNetCopy->GetBranch( mapOldNodeToNew[ pMixAncNode ],  mapOldNodeToNew[ pnodeOther ]);
        YW_ASSERT_INFO(pBrOtherCopy != NULL, "NULL002");
        pBrOtherCopy->DetachFromSrc();
        pBrMixCopy->AttachToSrc( mapOldNodeToNew[ pMixAncNode ], mapOldNodeToNew[ pnodeMixDesc ] );
        pBrOtherCopy->AttachToSrc( mapOldNodeToNew[ pnodeMix ], mapOldNodeToNew[pnodeOther] );
        
        GenealogicalNetworkNode *pnodeMixDescNew = mapOldNodeToNew[pnodeMixDesc];
        set<GenealogicalNetworkBranch *> setAncBrs;
        YW_ASSERT_INFO( pnodeMixDescNew != NULL, "mix node new: cannot be NULL" );
        pNetCopy->GetAncesBranchesFrom(pnodeMixDescNew, setAncBrs);
        
        // ensure the network is OK
        pNetCopy->RemoveTrivialCycles();
        
#if 0
        cout << "*****************Constructed one neighbor network by attaching to the root sib: ";
        pNetCopy->Dump();
        //pNetCopy->DumpMargTrees();
#endif
        
        pNetCopy->RemoveIsolatedNodes();
        YW_ASSERT_INFO( pNetCopy->CheckSingleRoot(), "Bad net: more than 1 roots3" );
        listNgbrNetsDist1.push_back(pNetCopy);
        
        if( pListCriticalBrs != NULL )
        {
            // in case some branches are not there...
            pNetCopy->PruneNonexistBrs(setAncBrs);
            pListCriticalBrs->push_back(setAncBrs);
        }
    }
//cout << "Done for mixing branch: num of networks so far: " << listNgbrNetsDist1.size() << endl;
}

void GenealogicalNetwork :: RetriveNgbrNetsOneEvtChangeMixEdge( GenealogicalNetworkBranch *pMixBrToChange, vector<GenealogicalNetwork *> &listNgbrNetsDist1, int outgroup, vector<set<GenealogicalNetworkBranch *> > *pListCriticalBrs ) const
{
    // change this mixing branch to create a new admix node below it;
    // that is, create a new admixture node
    GenealogicalNetworkNode *pnodeMix=pMixBrToChange->GetDestNode();
    YW_ASSERT_INFO(pnodeMix->IsMixNode(), "Fatal error");
    GenealogicalNetworkNode *pnodeMixSrc=pMixBrToChange->GetSrcNode();
    GenealogicalNetworkNode *pnodeMixDest = pnodeMix->GetChildSingle();
    GenealogicalNetworkBranch *pmixBr1=pnodeMix->GetAnces1();
    GenealogicalNetworkBranch *pmixBr2=pnodeMix->GetAnces2();
    YW_ASSERT_INFO(pmixBr1!=NULL && pmixBr2!=NULL, "Cannot be null ancestors at mix");
    GenealogicalNetworkBranch *pmixBrOther = pmixBr1;
    if( pmixBr1 == pMixBrToChange )
    {
        pmixBrOther = pmixBr2;
    }
    GenealogicalNetworkNode *pMixAncOther = pmixBrOther->GetSrcNode();
    
    // find all places to add the new mix: must be away from the current mixing location
    vector<GenealogicalNetworkNode *> listDecNodesNewMixAnc, listDecNodesNewMixDesc;

    // try all branches
    vector<GenealogicalNetworkBranch *> setAllBranches;
    GetAllBranchesOrdered2( setAllBranches);
    for(unsigned int i=0; i<setAllBranches.size(); ++i)
    {
        GenealogicalNetworkBranch *pbr = setAllBranches[i];
        GenealogicalNetworkNode *pbrSrc = pbr->GetSrcNode();
        GenealogicalNetworkNode *pbrDest = pbr->GetDestNode();
        if( pbr->IsAdjacentToOutgroup(outgroup) )
        {
            continue;
        }
        if( pbrDest->IsAncestralTo(pnodeMixSrc) ||  pnodeMix->IsAncestralTo(pbrSrc) || pbrDest == pnodeMix )
        {
            // this would lead to cycle...
            continue;
        }
        listDecNodesNewMixAnc.push_back(pbrSrc);
        listDecNodesNewMixDesc.push_back(pbrDest);
    }

    //
    int numNetsAdded = 0;
    for(int i=0; i<(int)listDecNodesNewMixAnc.size(); ++i)
    {
        GenealogicalNetworkNode *pMixAncNodeNew = listDecNodesNewMixAnc[i];
        GenealogicalNetworkNode *pMixDestNodeNew = listDecNodesNewMixDesc[i];
        
        map<GenealogicalNetworkNode *, GenealogicalNetworkNode*> mapOldNodeToNew;
        GenealogicalNetwork *pNetCopy = Copy( &mapOldNodeToNew );
        YW_ASSERT_INFO(mapOldNodeToNew.find( pnodeMix ) != mapOldNodeToNew.end(), "Fail to find3");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pnodeMixSrc ) != mapOldNodeToNew.end(), "Fail to find4");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pnodeMixDest ) != mapOldNodeToNew.end(), "Fail to find4.1");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pMixAncOther ) != mapOldNodeToNew.end(), "Fail to find4.2");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pMixAncNodeNew ) != mapOldNodeToNew.end(), "Fail to find4.3");
        YW_ASSERT_INFO(mapOldNodeToNew.find( pMixDestNodeNew ) != mapOldNodeToNew.end(), "Fail to find4.4");
//cout << "1.1\n";
        // detach mix desc branch
        GenealogicalNetworkBranch *pBrMixDescCopy = pNetCopy->GetBranch( mapOldNodeToNew[ pnodeMix ],  mapOldNodeToNew[ pnodeMixDest ]);
        YW_ASSERT_INFO(pBrMixDescCopy != NULL, "NULL001");
        pBrMixDescCopy->DetachFromSrc();
//cout << "1.2\n";
        // orig target branch, now enter from the other side of the mixing node
        GenealogicalNetworkBranch *pBrTargetCopy = pNetCopy->GetBranch( mapOldNodeToNew[ pMixAncNodeNew ],  mapOldNodeToNew[ pMixDestNodeNew ]);
        YW_ASSERT_INFO(pBrTargetCopy!=NULL, "NULL003");
        pBrTargetCopy->DetachFromSrc();
//cout << "1.3\n";
        // the other mix branch now connects to mix dest
        GenealogicalNetworkBranch *pBrMixAncCopy = pNetCopy->GetBranch( mapOldNodeToNew[ pMixAncOther ],  mapOldNodeToNew[ pnodeMix ]);
        YW_ASSERT_INFO(pBrMixAncCopy != NULL, "NULL001");
        pBrMixAncCopy->DetachFromSrc();
//cout << "1.4\n";
        
        // now re-attach
        pBrMixDescCopy->AttachToSrc( mapOldNodeToNew[ pMixAncOther ], mapOldNodeToNew[ pnodeMixDest ] );
//cout << "1.5\n";
        pBrTargetCopy->AttachToSrc( mapOldNodeToNew[pnodeMix], mapOldNodeToNew[pMixDestNodeNew] );
//cout << "1.6\n";
        pBrMixAncCopy->AttachToSrc( mapOldNodeToNew[pMixAncNodeNew], mapOldNodeToNew[pnodeMix]  );
//cout << "1.7\n";
//cout << "Network so far: ";
//pNetCopy->Dump();
        
        // ensure the network is OK
        pNetCopy->RemoveTrivialCycles();
//cout << "1.7.1\n";
        GenealogicalNetworkNode *pnodeMixDescNew = mapOldNodeToNew[pnodeMixDest];
        set<GenealogicalNetworkBranch *> setAncBrs;
        YW_ASSERT_INFO( pnodeMixDescNew != NULL, "mix node new: cannot be NULL" );
        pNetCopy->GetAncesBranchesFrom(pnodeMixDescNew, setAncBrs);
//cout << "1.8\n";
#if 0
        cout << "*****************Constructed one neighbor network by changing an admix edge: ";
        //pNetCopy->Dump();
        pNetCopy->DumpMargTrees(false);
#endif
        
        pNetCopy->RemoveIsolatedNodes();
        YW_ASSERT_INFO( pNetCopy->CheckSingleRoot(), "Bad net: more than 1 roots3" );
        listNgbrNetsDist1.push_back(pNetCopy);
//cout << "1.9\n";
        if( pListCriticalBrs != NULL )
        {
            // in case some branches are not there...
            pNetCopy->PruneNonexistBrs(setAncBrs);
            pListCriticalBrs->push_back(setAncBrs);
        }
//cout << "1.10\n";
        ++numNetsAdded;
    }
//cout << "Done for changing admixture branch: num of networks added: " << numNetsAdded << endl;
}

void GenealogicalNetwork :: RetriveNgbrNetsOneNewMixForEdge( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetwork *> &listNetsOneNewMix ) const
{
//cout << "GenealogicalNetwork :: RetriveNgbrNetsOneNewMixForEdge: current net: ";
//this->DumpMargTrees(false);
//cout << "Edge: src: " << pBranchToChange->GetSrcNode()->GetTaxonId() << ", dest: " << pBranchToChange->GetDestNode()->GetTaxonId() << endl;
    
    // don't work on leaf branch (unless its leaf branch)
    if(pBranchToChange->IsLeafBranch()==false)
    {
        return;
    }
    
    // change where the branch attaches to and return all the possible networks
    vector<GenealogicalNetworkBranch *> setBrsToReattach;
    FindReattachBranches( pBranchToChange, setBrsToReattach, false );   // for new mixing, avoid NNI: need to search for more choices
//cout << "setBrsToReattach: size " << setBrsToReattach.size() << endl;
    // YW: 02/02/16: the number of mix nodes of new network should work better
    int numNetNodes = GetNumMixNodes();
    
    for( vector<GenealogicalNetworkBranch *> :: iterator it = setBrsToReattach.begin(); it != setBrsToReattach.end(); ++it )
    {
#if 0
        cout << "Now reattach to this branch: ";
        (*it)->Dump();
#endif
        map<GenealogicalNetworkNode *, GenealogicalNetworkNode*> mapOldNodeToNew;
        GenealogicalNetwork *pNetCopy = Copy( &mapOldNodeToNew );
#if 0
        cout << "Network copy: ";
        pNetCopy->DumpMargTrees(false);
        //pNetCopy->Dump();
#endif
        //
        GenealogicalNetworkNode *pnodeDest = pBranchToChange->GetDestNode();
        GenealogicalNetworkNode *pnodeSrc = pBranchToChange->GetSrcNode();
        YW_ASSERT_INFO(mapOldNodeToNew.find(pnodeDest) != mapOldNodeToNew.end(), "Fail to find1");
        YW_ASSERT_INFO(mapOldNodeToNew.find(pnodeSrc) != mapOldNodeToNew.end(), "Fail to find2");
        GenealogicalNetworkNode *pnodeSrcCopy = mapOldNodeToNew[pnodeSrc];
        GenealogicalNetworkNode *pnodeDestCopy = mapOldNodeToNew[pnodeDest];
        GenealogicalNetworkBranch *pBrSrc = pNetCopy->GetBranch( pnodeSrcCopy,  pnodeDestCopy);
        YW_ASSERT_INFO(mapOldNodeToNew.find((*it)->GetSrcNode() ) != mapOldNodeToNew.end(), "Fail to find3");
        YW_ASSERT_INFO(mapOldNodeToNew.find((*it)->GetDestNode() ) != mapOldNodeToNew.end(), "Fail to find4");
        GenealogicalNetworkBranch *pBrDestCopy = pNetCopy->GetBranch( mapOldNodeToNew[(*it)->GetSrcNode()],  mapOldNodeToNew[(*it)->GetDestNode()]);
        YW_ASSERT_INFO(pBrDestCopy != NULL, "Fail to find5");
        pNetCopy->AddMixBtwBranches( pBrSrc, pBrDestCopy );
        
        // ensure the network is OK
        pNetCopy->RemoveTrivialCycles();
        pNetCopy->RemoveIsolatedNodes();
        YW_ASSERT_INFO( pNetCopy->CheckSingleRoot(), "Bad net: more than 1 roots4" );
        
#if 0
        cout << "*Constructed one neighbor network: ";
        pNetCopy->DumpMargTrees(false);
        pNetCopy->Dump();
#endif
        
        // YW: 02/02/16: only accept the new candidate when the number of mix nodes of new network is larger
        // YW: is this really the good choice? Not sure yet
        if( pNetCopy->GetNumMixNodes() > numNetNodes )
        {
            listNetsOneNewMix.push_back(pNetCopy);
        }
        else
        {
            delete pNetCopy;
        }
    }

}

void GenealogicalNetwork :: FindReattachBranches( GenealogicalNetworkBranch *pBranchToChange, vector<GenealogicalNetworkBranch *> &setBranchesToAttach, bool fUseNNIMode ) const
{
    YW_ASSERT_INFO(pBranchToChange!=NULL, "NULL113");
    // if this branch's source node is mixing, don't allow it to change
    if( pBranchToChange->IsSrcMixing() == true )
    {
//cout << "Already mixing..\n";
        return;
    }
    
    
    // first, cannot use any edges that are descendent to the current edge
    vector<GenealogicalNetworkBranch *> setAllBranches;
    //GetAllBranches( setAllBranches );
    // Use ordered branches to remove randomness
    //GetAllBranchesBottomUp( setAllBranches );
    // YW: May 25, 2024
    GetAllBranchesOrdered2( setAllBranches );
    vector<GenealogicalNetworkBranch *> setUsableBrs;
    for( vector<GenealogicalNetworkBranch *> :: iterator it = setAllBranches.begin(); it != setAllBranches.end(); ++it)
    {
        // also don't allow if a branch share nodes with the changed edge
        // and do not allow mixing edges
        GenealogicalNetworkNode *pNodeSrc = pBranchToChange->GetSrcNode();
        GenealogicalNetworkNode *pNodeDest = pBranchToChange->GetDestNode();
        YW_ASSERT_INFO(*it != NULL, "NULL111");
        // are we in NNI mode
        if( fUseNNIMode == false )
        {
            //if( pBranchToChange->IsAncestralTo(*it) == false && (*it)->IsMixing() == false && (*it)->GetSrcNode() != pNodeSrc && (*it)->GetDestNode() != pNodeSrc && (*it)->GetSrcNode() != pNodeDest && (*it)->GetDestNode() != pNodeDest   )
            if( pBranchToChange->IsAncestralTo(*it) == false && (*it)->GetSrcNode() != pNodeSrc && (*it)->GetDestNode() != pNodeSrc && (*it)->GetSrcNode() != pNodeDest && (*it)->GetDestNode() != pNodeDest   )
            //if( pBranchToChange->IsAncestralTo(*it) == false && (*it)->IsMixing() == false && (*it)->GetSrcNode() != pNodeSrc && (*it)->GetDestNode() != pNodeSrc  )
            {
                //
                setUsableBrs.push_back( *it );
            }
        }
        else
        {
            // only allow edges that are sibling of the current edge
            GenealogicalNetworkNode *pSrcOther = (*it)->GetSrcNode();
            //if( pBranchToChange->IsAncestralTo(*it) == false && pNodeSrc != NULL && pNodeSrc->IsMixNode() == false && pSrcOther != NULL  && pSrcOther->IsMixNode() == false &&  pNodeSrc->GetParentSingle() != NULL && pNodeSrc->GetParentSingle() == pSrcOther->GetParentSingle()  && (*it)->GetSrcNode() != pNodeSrc && (*it)->GetDestNode() != pNodeSrc && (*it)->GetSrcNode() != pNodeDest && (*it)->GetDestNode() != pNodeDest  )
            if( pBranchToChange->IsAncestralTo(*it) == false && pNodeSrc != NULL && pNodeSrc->IsMixNode() == false && pSrcOther != NULL  && pSrcOther->IsMixNode() == false &&  pNodeSrc->GetParentSingle() != NULL && pNodeSrc->GetParentSingle() == pSrcOther->GetParentSingle()  && (*it)->GetSrcNode() != pNodeSrc && (*it)->GetDestNode() != pNodeSrc && (*it)->GetSrcNode() != pNodeDest && (*it)->GetDestNode() != pNodeDest  )
            {
                setUsableBrs.push_back(*it);
            }
        }
    }
    setBranchesToAttach = setUsableBrs;
}

GenealogicalNetwork * GenealogicalNetwork :: Copy( map<GenealogicalNetworkNode *, GenealogicalNetworkNode *> *pMapOldNodeToNew ) const
{
//cout << "-- copying...\n";
    // Make a copy of the current network
    // first create all the nodes
    GenealogicalNetwork *pNetNew = new GenealogicalNetwork();
    
    map<GenealogicalNetworkNode *, GenealogicalNetworkNode *> mapOldNodeToNew;
    map<GenealogicalNetworkNode *, GenealogicalNetworkNode *> mapNewNodeToOld;
    set<GenealogicalNetworkNode *> setNodesDone;
    //for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it )
    for( vector<GenealogicalNetworkNode *> :: const_iterator it = listNodesAdded.begin(); it != listNodesAdded.end(); ++it )
    {
        if( setNodesDone.find(*it) != setNodesDone.end() || setNodes.find(*it) == setNodes.end())
        {
            continue;
        }
        setNodesDone.insert(*it);
//cout << "  copy node: ";
//(*it)->Dump();
//cout << endl;
        //
        GenealogicalNetworkNode *pnodeNew = new GenealogicalNetworkNode();
        pnodeNew->SetMixRatio( (*it)->GetMixRatio() );
        pnodeNew->SetTaxonId( (*it)->GetTaxonId() );
        pnodeNew->SetTaxonStrUser( (*it)->GetTaxonStrUser() );
        pNetNew->AddNode( pnodeNew );
        mapOldNodeToNew.insert( map<GenealogicalNetworkNode *, GenealogicalNetworkNode *> :: value_type( *it, pnodeNew ) );
        mapNewNodeToOld[pnodeNew] = *it;
    }
    // now create all the edges
    for( set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
    {
        GenealogicalNetworkNode *pSrc = (*it)->GetSrcNode();
        GenealogicalNetworkNode *pDest = (*it)->GetDestNode();
        GenealogicalNetworkBranch *pBrOld = GetBranch( (*it)->GetSrcNode(), (*it)->GetDestNode() );
        double lenOld = pBrOld->GetLength();
        YW_ASSERT_INFO( mapOldNodeToNew.find(pSrc) != mapOldNodeToNew.end(), "Fail to find" );
        YW_ASSERT_INFO( mapOldNodeToNew.find(pDest) != mapOldNodeToNew.end(), "Fail to find" );
        pNetNew->AddBranch( mapOldNodeToNew[pSrc], mapOldNodeToNew[pDest], lenOld );
    }
    
    // fix the order of left/right parents for admix nodes so that they match
    for( set<GenealogicalNetworkNode *> :: iterator it = pNetNew->setNodes.begin(); it != pNetNew->setNodes.end(); ++it )
    {
        GenealogicalNetworkNode *pn = *it;
        // maintain it if it is admix
        if( pn->IsMixNode() )
        {
            GenealogicalNetworkBranch *pbr1 = pn->GetAnces1();
            GenealogicalNetworkBranch *pbr2 = pn->GetAnces2();
            YW_ASSERT_INFO( pbr1 != NULL && pbr2 != NULL, "WRONG!!!" );
            GenealogicalNetworkNode *pbr1Src = pbr1->GetSrcNode();
            GenealogicalNetworkNode *pbr2Src = pbr2->GetSrcNode();
            map<GenealogicalNetworkNode *, GenealogicalNetworkNode *> :: iterator it22 = mapNewNodeToOld.find(pbr1Src);
            YW_ASSERT_INFO( it22 != mapNewNodeToOld.end(), "Fail to find 3344" );
            GenealogicalNetworkNode *pbr1SrcOld = it22->second;
            it22 = mapNewNodeToOld.find(pbr2Src);
            YW_ASSERT_INFO( it22 != mapNewNodeToOld.end(), "Fail to find 3344" );
            GenealogicalNetworkNode *pbr2SrcOld = it22->second;

            it22 = mapNewNodeToOld.find(pn);
            YW_ASSERT_INFO( it22 != mapNewNodeToOld.end(), "Fail to find 3344" );
            GenealogicalNetworkNode *pnOld = it22->second;
            GenealogicalNetworkBranch *pbr1Old = pnOld->GetAnces1();
            GenealogicalNetworkBranch *pbr2Old = pnOld->GetAnces2();
            YW_ASSERT_INFO( pbr1Old != NULL && pbr2Old != NULL, "WRONG!!!" );
            GenealogicalNetworkNode *pbr1OldSrc = pbr1Old->GetSrcNode();
            GenealogicalNetworkNode *pbr2OldSrc = pbr2Old->GetSrcNode();
            // swap if the old and new nodes in parents are swapped
            bool f1 = (pbr1SrcOld == pbr1OldSrc) && (pbr2SrcOld == pbr2OldSrc);
            bool f2 = (pbr1SrcOld == pbr2OldSrc) && (pbr2SrcOld == pbr1OldSrc);
            YW_ASSERT_INFO(f1 || f2, "FATAL ERROR: nodes mismatch");
            if( f1 == false )
            {
                // swap now
                pn->SwapAncesBrs();    // swap the left and right branches
            }
        }
    }
    
    if(pMapOldNodeToNew != NULL)
    {
        *pMapOldNodeToNew = mapOldNodeToNew;
    }
    
    return pNetNew;
}

GenealogicalNetworkNode * GenealogicalNetwork :: ReattachBranchTo( GenealogicalNetworkBranch *pBrSrc, GenealogicalNetworkBranch *pBrDest )
{
    // return the newly added parent node
#if 0
cout << "ReattachBranchTo: pBrSrc: ";
pBrSrc->Dump();
cout << "   , pBrDest: ";
pBrDest->Dump();
cout << endl;
#endif
    // attach the pbrsrc to be as sibling as pBrDest
    GenealogicalNetworkNode *pNodeSrcOfOrig = pBrSrc->GetSrcNode();
    GenealogicalNetworkNode *pNodeDestOfOrig = pBrSrc->GetDestNode();
    GenealogicalNetworkNode *pNodeSrcOfSib = pBrDest->GetSrcNode();
    GenealogicalNetworkNode *pNodeDestOfSib = pBrDest->GetDestNode();
    
    double brSrcOrigLen = pBrSrc->GetLength();
    
    // attach in the middle of the new edge
    pBrSrc->DetachFromSrc();
    RemoveBranch( pNodeSrcOfOrig, pNodeDestOfOrig, true );
//cout << "After detach branch pBrSrc: " << endl;
    CleanupAt( pNodeSrcOfOrig );
#if 0
cout << "After cleanup by detaching pBrSrc: network becomes: ";
this->Dump();
#endif
    GenealogicalNetworkBranch *pBrToAttach = GetBranch(pNodeSrcOfSib, pNodeDestOfSib);
    double brOld = pBrToAttach->GetLength();
    pBrDest->DetachFromSrc();
    RemoveBranch( pNodeSrcOfSib, pNodeDestOfSib, true );
    GenealogicalNetworkBranch *pBrNew1=NULL;
    GenealogicalNetworkBranch *pBrNew2=NULL;
    GenealogicalNetworkNode *pnParNew = pNodeDestOfOrig->AddToSib( pNodeDestOfSib, brSrcOrigLen, 0.5*brOld, pBrNew1, pBrNew2 );
    this->AddNode( pnParNew );
    this->AddBranch( pNodeSrcOfSib, pnParNew, 0.5*brOld );
    YW_ASSERT_INFO( pBrNew1!=NULL&&pBrNew2 !=NULL, "Fail to add" );
    this->setBranches.insert( pBrNew1 );
    this->setBranches.insert( pBrNew2 );
    return pnParNew;
}

void GenealogicalNetwork :: ReattachBranchAsRootSib(GenealogicalNetworkBranch *pBrSrc)
{
#if 0
    cout << "ReattachBranchAsRootSib: pBrSrc: ";
    pBrSrc->Dump();
    cout << endl;
#endif
    // attach the pbrsrc to be as sibling as pBrDest
    GenealogicalNetworkNode *pNodeSrcOfOrig = pBrSrc->GetSrcNode();
    GenealogicalNetworkNode *pNodeDestOfOrig = pBrSrc->GetDestNode();
    GenealogicalNetworkNode *pNodeRootCur = GetRoot();
    
    double brSrcOrigLen = pBrSrc->GetLength();
    
    // attach in the middle of the new edge
    pBrSrc->DetachFromSrc();
    RemoveBranch( pNodeSrcOfOrig, pNodeDestOfOrig, true );
    //cout << "After detach branch pBrSrc: " << endl;
    CleanupAt( pNodeSrcOfOrig );
#if 0
    cout << "After cleanup by detaching pBrSrc: network becomes: ";
    this->Dump();
#endif
    // new branch to root is set to a default value of 0.1
    const double LEN_DEF_ROOT_BRANCH = 0.1;
    double brOld = LEN_DEF_ROOT_BRANCH;
    GenealogicalNetworkBranch *pBrNew1=NULL;
    GenealogicalNetworkBranch *pBrNew2=NULL;
    GenealogicalNetworkNode *pNodeRootNew = pNodeDestOfOrig->AddToSib( pNodeRootCur, brSrcOrigLen, brOld, pBrNew1, pBrNew2 );
    this->AddNode( pNodeRootNew );
    YW_ASSERT_INFO( pBrNew1!=NULL&&pBrNew2 !=NULL, "Fail to add" );
    this->setBranches.insert( pBrNew1 );
    this->setBranches.insert( pBrNew2 );
}

void GenealogicalNetwork :: AddMixBtwBranches( GenealogicalNetworkBranch *pBrMixSink, GenealogicalNetworkBranch *pBrMixSource )
{
    //
#if 0
    cout << "AddMixBtwBranches: pBrMixSink: ";
    pBrMixSink->Dump();
    cout << "   , pBrMixSource: ";
    pBrMixSource->Dump();
    cout << endl;
#endif
    // attach the pbrsrc to be as sibling as pBrDest
    GenealogicalNetworkNode *pNodeSourceOfSink = pBrMixSink->GetSrcNode();
    GenealogicalNetworkNode *pNodeDestOfSink = pBrMixSink->GetDestNode();
    GenealogicalNetworkNode *pNodeSrcOfSource = pBrMixSource->GetSrcNode();
    GenealogicalNetworkNode *pNodeDestOfSource = pBrMixSource->GetDestNode();
    
    double brSrcSinkLen = pBrMixSink->GetLength();
    double brSrcSourceLen = pBrMixSource->GetLength();
    
    // attach in the middle of the new edge
    pBrMixSink->DetachFromSrc();
    RemoveBranch( pNodeSourceOfSink, pNodeDestOfSink, true );
    GenealogicalNetworkNode *pnMiddleNew = new GenealogicalNetworkNode;
    const double DEF_MIX_RATIO = 0.5;
    pnMiddleNew->SetMixRatio(DEF_MIX_RATIO);
    this->AddNode( pnMiddleNew );
    this->AddBranch( pNodeSourceOfSink, pnMiddleNew, 0.5*brSrcSinkLen );
    this->AddBranch( pnMiddleNew, pNodeDestOfSink, 0.5*brSrcSinkLen );

#if 0
    cout << "After cleanup by detaching pBrSrc: network becomes: ";
    this->Dump();
#endif
    //GenealogicalNetworkBranch *pBrToAttach = GetBranch(pNodeSrcOfSource, pNodeDestOfSource);
    pBrMixSource->DetachFromSrc();
    RemoveBranch( pNodeSrcOfSource, pNodeDestOfSource, true );
    GenealogicalNetworkBranch *pBrNew1=NULL;
    GenealogicalNetworkBranch *pBrNew2=NULL;
    GenealogicalNetworkNode *pnParNew = pnMiddleNew->AddToSib( pNodeDestOfSource, brSrcSinkLen, 0.5*brSrcSourceLen, pBrNew1, pBrNew2 );
    this->AddNode( pnParNew );
    this->AddBranch( pNodeSrcOfSource, pnParNew, 0.5*brSrcSourceLen );
    YW_ASSERT_INFO( pBrNew1!=NULL&&pBrNew2 !=NULL, "Fail to add" );
    this->setBranches.insert( pBrNew1 );
    this->setBranches.insert( pBrNew2 );

}

void GenealogicalNetwork :: MapNetCodedEdgeToMTBrs(const GenealogicalNetworkMTreeCode &ncode, const MarginalTree &margTree, std::map<GenealogicalNetworkBranch *, int> &mapChosenEdgeToMargBr ) const
{
//cout << "----net code: ";
//ncode.Dump();
    //
    mapChosenEdgeToMargBr.clear();
    
    // get all chosen edges
    set<GenealogicalNetworkBranch *> setChosenBrs;
    ncode.GetChosenBrs(setChosenBrs);
    // add in all tree edges of the net
    for(auto x: this->setBranches)
    {
        if( x->IsMixing() == false )
        {
            setChosenBrs.insert(x);
        }
    }
    
    // collect all splits associated with these chosen brs
    vector<GenealogicalNetworkBranch *> setAllBranches;
    GetAllBranchesBottomUp( setAllBranches );
    map<GenealogicalNetworkBranch *, set<int>  > mapBrSplit;
    for(unsigned int i=0; i<setAllBranches.size(); ++i)
    {
        GenealogicalNetworkBranch *pb = setAllBranches[i];
        if( setChosenBrs.find(pb) == setChosenBrs.end())
        {
            continue;
        }
        if( pb->IsLeafBranch() )
        {
            //
            int tid = pb->GetDestNode()->GetTaxonId();
            mapBrSplit[pb].insert(tid);
        }
        else
        {
            // get all descendenats of this branch
            GenealogicalNetworkBranch *pd1 = pb->GetDestNode()->GetLeftDesc();
            if( pd1 != NULL  )
            {
                UnionSets( mapBrSplit[pb], mapBrSplit[pd1] );
            }
            GenealogicalNetworkBranch *pd2 = pb->GetDestNode()->GetRightDesc();
            if( pd2 != NULL  )
            {
                UnionSets( mapBrSplit[pb], mapBrSplit[pd2] );
            }
#if 0
cout << "Branch: ";
pb->Dump();
cout << "  clade: ";
DumpIntSet(mapBrSplit[pb]);
#endif
        }
    }
    
    // get all clusters by marginal tree, for each
    vector<set<int> > listSplits;
    margTree.FindAllSplits(listSplits);
//cout << "***Set of clades in mtree: \n";
//for(unsigned int i=0; i<listSplits.size(); ++i)
//{
//DumpIntSet(listSplits[i]);
//}
    map<set<int>, int> mapSplitToBr;
    for(unsigned int i=0; i<listSplits.size(); ++i)
    {
        mapSplitToBr[listSplits[i]] = i;
    }
    
    // create the mapping
    for(auto x: setChosenBrs)
    {
        set<int> clade = mapBrSplit[x];
        int mid = -1;
        if( clade.size() > 0 )
        {
            YW_ASSERT_INFO( mapSplitToBr.find(clade) != mapSplitToBr.end(), "Fail to find clade" );
            mid = mapSplitToBr[clade];
        }
        mapChosenEdgeToMargBr[x] = mid;
    }
}

void GenealogicalNetwork :: CleanupAt( GenealogicalNetworkNode *pNodeChanged )
{
    if( pNodeChanged == NULL )
    {
        return;
    }
#if 0
cout << "Cleanu: pNodeChanged: ";
pNodeChanged->Dump();
cout << endl;
#endif
    // after changing the network, clean up the structure at this node
    // if this node has no taxon and there is no children
    GenealogicalNetworkNode *pNodeAnc1 = NULL;
    if( pNodeChanged->GetAnces1() != NULL )
    {
        pNodeAnc1 = pNodeChanged->GetAnces1()->GetSrcNode();
    }
    GenealogicalNetworkNode *pNodeAnc2 = NULL;
    if( pNodeChanged->GetAnces2() != NULL )
    {
        pNodeAnc2 = pNodeChanged->GetAnces2()->GetSrcNode();
    }
    
    if( pNodeChanged->IsLeaf() == true )
    {
        if( pNodeChanged->HasTaxon() == false )
        {
//cout << "Cleanup at a leaf\n";
            // erase it and then cleanup each of its ancestor
            RemoveNode(pNodeChanged);
            CleanupAt( pNodeAnc1 );
            CleanupAt( pNodeAnc2 );
        }
    }
    else if( pNodeChanged->GetNumChildren() == 1 )
    {
        GenealogicalNetworkNode *pChildSingle = pNodeChanged->GetChildSingle();
#if 0
cout << "In cleanup: pChildSingle: ";
pChildSingle->Dump();
cout << endl;
#endif
        // if there are two paernts, don't do anything
        if( pNodeChanged->GetNumParents() <= 1 )
        {
            GenealogicalNetworkBranch *pBrOld = GetBranch(pNodeChanged, pChildSingle);
            YW_ASSERT_INFO(pBrOld!=NULL, "NULL");
            double brLenOld = pBrOld->GetLength();
            
            pNodeChanged->RemoveChild(pChildSingle);
            RemoveBranch( pNodeChanged, pChildSingle, true);
            //cout << "Remove child1\n";
            
            // remove it
            GenealogicalNetworkNode *pParOrig = pNodeChanged->GetParentSingle();
            if( pParOrig != NULL )
            {
                GenealogicalNetworkBranch *pBrOld2 = GetBranch(pParOrig, pNodeChanged);
                YW_ASSERT_INFO(pBrOld2 != NULL, "NULL");
                brLenOld += pBrOld2->GetLength();
                pParOrig->RemoveChild( pNodeChanged );
//cout << "Remove child2\n";
                RemoveBranch( pParOrig, pNodeChanged, true );
                AddBranch( pParOrig, pChildSingle, brLenOld );
            }
            RemoveNode( pNodeChanged, false );
            //delete pNodeChanged;
        }
    }
//cout << "Cleanup: done\n";
}

void GenealogicalNetwork :: RemoveBranch( GenealogicalNetworkNode *pSrcNode, GenealogicalNetworkNode *pDestNode, bool fMem )
{
    //
    for( set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
    {
        if( (*it)->GetSrcNode() == pSrcNode && (*it)->GetDestNode() == pDestNode )
        {
            if( fMem == true )
            {
                delete *it;
            }
            setBranches.erase( *it);
            break;
        }
    }
}

void GenealogicalNetwork :: RemoveNode(GenealogicalNetworkNode *pNode, bool fRmBrs )
{
    // remove node and all attached edges
    if( fRmBrs)
    {
        set<GenealogicalNetworkBranch *> setBrsLeft;
        set<GenealogicalNetworkBranch *>setBrsRm;
        for( set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
        {
            if( (*it)->GetSrcNode() != pNode && (*it)->GetDestNode() != pNode )
            {
                setBrsLeft.insert(*it);
            }
            else
            {
                setBrsRm.insert(*it);
            }
        }
        setBranches = setBrsLeft;
        for( set<GenealogicalNetworkBranch *> :: iterator it = setBrsRm.begin(); it != setBrsRm.end(); ++it)
        {
            setBranches.erase(*it);
            delete *it;
        }
    }
    
    setNodes.erase(pNode);
    delete pNode;
}

void GenealogicalNetwork :: MapBackUserLabels(TaxaMapper &mapperTaxa)
{
    // process each node
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it )
    {
        //
        int idTaxaInternal = (*it)->GetTaxonId();
        if( idTaxaInternal < 0)
        {
            // this is an internal node so skip
            continue;
        }
        string strUserLbl = mapperTaxa.GetString( idTaxaInternal );
        (*it)->SetTaxonStrUser( strUserLbl );
    }
}

void GenealogicalNetwork :: MapBackUserLabelsByUserLabels(TaxaMapper &mapperTaxa)
{
    // ASSUME: the labels are in integers format
    // process each node
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it )
    {
        //
        string idTaxaInternalLbl = (*it)->GetTaxonStrUser();
        if( idTaxaInternalLbl.length() == 0)
        {
            // this is an internal node so skip
            continue;
        }
        int idTaxaInternal = atoi(idTaxaInternalLbl.c_str() );
        string strUserLbl = mapperTaxa.GetString( idTaxaInternal );
        (*it)->SetTaxonStrUser( strUserLbl );
    }
}

void GenealogicalNetwork :: OutputGML ( const char *inFileName )
{
    // map ID to id from 1, 2, 3, and so on
    int idToUse = 1;
    map<GenealogicalNetworkNode *,int> mapIDToBase1;
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it )
    {
        GenealogicalNetworkNode *pnodeCur = *it;
        mapIDToBase1[pnodeCur] = idToUse++;
    }
    
	// Now output a file in GML format
	// First create a new name
	string name = inFileName;
    //cout << "num edges = " << listEdges.size() << endl;
    
	// Now open file to write out
	ofstream outFile( name.c_str() );
    
	// First output some header info
	outFile << "graph [\n";
	outFile << "comment ";
	OutputQuotedString(outFile, "Population network, automatically generated by Graphing tool");
	outFile << "\ndirected  1\n";
	outFile << "id  1\n";
	outFile << "label ";
	OutputQuotedString ( outFile, "Population network....");
    outFile << endl;
    //outFile << endl;
    
	// Now output all the vertices
    //	int i;

    //cout << "a.1.1\n";
    for( set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it )
	{
        GenealogicalNetworkNode *pnodeCur = *it;

		outFile << "node [\n";
        
		outFile << "id " <<  mapIDToBase1[pnodeCur]  << endl;
		outFile << "label ";
 		string nameToUse = "";
        string strLblUse = pnodeCur->GetTaxonStrUser();
        if( strLblUse.length() == 0 )
        {
            int taxonId = pnodeCur->GetTaxonId();
            if( taxonId >= 0 )
            {
                char buf[100];
                snprintf(buf, 100,"%d", taxonId);
                nameToUse = buf;
            }
            else
            {
                // for mixing node, output mixing ratio
                if( pnodeCur->IsMixNode() == true )
                {
                    char buf[100];
                    snprintf(buf, 100,"Mix:%.2f", pnodeCur->GetMixRatio() );
                    nameToUse = buf;
                }
            }
        }
        else
        {
            nameToUse = strLblUse;
        }

        const char *name = nameToUse.c_str();
        OutputQuotedString (outFile,  name  );
		outFile << endl;
        
        // See if we need special shape for mixing node
        if( pnodeCur->IsMixNode() == true )
        {
            outFile << "vgj [ \n shape  ";
            OutputQuotedString( outFile, "Rectangle");
    		outFile << "\n]\n";
        }
        else
        {
		    outFile << "defaultAtrribute   1\n";
        }
        
		outFile << "]\n";
	}
    //cout << "a.1.3\n";
    
	// Now output all the edges, by again starting from root and output all nodes
    for( set<GenealogicalNetworkBranch *> :: iterator it = setBranches.begin(); it != setBranches.end(); ++it)
    {
        GenealogicalNetworkNode *pnsrc = (*it)->GetSrcNode();
        GenealogicalNetworkNode *pndest = (*it)->GetDestNode();
        //cout << "Output an edge \n";
        outFile << "edge [\n";
        outFile << "source " << mapIDToBase1[pnsrc] << endl;
        outFile << "target  " << mapIDToBase1[pndest] << endl;
        outFile << "label " ;
        
        // output branch length
        double len = (*it)->GetLength();
        char buf[100];
        snprintf(buf, 100, "%1.3f", len);
        
        OutputQuotedString( outFile, buf);
        outFile << "\n";
        outFile << "]\n";
	}
    
    
	// Finally quite after closing file
	outFile << "\n]\n";
	outFile.close();
}

void GenealogicalNetwork :: DumpMargTrees(bool fConv) const
{
    map<string,string> mapTaxonIdToUserId;
    CreateMapTaxonIdToUserId( mapTaxonIdToUserId );
//cout << "mapTaxonIdToUserId: ";
//for( map<string,string> ::iterator it = mapTaxonIdToUserId.begin(); it!= mapTaxonIdToUserId.end(); ++it )
//{
//cout << it->first   << " to "  << it->second << "  ";
//}
    
    // output all marginal trees
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTreesWithFreq;
    RetriveAllMarginalTrees( mapMargTreesWithFreq );
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo>::const_iterator it= mapMargTreesWithFreq.begin(); it != mapMargTreesWithFreq.end(); ++it )
    {
        // format: [freq] newick-format-tree
        cout << "[" << it->second.GetFreq() << "] ";
        string strNWTree = it->second.GetMargTree()->GetNewick() ;
//cout << "Before convert: tree is " << strNWTree << endl;
        if( fConv == true )
        {
            NewickUtils:: UpdateLabells( strNWTree, mapTaxonIdToUserId );
        }
//cout << "AFTER convert: tree is " << strNWTree << endl;
        cout << strNWTree << endl;
        it->second.FreeTree();
    }
    // also dump out admixture
    //DumpAdmixNodes();
    // free
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo>::const_iterator it= mapMargTreesWithFreq.begin(); it != mapMargTreesWithFreq.end(); ++it )
    {
        delete it->second.GetMargTree();
    }
}

void GenealogicalNetwork :: DumpAdmixNodes() const
{
    //
    // find all the admixture nodes
    vector<GenealogicalNetworkNode *> listMixNodes;
    GetMixNodes( listMixNodes );
    for(int i=0; i<(int)listMixNodes.size(); ++i)
    {
        cout << "Admixture population: ";
        set<string> setLvTaxaIds;
        listMixNodes[i]->GetAllLeafTaxonIdsUnder(setLvTaxaIds);
        if( setLvTaxaIds.size() == 1 )
        {
            cout << *( setLvTaxaIds.begin() ) << endl;
        }
        else if(setLvTaxaIds.size() > 1 )
        {
            cout << "(";
            for(set<string> :: iterator it2 = setLvTaxaIds.begin(); it2 != setLvTaxaIds.end(); ++it2)
            {
                if( it2 != setLvTaxaIds.begin() )
                {
                    cout << ",";
                }
                cout << *it2;
            }
            cout << ")\n";
        }
    }
}
void GenealogicalNetwork :: DumpAdmixNodeIds() const
{
    //
    // find all the admixture nodes
    vector<GenealogicalNetworkNode *> listMixNodes;
    GetMixNodes( listMixNodes );
    for(int i=0; i<(int)listMixNodes.size(); ++i)
    {
        cout << "Admixture population: ";
        set<int> setLvTaxaIds;
        listMixNodes[i]->GetAllLeafTaxonIdsUnder(setLvTaxaIds);
        if( setLvTaxaIds.size() == 1 )
        {
            cout << *( setLvTaxaIds.begin() ) << endl;
        }
        else if(setLvTaxaIds.size() > 1 )
        {
            cout << "(";
            for(set<int> :: iterator it2 = setLvTaxaIds.begin(); it2 != setLvTaxaIds.end(); ++it2)
            {
                if( it2 != setLvTaxaIds.begin() )
                {
                    cout << ",";
                }
                cout << *it2;
            }
            cout << ")\n";
        }
    }
}


// Read from file
// assume (i) no space in the tree (the first row),
// (ii) each of the following rows specifies a migration event
void GenealogicalNetwork :: ReadFromFileTreeMixFormat( const string &strFileTreemix, TaxaMapper *pMap )
{
//cout << "ReadFromFileTreeMixFormat: taxaMap: ";
//pMap->Dump();
    
    // first read from file
    ifstream inFile(strFileTreemix);
    if( inFile.is_open() == false )
    {
        throw string("ReadFromFileTreeMixFormat: input file: does not exist");
        exit(1);
    }
    // read in the tree
    vector< pair<set<int>, set<int> > > listMigInfo;
    vector< pair<double, double > > listMigPosInfo;
    string strLine;
    int numLines = 0;
    while( std::getline(inFile, strLine) )
    {
        //
        if( numLines == 0 )
        {
cout << "Reading tree from: " << strLine << endl;
            // read in the tree
            //remove_if(strLine.begin(), strLine.end(), isspace);
            std::string::iterator end_pos = std::remove(strLine.begin(), strLine.end(), ' ');
            strLine.erase(end_pos, strLine.end());
            int numLeaves = GetNewickNumLeaves(strLine);
cout << "Num of leaves: " << numLeaves << endl;
            PhylogenyTreeBasic phTree;
            //TaxaMapper taxMap;
            phTree.ConsOnNewick( strLine, numLeaves, false, pMap );
//cout << "Tree: ";
            string strCons;
            phTree.ConsNewick(strCons, false, 1.0, true);
cout << "Constructed tree: " << strCons << endl;
//phTree.Dump();
            
            // init net with this tree newick
            InitWithTree( strCons, numLeaves );
        }
        else
        {
            // now read the migration events
cout << "Read one migration: " << strLine << endl;
            std::istringstream iss(strLine);
            string strJunk;
            iss >> strJunk >> strJunk >> strJunk >> strJunk;
            string msrc, mdest;
            iss >> msrc >> mdest;
            double plen1 = 0.0, plen2 = 0.0;
            std::string::size_type pos1 = msrc.find_last_of(':');
            std::string::size_type pos2 = msrc.find_last_of(')');
            if (pos1 != std::string::npos && (pos2 == std::string::npos || pos1 > pos2))
            {
                // get mig pos info
                string strps = msrc.substr(pos1+1);
                plen1 = std::stod(strps);
                
                msrc = msrc.substr(0, pos1);

            }
            pos1 = mdest.find_last_of(':');
            pos2 = msrc.find_last_of(')');
            if (pos1 != std::string::npos && (pos2 == std::string::npos || pos1 > pos2))
            {
                // get mig pos info
                string strps = mdest.substr(pos1+1);
                plen2 = std::stod(strps);
                
                mdest = mdest.substr(0, pos1);
            }
            cout << "msrc: " << msrc << ", mdest: " << mdest << ", plen1:" << plen1 << ", plen2:" << plen2 << endl;
            
            set<string> setStrs1;
            GetTaxaFromNWStr(msrc, setStrs1);
            set<int> ss1;
            for(auto x: setStrs1)
            {
cout << "src: " << x << ", mapped id: " << pMap->GetId(x) << endl;
                ss1.insert( pMap->GetId(x) );
            }
cout << "ss1: ";
DumpIntSet(ss1);
            set<string> setStrs2;
            GetTaxaFromNWStr(mdest, setStrs2);
            set<int> ss2;
            for(auto x: setStrs2)
            {
cout << "dest: " << x << ", mapped id: " << pMap->GetId(x) << endl;
                ss2.insert( pMap->GetId(x) );
            }
cout << "ss2: ";
DumpIntSet(ss2);
            listMigInfo.push_back( std::make_pair(ss1, ss2)  );
            listMigPosInfo.push_back(std::make_pair(plen1, plen2));
        }
        ++numLines;
    }
    
    // now adding migration edges
    for(unsigned int i=0; i<listMigInfo.size(); ++i)
    {
        // first collect subtree info
        map< set<int>, map<double, GenealogicalNetworkBranch *> > mapSetLinsST;
        set<GenealogicalNetworkBranch *> setAllBranches;
        GetAllBranches( setAllBranches );
        for(auto x: setAllBranches)
        {
            double brl = x->GetLength();
            GenealogicalNetworkNode *pnt = x->GetDestNode();
            set<int> ss;
            pnt->GetAllLeafTaxonIdsUnder(ss);
            cout << "Set of taxon ids for node: ";
            DumpIntSet(ss);
            mapSetLinsST[ss][brl] = x;
        }
        
        set<int> ssSrc = listMigInfo[i].first;
        set<int> ssDest = listMigInfo[i].second;
        double bSrc = listMigPosInfo[i].first;
        double bDest = listMigPosInfo[i].second;
        auto it1 = mapSetLinsST.find(ssSrc);
        auto it2 = mapSetLinsST.find(ssDest);
        if( it1 == mapSetLinsST.end() || it2 == mapSetLinsST.end() )
        {
            YW_ASSERT_INFO(false, "Double check the network file: wrong migration source or destination");
        }
        // find the branch to be srce/dest
        GenealogicalNetworkBranch *brSrcToUse = NULL;
        const double LEN_SMALL = 0.00000001;
        double lenSrc = LEN_SMALL;
        for(auto y : it1->second)
        {
            if( bSrc <= y.first)
            {
                brSrcToUse = y.second;
            }
            else
            {
                lenSrc = y.first;
            }
        }
        if( brSrcToUse == NULL )
        {
            // if not found, just use the first one
            brSrcToUse = it1->second.begin()->second;
        }
        GenealogicalNetworkBranch *brDestToUse = NULL;
        double lenDest = LEN_SMALL;
        for(auto y : it2->second)
        {
            if( bDest <= y.first)
            {
                brDestToUse = y.second;
            }
            else
            {
                lenDest = y.first;
            }
        }
        if( brDestToUse == NULL )
        {
            // if not found, just use the first one
            brDestToUse = it2->second.begin()->second;
        }
        //
        AddMixBtwBranches( brDestToUse, brSrcToUse );
        // update info on
        GenealogicalNetworkNode *pnSrcOrig = brSrcToUse->GetDestNode();
        GenealogicalNetworkNode *pnSrcOrigPar = pnSrcOrig->GetParentSingle();
        YW_ASSERT_INFO( pnSrcOrigPar != NULL, "Cannot be null" );
        GenealogicalNetworkBranch *pnSrcOrigParBr = pnSrcOrigPar->GetParentBrSingle();
        int round = 1;
        while( it1->second.find(lenSrc) != it1->second.end() )
        {
            lenSrc += LEN_SMALL/(++round);
        }
        it1->second[lenSrc] = pnSrcOrigParBr;
        
        GenealogicalNetworkNode *pnDestOrig = brDestToUse->GetDestNode();
        GenealogicalNetworkNode *pnDestOrigPar = pnDestOrig->GetParentSingle();
        YW_ASSERT_INFO( pnDestOrigPar != NULL, "Cannot be null" );
        GenealogicalNetworkBranch *pnDestOrigParBr = pnDestOrigPar->GetAnces1();
        if( pnDestOrigParBr->GetSrcNode() == pnSrcOrigPar )
        {
            // use the other branch
            pnDestOrigParBr = pnDestOrigPar->GetAnces2();
        }
        round = 1;
        while( it2->second.find(lenDest) != it2->second.end() )
        {
            lenDest += LEN_SMALL/(++round);
        }
        it2->second[lenDest] = pnDestOrigParBr;
    }
    
    // dump out
    cout << "Afte reading from file: network is: ";
    Dump();
    cout << "Contained trees: \n";
    this->DumpMargTrees(false);
}


//
void GenealogicalNetwork :: CreateMapTaxonIdToUserId( map<string,string> &mapTaxonIdToUserId ) const
{
    //
    for( set<GenealogicalNetworkNode *> :: const_iterator it = setNodes.begin(); it != setNodes.end(); ++it )
	{
        GenealogicalNetworkNode *pnodeCur = *it;
        char buf[100];
        snprintf(buf, 100, "%d", pnodeCur->GetTaxonId() );
        string bufstr = buf;
        string taxonStrUser = pnodeCur->GetTaxonStrUser();
        mapTaxonIdToUserId.insert( map<string,string>  :: value_type(bufstr, taxonStrUser) );
    }
}

void GenealogicalNetwork :: GetAllTaxa(set<int> &setTaxa) const
{
    vector<GenealogicalNetworkNode *> listLeaves;
    GetLeaves( listLeaves );
    setTaxa.clear();
    for(int i=0; i<(int)listLeaves.size(); ++i)
    {
        setTaxa.insert( listLeaves[i]->GetTaxonId() );
    }
}

void GenealogicalNetwork :: GetListofNodesTopdown( vector<GenealogicalNetworkNode *> &listNodes ) const
{
#if 0
    listNodes.clear();
    for(set<GenealogicalNetworkNode *> :: iterator it = setNodes.begin(); it != setNodes.end(); ++it )
    {
        GenealogicalNetworkNode *pn = *it;
        // place it in the first position (from the right) that it is not ancestral to
        vector<GenealogicalNetworkNode *> :: iterator itCur = listNodes.end();
        while( true )
        {
            if( itCur == listNodes.begin() )
            {
                break;
            }
            --itCur;
            GenealogicalNetworkNode *pnPre = *itCur;
            if( pnPre->IsAncestralTo( pn ) == true )
            {
                break;
            }
        }
        listNodes.insert( itCur, pn  );
    }
#endif
    listNodes.clear();
    PopulateVecBySetGen( listNodes, this->setNodes);
    // now sort
    //std::sort(listNodes.begin(), listNodes.end(), [](GenealogicalNetworkNode *pa, GenealogicalNetworkNode *pb) -> bool
    //{
    //    return pa->IsAncestralTo(pb);
    //});
    // sort
    for(unsigned int i=0; i<listNodes.size(); ++i)
    {
        for(unsigned int j=i+1; j<listNodes.size(); ++j)
        {
            if( listNodes[j]->IsAncestralTo(listNodes[i]) )
            {
                // swap
                GenealogicalNetworkNode *ptemp = listNodes[i];
                listNodes[i] = listNodes[j];
                listNodes[j] = ptemp;
            }
        }
    }
}

static bool GNCmpPreferSmall(const pair<set<int>, GenealogicalNetworkNode *> &s1, const pair<set<int>, GenealogicalNetworkNode *> &s2)
{
    if( s1.first.size() < s2.first.size() )
    {
        return true;
    }
    else if( s1.first.size() > s2.first.size() )
    {
        return false;
    }
    else
    {
        return s1 < s2;
    }
}

void GenealogicalNetwork :: GetListofNodesOrdered( vector<GenealogicalNetworkNode *> &listNodes ) const
{
    listNodes.clear();
    // order all leaf nodes by id
    map<GenealogicalNetworkNode *, set<int> > mapNodesIds;
    vector<GenealogicalNetworkNode *> listLeaves;
    GetLeaves( listLeaves );
    for(auto x : setNodes)
    {
        set<GenealogicalNetworkNode *> setLeaves;
        x->GetAllLeavesUnder(setLeaves);
        set<int> ss;
        for(auto y: setLeaves)
        {
            ss.insert( y->GetTaxonId() );
        }
        mapNodesIds[x] = ss;
    }
    // sort based on the set of ids
    vector<pair<set<int>, GenealogicalNetworkNode *> > listPPs;
    for(auto x : mapNodesIds)
    {
        listPPs.push_back( std::make_pair(x.second, x.first) );
    }
    std::sort(listPPs.begin(), listPPs.end(), GNCmpPreferSmall);
    for(unsigned int i = 0; i< listPPs.size(); ++i)
    {
        listNodes.push_back(listPPs[i].second);
    }
}

void GenealogicalNetwork :: FindAllPathsBtwAllPairTaxa(map<pair<int,int>, set<set<GenealogicalNetworkBranch *> > > &mapPathsBtwPairTaxa) const
{
    // approach: for each node v, keep track for all nodes v1 that is above (ancestral to, i.e. backwwards in time) v the set of paths from v to v1; then for each node u and two leaves v1 and v2, there is a LCA path between v1 and v2 (i.e., v1 to LCA(v1,v2) plus v2 to LCA(v1,v2)) with the condition there is no common edges between these two paths
    mapPathsBtwPairTaxa.clear();
    // [v, [v1, <set of paths from v up to v1> ] ]: v1 is ancestral to v
    map<GenealogicalNetworkNode *, map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > > mapNodeToAncesPaths;
    // init the root
    //set<set<GenealogicalNetworkBranch *> > setTriv;
    //set<GenealogicalNetworkBranch *> setTriv2;
    //setTriv.insert(setTriv2);
    //mapNodeToAncesPaths[GetRoot()][GetRoot()] = setTriv;
    // start from the root and going downwards
    vector<GenealogicalNetworkNode *> listNodes;
    GetListofNodesTopdown( listNodes );
#if 0
cout << "List nodes top down: ";
for(auto x : listNodes)
{
x->Dump();
}
#endif
    YW_ASSERT_INFO(listNodes[0] == GetRoot(), "The first node must be the root");
    for(unsigned int i=0; i<listNodes.size(); ++i)
    {
        //
        GenealogicalNetworkNode *pnCurr = listNodes[i];
#if 0
cout << "pnCurr: ";
pnCurr->Dump();
cout << endl;
#endif
        
        // To begin with, insert a path of no dge from pnCurr to self
        set<GenealogicalNetworkBranch *> setBrs1;
        mapNodeToAncesPaths[pnCurr][pnCurr].insert(setBrs1);
        
        // get its two ancestors
        vector<GenealogicalNetworkNode *> listAncesCurr;
        vector<GenealogicalNetworkBranch *> listAncesCurrBrs;
        GenealogicalNetworkBranch *pBrAnc1 = pnCurr->GetAnces1();
        if( pBrAnc1 != NULL )
        {
            listAncesCurrBrs.push_back(pBrAnc1);
            GenealogicalNetworkNode *pAnc1 = pBrAnc1->GetSrcNode();
            listAncesCurr.push_back(pAnc1);
        }
        GenealogicalNetworkBranch *pBrAnc2 = pnCurr->GetAnces2();
        if( pBrAnc2 != NULL )
        {
            listAncesCurrBrs.push_back( pBrAnc2 );
            GenealogicalNetworkNode *pAnc2 = pBrAnc2->GetSrcNode();
            listAncesCurr.push_back(pAnc2);
        }
        for(unsigned int jj=0; jj<listAncesCurr.size(); ++jj)
        {
            GenealogicalNetworkNode *pna = listAncesCurr[jj];
#if 0
cout << "pna: ";
pna->Dump();
cout << endl;
#endif
            GenealogicalNetworkBranch *pnbr = listAncesCurrBrs[jj];
            map<GenealogicalNetworkNode *, map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > > :: iterator itg = mapNodeToAncesPaths.find(pna);
            YW_ASSERT_INFO( itg != mapNodeToAncesPaths.end(), "Fail to find 12345");
            for( map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > :: iterator itg2 = itg->second.begin(); itg2 != itg->second.end(); ++itg2 )
            {
                GenealogicalNetworkNode *pnAnc = itg2->first;
                for( set<set<GenealogicalNetworkBranch *> > :: iterator itg3 = itg2->second.begin(); itg3 != itg2->second.end(); ++itg3 )
                {
                    set<GenealogicalNetworkBranch *> setBrs = *itg3;
                    // add the branch from
                    setBrs.insert( pnbr );
                    // add to the record
                    mapNodeToAncesPaths[pnCurr][pnAnc].insert( setBrs );
                }
            }
        }
    }
    // get all the leaves
    vector<GenealogicalNetworkNode *> listLeaves;
    GetLeaves(listLeaves);
    // consider all pairs of leaves
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        GenealogicalNetworkNode *pn1 = listLeaves[i];
        map<GenealogicalNetworkNode *, map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > > :: iterator it1 = mapNodeToAncesPaths.find(pn1);
        YW_ASSERT_INFO( it1 != mapNodeToAncesPaths.end(), "Fail to find 345");
        set<GenealogicalNetworkNode *> setNodeAnc1;
        for( map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > :: iterator it11 = it1->second.begin(); it11 != it1->second.end(); ++it11 )
        {
            setNodeAnc1.insert(it11->first);
        }
        for(unsigned int j=i+1; j<listLeaves.size(); ++j)
        {
            GenealogicalNetworkNode *pn2 = listLeaves[j];
            map<GenealogicalNetworkNode *, map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > > :: iterator it2 = mapNodeToAncesPaths.find(pn2);
            YW_ASSERT_INFO( it2 != mapNodeToAncesPaths.end(), "Fail to find 346");
            
            set<GenealogicalNetworkNode *> setNodeAnc2;
            for( map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > :: iterator it21 = it2->second.begin(); it21 != it2->second.end(); ++it21 )
            {
                setNodeAnc2.insert(it21->first);
            }
            
            // find the common node
            set<GenealogicalNetworkNode *> setJoin;
            JoinSetsGen(setNodeAnc1, setNodeAnc2, setJoin);
            for(set<GenealogicalNetworkNode *> :: iterator it33 = setJoin.begin(); it33 != setJoin.end(); ++it33)
            {
                GenealogicalNetworkNode *pna2 = *it33;
                map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > :: iterator it41 = it1->second.find(pna2);
                YW_ASSERT_INFO( it41 != it1->second.end(), "Fail to find 51");
                map<GenealogicalNetworkNode *, set<set<GenealogicalNetworkBranch *> > > :: iterator it42 = it2->second.find(pna2);
                YW_ASSERT_INFO( it42 != it2->second.end(), "Fail to find 51");
                
                for(set<set<GenealogicalNetworkBranch *> > :: iterator it51 = it41->second.begin(); it51 != it41->second.end(); ++it51 )
                {
                    for(set<set<GenealogicalNetworkBranch *> > :: iterator it52 = it42->second.begin(); it52 != it42->second.end(); ++it52 )
                    {
                        // ensure there is no overlap
                        set<GenealogicalNetworkBranch *> setJoin2;
                        JoinSetsGen(*it51, *it52, setJoin2 );
                        if( setJoin2.size() == 0 )
                        {
                            // found a valid LCA path
                            set<GenealogicalNetworkBranch *> setBrCombo = *it51;
                            UnionSetsGen( setBrCombo, *it52 );
                            int pid1 = pn1->GetTaxonId();
                            int pid2 = pn2->GetTaxonId();
                            pair<int,int> pp(pid1, pid2);
                            OrderInt( pp.first, pp.second );
                            mapPathsBtwPairTaxa[pp].insert( setBrCombo );
                        }
                    }
                }
            }
        }
    }
}

bool GenealogicalNetwork :: IsOutgroup(int taxon) const
{
    // is taxon outgroup?
    GenealogicalNetworkNode *pn = GetLeafWithTaxon(taxon);
    YW_ASSERT_INFO(pn != NULL, "Fail to find the leaf with the taxon");
    GenealogicalNetworkBranch *pb1 = pn->GetAnces1();
    GenealogicalNetworkBranch *pb2 = pn->GetAnces2();
    if( pb2 != NULL )
    {
        return false;
    }
    YW_ASSERT_INFO(pb1 != NULL, "cannot be null");
    return pb1->GetSrcNode()->IsRoot();
}

void GenealogicalNetwork :: GetAncesBranchesFrom(GenealogicalNetworkNode *pNode, set<GenealogicalNetworkBranch *> &setAncBrs) const
{
//return;
    // assume: network has no cycle
    // recursively finds; don't clear passed anc brs set
    if( pNode == NULL )
    {
        return;
    }
    GenealogicalNetworkBranch *pAncBr1 = pNode->GetAnces1();
    if( pAncBr1 != NULL )
    {
        setAncBrs.insert(pAncBr1);
        GetAncesBranchesFrom( pAncBr1->GetSrcNode(), setAncBrs );
    }
    GenealogicalNetworkBranch *pAncBr2 = pNode->GetAnces2();
    if( pAncBr2 != NULL )
    {
        setAncBrs.insert(pAncBr2);
        GetAncesBranchesFrom( pAncBr2->GetSrcNode(), setAncBrs );
    }
}

void GenealogicalNetwork :: PruneNonexistBrs(set<GenealogicalNetworkBranch *> &setAncBrs) const
{
//return;
//cout << "PruneNonexistBrs: input setAncBrs size: " << setAncBrs.size() << endl;
    for(auto it = setAncBrs.begin(); it != setAncBrs.end();)
    {
        if(setBranches.find(*it) == setBranches.end() || setNodes.find((*it)->GetSrcNode() ) == setNodes.end() || setNodes.find((*it)->GetDestNode() ) == setNodes.end()  ) {
            setAncBrs.erase(it++);
        }
        else {
            ++it;
        }
    }
//cout << "After pruning: size is: " << setAncBrs.size() << endl;
}

void GenealogicalNetwork :: GetTopologySignature(std::map<string,int> &margTreeCnts) const
{
    margTreeCnts.clear();
    // get all marg trees
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTreesWithFreq;
    RetriveAllMarginalTrees( mapMargTreesWithFreq );
    for(auto it : mapMargTreesWithFreq)
    {
        string strNW = it.second.GetMargTree()->GetNewickSorted(false);
        ++margTreeCnts[strNW];
        
        // free
        delete it.second.GetMargTree();
    }
}
