//
//  AGGeneTreeProb.cpp
//  
//
//  Created by Yufeng Wu on 5/15/24.
//

#include "AGGeneTreeProb.hpp"
#include "GenealogicalNetwork.h"
#include "PhylogenyTreeBasic.h"
#include "GeneSpeciesTreeProb.h"
#include "ApproxGeneTreeProb2.h"
#include "AGFastLenOpt.hpp"
#include "Utils3.h"
#include <cmath>
#include <thread>
#include <mutex>
#include "UtilsNumerical.h"
#include "GenealogicalNetworkProbStore.h"


//***********************************************************************************
// storing computed prob of networks

AGGeneTreeProbDepot & AGGeneTreeProbDepot :: Instance()
{
    static AGGeneTreeProbDepot inst;
    return inst;
}

// query whether this tree has been computed or not
bool AGGeneTreeProbDepot :: QueryProbForNetTree(const std::string &netTreeNW, std::vector<double> &listGTProbs) const
{
return false;
    AGGeneTreeProbDepot *pthis = const_cast<AGGeneTreeProbDepot *>(this);
    ++pthis->numTotQuery;
    listGTProbs.clear();
    auto it = mapNetTreeGTProbs.find(netTreeNW);
    if( it == mapNetTreeGTProbs.end() )
    {
        return false;
    }
    ++pthis->numQueryCached;
    listGTProbs = it->second;
    YW_ASSERT_INFO( listGTProbs.size() > 0, "Cannot be empty01" );
    return true;
}

// store a newly computed tree
void AGGeneTreeProbDepot :: AppendProb( const std::string &netTreeNW, const std::vector<double> &listGTProbs )
{
return;
    ++numStores;
    YW_ASSERT_INFO(listGTProbs.size() >0, "Cannot be empty00" );
    auto it = mapNetTreeGTProbs.find(netTreeNW);
    if( it == mapNetTreeGTProbs.end() )
    {
//cout << "AppendProb: tree:" << netTreeNW << ", probs=";
//DumpDoubleVec(listGTProbs);
        // update
        mapNetTreeGTProbs[netTreeNW] = listGTProbs;
    }
}

void AGGeneTreeProbDepot :: DumpStats() const
{
    cout << "Cache probability computation: [" <<numQueryCached << "] cache hits among total " << numTotQuery << " queries," << " num of storing operations: " << numStores << endl;
#if 0
    cout << "List of species trees: \n";
    for(auto x: mapNetTreeGTProbs)
    {
        cout << x.first << endl;
    }
#endif
}

// clear out
void AGGeneTreeProbDepot :: Clear()
{
    mapNetTreeGTProbs.clear();
}

//***********************************************************************************
// Calc prob of genealogical networks: enhanced

AGGeneTreeProb2 :: AGGeneTreeProb2(GenealogicalNetwork &netIn, vector<PhylogenyTreeBasic *> &listGeneTreesIn, int tidIn) : network(netIn), listGeneTrees(listGeneTreesIn), fFresh(false), fWeightUnchanged(false), probLastComputed(MAX_NEG_DOUBLE_VAL), threadIndex(tidIn)
{
    //
//cout << "AGGeneTreeProb2: init: net: ";
//network.Dump();
    Init();
//cout << "Init done.\n";
}

AGGeneTreeProb2 :: ~AGGeneTreeProb2()
{
    Clear();
}

void AGGeneTreeProb2 :: Init()
{
    Clear();
    
    // init the prob computing
    set<GenealogicalNetworkMTreeCode > setMTCodes;
    network.RetriveAllNetCode( setMTCodes );
    int indexNext = 0;
    for( set<GenealogicalNetworkMTreeCode > :: iterator it = setMTCodes.begin(); it != setMTCodes.end(); ++it )
    {
        //
        mapNetCodeIndex.insert( map<GenealogicalNetworkMTreeCode, int> :: value_type( *it, indexNext++ ) );
    }
//cout << "Num of network codes: " << mapNetCodeIndex.size() << endl;
    
    // create the prob computing tools
    CreateProbCalc();
    
    // init cache to garbage
    listCachedProbs.resize( listGeneTrees.size() );
    for(int i=0; i<(int)listCachedProbs.size(); ++i)
    {
        for(int j=0; j<(int)listMargTreesInNet.size(); ++j )
        {
            listCachedProbs[i].push_back(0.0);
        }
    }
    
//cout << "now init marginal trees...\n";
    // initialize the marginal tree branch vs net branch
    for( set<GenealogicalNetworkMTreeCode > :: iterator it = setMTCodes.begin(); it != setMTCodes.end(); ++it )
    {
        auto it2 = mapNodeToMargTree.find(*it);
        YW_ASSERT_INFO(it2 != mapNodeToMargTree.end(), "Wrong222");
        MarginalTree *mt = it2->second;
        
        network.MapNetCodedEdgeToMTBrs(*it, *mt, mapMargTreeNetBr[*it] );
    }
//cout << "AGGeneTreeProb2: init done.\n";
}

double AGGeneTreeProb2 :: GetCachedProbValFor(int gtIndex, int mtIndex)
{
    return listCachedProbs[gtIndex][mtIndex];
}
void AGGeneTreeProb2 :: SetCachedProbValFor(int gtIndex, int mtIndex, double val)
{
    //
    listCachedProbs[gtIndex][mtIndex] = val;
}

void AGGeneTreeProb2 :: Clear()
{
    for(int i=0; i<(int)listGTProbPtrs.size(); ++i )
    {
        delete listGTProbPtrs[i];
        listGTProbPtrs[i] = NULL;
    }
    for(int i=0; i<(int)listMargTreesInNet.size(); ++i)
    {
        delete listMargTreesInNet[i];
    }
    listGTProbPtrs.clear();
    listMargTreesInNet.clear();
    listGenetreeProbWts.clear();
    listMargTreeFresh.clear();
}

double AGGeneTreeProb2 :: CalcProb()
{
//cout << "AGGeneTreeProb2: CalcProb\n";
//#if 0
    if( IsFresh() )
    {
//cout << "Fresh: so compute cached\n";
        return CalcProbFromCached();
    }
    fWeightUnchanged = true;
    SetFresh(true);
//#endif
    
//cout << "here\n";
    //Init();
    
#if 0
    if( numThreads > 1 )
    {
        probLastComputed = CalcProbMultithread();
        return probLastComputed;
    }
#endif
    
    //
    //if( IsInit() == false )
    //{
    //    CreateProbCalc();
    //}
    //map<string, double> mapTreeTopoProb;
    
    YW_ASSERT_INFO(listGTProbPtrs.size() > 0, "Must have at least one tree");
    // get all the probs
    vector<vector<double> > listAllProbs(listGeneTrees.size());
    for(unsigned int i=0; i<listGeneTrees.size(); ++i)
    {
        listAllProbs[i].resize( listMargTreesInNet.size() );
    }
    
    // snap marg tree br lengths
    SnapMargTreeLens();
    //static int numTot = 0;
    //static int numCacheHits = 0;
    
    for(int j=0; j<(int)listGTProbPtrs.size(); ++j)
    {
        vector<double> listLogProbs;
        
        // check if it has been computed before
        string netTreeNW = listGTProbPtrs[j]->GetSpeciesTreeNW();
//cout << j << "-th netTreeNW: " << netTreeNW << endl;
        bool fPre = AGGeneTreeProbDepot::Instance().QueryProbForNetTree(netTreeNW, listLogProbs);
        //++numTot;
        if( fPre == false )
        {
            listGTProbPtrs[j]->CalcLogProb(listLogProbs);
//cout << "CalcProb: logprobs: ";
//DumpDoubleVec(listLogProbs);
            // save it
            AGGeneTreeProbDepot::Instance().AppendProb(netTreeNW, listLogProbs);
        }
        else
        {
            //++numCacheHits;
            //cout << "++Cache hit: " << numCacheHits << " out of " << numTot << endl;
        }
        
        for(unsigned int i=0; i<listLogProbs.size(); ++i)
        {
            YW_ASSERT_INFO( i<listAllProbs.size() && j<(int)listAllProbs[i].size(), "Overflow222" );
            listAllProbs[i][j] = listLogProbs[i];
            
            SetCachedProbValFor(i, j, listLogProbs[i]);
        }
    }
    
    // now calculate prob
    double res = 0.0;
    for(int i=0; i<(int)listAllProbs.size(); ++i)
    {
        // multiply by the weight for each
        for(int j=0; j<(int)listAllProbs[i].size(); ++j)
        {
            YW_ASSERT_INFO( j<listGenetreeProbWts.size(), "Overflow334" );
            double wtStep = log(listGenetreeProbWts[j]);
            listAllProbs[i][j] += wtStep;
        }
        
        double logprobstep = GetLogSumOfLogs( listAllProbs[i] );
        res += logprobstep;
//string strNW;
//listGeneTrees[i]->ConsNewick(strNW);
//cout << "Gene tree: " << strNW << ", prob: " << logprobstep << " components: ";
//DumpDoubleVec(listAllProbs[i]);
    }

//cout << "GenealogicalNetworkProb :: CalcProb: prob = " << res << endl;
    // after one round of computation: all marg trees become fresh
    //SetAllMargTreeFresh();
    probLastComputed = res;
    return res;
}

double AGGeneTreeProb2 :: CalcProbFromCached()
{
    if( fWeightUnchanged == true )
    {
        return probLastComputed;
    }
    fWeightUnchanged = true;
    // re-calculate since weight changed
    double res = 0.0;
    for(int i=0; i<(int)listGeneTrees.size(); ++i)
    {
        // get the tree topology to see if this topology has been computed before
        vector<double> listLogProbs;
        //double logprobMax = MAX_NEG_DOUBLE_VAL;
        for(int j=0; j<(int)listGenetreeProbWts.size(); ++j)
        {
            // is this fresh? if so, just use cached value
            double logponetree = GetCachedProbValFor(i, j);
            double wtStep = listGenetreeProbWts[j];
            logponetree += log(wtStep);
            //if(logprobMax < logponetree )
            //{
            //    logprobMax = logponetree;
            //}
            
            listLogProbs.push_back(logponetree);
            //cout << "Marginal tree: wt=" << wtStep << ", logponetree: " << logponetree << endl;
        }
        double logprobstep = GetLogSumOfLogs( listLogProbs );
        res += logprobstep;
    }
    probLastComputed = res;
    return res;
}

void AGGeneTreeProb2 :: CreateProbCalc()
{
    // retrieve marginal trees
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTreesWithFreq ;
    network.RetriveAllMarginalTrees( mapMargTreesWithFreq );
//cout << "Number of marginal trees: " << mapMargTreesWithFreq.size() << endl;
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: iterator it = mapMargTreesWithFreq.begin(); it != mapMargTreesWithFreq.end(); ++it )
    {
        listMargTreesInNet.push_back( it->second.GetMargTree() );
        listGenetreeProbWts.push_back( it->second.GetFreq() );
        mapNodeToMargTree[it->first] = it->second.GetMargTree();
    }
    
    for(unsigned int i=0; i<listMargTreesInNet.size(); ++i)
    {
//cout << "Creating prob calculator...i:" << i << endl;
//cout << "marginal tree: " << listMargTreesInNet[i]->GetNewickSorted(false) << endl;
        //listGTProbPtrs.push_back( new STApproxGeneTreeProb( *listMargTreesInNet[i], listGeneTrees ) );
        //listGTProbPtrs.push_back( new STApproxGeneTreeProb2( *listMargTreesInNet[i], listGeneTrees ) );
        listGTProbPtrs.push_back( new ApproxGeneTreeProbHeu2( *listMargTreesInNet[i], listGeneTrees, threadIndex ) );
#if 0
        SGSTApproxGeneTreeProb *pProb = new SGSTApproxGeneTreeProb( *listMargTreesInNet[i], listGeneTrees );
        //pProb->SnapSpeciesTreeLens();
        listGTProbPtrs.push_back(pProb);
#endif
    }
//cout << "Done: prob calculator\n";
}


void AGGeneTreeProb2 :: UpdateBranch( GenealogicalNetworkBranch *pBrChange, double lenDiff )
{
    fWeightUnchanged = false;
    
#if 0
cout << "UpdateBranch: to lenDiff= " << lenDiff << ", ";
if( pBrChange != NULL)
{
pBrChange->Dump();
}
else
{
cout << "  root branch";
}
cout << endl;
#endif
    
    // for each gene tree contained, if it is affected by this branch length change, update it
    for(auto x : mapMargTreeNetBr)
    {
        auto it = x.second.find(pBrChange);
        if( it == x.second.end() )
        {
            continue;
        }
        //YW_ASSERT_INFO(it != x.second.end(), "Fail to find the branch record");
        int mtBrInd = it->second;
        // TBD YW: don't allow change branch length of root branch
        if( mtBrInd >= 0 && mtBrInd < 2*GetNumTaxa()-2 )
        {
            auto it4 = mapNetCodeIndex.find(x.first);
            YW_ASSERT_INFO( it4 != mapNetCodeIndex.end(), "Fail to find555" );
            int tid = it4->second;
            YW_ASSERT_INFO(tid >=0 && tid<(int)listMargTreesInNet.size(), "Overflow336");
            // now update the marginal tree branch length by the difference
            auto it2 = mapNodeToMargTree.find(x.first);
            YW_ASSERT_INFO(it2 != mapNodeToMargTree.end(), "Fail to find the marginal tree");
            MarginalTree *pmt = it2->second;
            YW_ASSERT_INFO(pmt != NULL && mtBrInd >=0 && mtBrInd < pmt->GetTotNodesNum()-1, "Overflow337");
            double brLenOrig = pmt->GetEdgeLen(mtBrInd);
            double brLenNew = brLenOrig + lenDiff ;
            
            //listMargTreesInNet[tid]->SetBranchLen(mtBrInd, brLenNew);
            
            //double lenOrig = listMargTreesInNet[tid]->GetEdgeLen(mtBrInd);
            double lenConv = AGGTBranchLengthGrid::Instance().SnapToGrid(brLenNew);
            // 06/13/24: no longer snap branch lengths
            //double lenConv = brLenNew;
            listMargTreesInNet[tid]->SetBranchLen(mtBrInd, lenConv);
            
            // update branch length
            listGTProbPtrs[tid]->UpdateSTBrLen(mtBrInd, lenConv);
            
            //listMargTreesInNet[tid]->SetBranchLen(mtBrInd, brLenNew);
            //listGTProbPtrs[tid]->SnapSpeciesTreeLenAt(mtBrInd);
            
            //pmt->SetBranchLen( mtBrInd, brLenNew);
            SetFresh(false);
//cout << "set marg tree branch: " << mtBrInd << " to length: " << brLenNew << ", from length: " << brLenOrig << ", tree becomes: " << pmt->GetNewickSorted(true) << " tr address:" << (long)pmt << endl;
        }
    }
    
    // update each prob-calculator too
    for(unsigned int i=0; i<listGTProbPtrs.size(); ++i)
    {
        listGTProbPtrs[i]->OnUpdateSTBrLen();
    }
    
    //SetFresh(false);
    return;

#if 0
    // for each gene tree contained, if it is affected by this branch length change, update it
    for(auto x : mapMargTreeNetBr)
    {
        auto it = x.second.find(pBrChange);
        if( it == x.second.end() )
        {
            continue;
        }
        //YW_ASSERT_INFO(it != x.second.end(), "Fail to find the branch record");
        int mtBrInd = it->second;
        if( mtBrInd >= 0 )
        {
            // now update the marginal tree branch length by the difference
            auto it2 = mapNodeToMargTree.find(x.first);
            YW_ASSERT_INFO(it2 != mapNodeToMargTree.end(), "Fail to find the marginal tree");
            MarginalTree *pmt = it2->second;
            double brLenOrig = pmt->GetEdgeLen(mtBrInd);
            double brLenNew = brLenOrig + lenDiff ;
            pmt->SetBranchLen( mtBrInd, brLenNew);
            
            auto it3 = mapNetCodeIndex.find(x.first);
            YW_ASSERT_INFO(it3 != mapNetCodeIndex.end(), "Fail to find marginal tree index");
            int mtIndex = it3->second;
            // update branch length
            for(unsigned int ii=0; ii<listGenetreeProbPtrs.size(); ++ii)
            {
                // set branch length
                listGenetreeProbPtrs[ii][mtIndex]->SetSTBranchLen(mtBrInd, brLenNew);
            }
            
            vector<double> vecLensOrig;
            SnapMargTreeLensFor(0, mtIndex, vecLensOrig);
            
            // check to see if this has already been computed before
            vector<double> listLogProbsStore;
            string netTreeNW = listGenetreeProbPtrs[0][mtIndex]->GetSpeciesTreeNW();
            //cout << "netTreeNW: " << netTreeNW << endl;
            RestoreBrLenFor(0, mtIndex, vecLensOrig);
            bool fPre = AGGeneTreeProbDepot::Instance().QueryProbForNetTree(netTreeNW, listLogProbsStore);
            //bool fPre = false;
            if( fPre == true )
            {
                // no need to update
                //pmt->SetBranchLen( mtBrInd, brLenOrig);
                // set prob of trees for all pre-computed trees
                for(unsigned int ii=0; ii<listLogProbsStore.size(); ++ii)
                {
                    // save it
                    SetCachedProbValFor(ii, mtIndex, listLogProbsStore[ii]);
                }
                
                continue;
            }
            

                for(unsigned int tr = 0; tr<listGenetreeProbPtrs.size(); ++tr)
                {
                    YW_ASSERT_INFO(mtIndex < listGenetreeProbPtrs[tr].size(), "Overflow 222" );
                    // re-claculate by changing the
                    map<int,set< LineageConfig,LTLinCfg > > origLinCfgs;
                    int threadID = 0;
                    double prNew = listGenetreeProbPtrs[tr][mtIndex]->TestNewBranchLenAt(threadID, mtBrInd, brLenNew, origLinCfgs, true);
                    //cout << "Updated branch: log-prob: " << prNew << ", tree:" << tr << ", marginal tree index:" << mtIndex << endl;
                    // save it
                    SetCachedProbValFor(tr, mtIndex, prNew);
                    
                    listLogProbsStore.push_back(prNew);
                }

            
            // save it
            AGGeneTreeProbDepot::Instance().AppendProb(netTreeNW, listLogProbsStore);
        }
    }
//cout << "Done with UpdateBranch\n";
    
//    SetFresh(false);
#endif
}

void AGGeneTreeProb2 :: UpdateWts()
{
    //
    listGenetreeProbWts.clear();
    network.GetWtsForAllMargTrees(listGenetreeProbWts);
    // set weight changed
    fWeightUnchanged = false;
    // set fresh so not to recompute
    SetFresh(true);
    //SetFresh(false);
}

int AGGeneTreeProb2 :: GetIndexForNetcode( const GenealogicalNetworkMTreeCode &ncode )
{
    //
    YW_ASSERT_INFO( mapNetCodeIndex.find(ncode) != mapNetCodeIndex.end(), "Fail to find1" );
    return mapNetCodeIndex[ncode];
}

void AGGeneTreeProb2 :: SnapMargTreeLens()
{
return;
    // snap all marg tree branch lengths to standard ones
    int numEdges = GetNumMargTreeEdges();
    //for(unsigned int i=0; i<listGenetreeProbPtrs.size(); ++i)
    //{
    for(unsigned int j=0; j<listMargTreesInNet.size(); ++j)
    {
        //
        for(int k=0; k<numEdges; ++k)
        {
            double lenOrig = listMargTreesInNet[j]->GetEdgeLen(k);
            double lenConv = AGGTBranchLengthGrid::Instance().SnapToGrid(lenOrig);
            listMargTreesInNet[j]->SetBranchLen(k, lenConv);
        }
    }
    //}
}

void AGGeneTreeProb2 :: SnapMargTreeLensFor(int gtIndex, int mtIndex, std::vector<double> &listOrigLens)
{
return;
    //
    // snap all marg tree branch lengths to standard ones
    int numEdges = GetNumMargTreeEdges();
    listOrigLens.resize(numEdges);
    //
    for(int k=0; k<numEdges; ++k)
    {
        double lenOrig = listMargTreesInNet[mtIndex]->GetEdgeLen(k);
        listOrigLens[k] = lenOrig;
        double lenConv = AGGTBranchLengthGrid::Instance().SnapToGrid(lenOrig);
        listMargTreesInNet[mtIndex]->SetBranchLen(k, lenConv);
    }
}
void AGGeneTreeProb2 :: RestoreBrLenFor(int gtIndex, int mtIndex, const std::vector<double> &listOrigLens)
{
return;
    //
    int numEdges = GetNumMargTreeEdges();
    //
    for(int k=0; k<numEdges; ++k)
    {
        double lenOrig = listOrigLens[k];
        listMargTreesInNet[mtIndex]->SetBranchLen(k, lenOrig);
    }
}

int AGGeneTreeProb2 :: GetNumMargTreeEdges() const
{
    YW_ASSERT_INFO(listMargTreesInNet.size()>0, "Wrong: not inited");
    return listMargTreesInNet[0]->GetTotNodesNum()-1;
}

int AGGeneTreeProb2 :: GetNumTaxa() const
{
    //
    return network.GetNumLeaves();
}
