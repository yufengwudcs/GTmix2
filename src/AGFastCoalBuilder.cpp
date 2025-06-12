//
//  AGFastCoalBuilder.cpp
//  
//
//  Created by Yufeng Wu on 3/18/24.
//

#include "AGFastCoalBuilder.hpp"
#include "AGFastLenOpt.hpp"
#include "PhylogenyTreeBasic.h"
#include "MarginalTree.h"
#include "Utils3.h"
#include "Utils4.h"
#include "DeepCoalescence.h"
#include <cmath>
#include <string>
#include <thread>
#include <mutex>
#include "AGGeneTreeProb.hpp"
#include "ApproxGeneTreeProb2.h"
#include "AGHeuristicSearch.hpp"

//const double MIN_INC_RATIO = 1.00001;
const double MIN_INC_RATIO = 1.001;
static bool fFastSearch = false;

extern int GetNumThreads();

extern void GetLeavesSpeciesFromGT( PhylogenyTreeBasic &treeGene, set<int> &species );

#if 0
// a useful function
void GetLeavesSpeciesFromGT( PhylogenyTreeBasic &treeGene, set<int> &species )
{
        // get all the species/labels in gene tree
        //
        vector<string> listLeafLabels;
        treeGene.GetAllLeafLabeles(listLeafLabels);
        // fill in the set
        for(int i=0; i<(int)listLeafLabels.size(); ++i)
        {
                // convert to integer
                int id = -1;
                sscanf( listLeafLabels[i].c_str(), "%d", &id );
                species.insert(id);
        }
}
#endif

//////////////////////////////////////////////////////////////////////////////////
// Store processed network topologies


AGProcessedNetsDepot& AGProcessedNetsDepot :: Instance()
{
    static AGProcessedNetsDepot inst;
    return inst;
}

bool AGProcessedNetsDepot :: IsNetProcessed( GenealogicalNetwork *pNet)
{
return false;
    map<int,int> sigNet;
    ConsNetSignature(pNet, sigNet);
    bool f = false;
    mut.lock();
    f = mapNettoId.find(sigNet) != mapNettoId.end();
    mut.unlock();
    return f;
}

void AGProcessedNetsDepot :: MarkNetProcessed(GenealogicalNetwork *pNet)
{
return;
    map<int,int> sigNet;
    ConsNetSignature(pNet, sigNet);
    mut.lock();
    mapNettoId.insert(sigNet);
    mut.unlock();
}

void AGProcessedNetsDepot :: SkipNet(GenealogicalNetwork *pNet)
{
    ++numSkipped;
}

void AGProcessedNetsDepot :: DumpStats() const
{
    cout << "AGProcessedNetsDepot: total number of novel networks processed: " << mapNettoId.size() << ", number of skipped networks: " << numSkipped << endl;
}

AGProcessedNetsDepot :: AGProcessedNetsDepot(): numSkipped(0)
{
}

void AGProcessedNetsDepot :: ConsNetSignature(GenealogicalNetwork *pNet, map<int,int> &sigNet)
{
    map<string,int> mapMargTree;
    pNet->GetTopologySignature(mapMargTree);
    sigNet.clear();
    for(auto x: mapMargTree)
    {
        int tid = STSubtreeDepot::Instance().GetSTIdFor(x.first);
        sigNet[tid] = x.second;
    }
}


//***********************************************************************************
// Search for optimal network

const int DEF_MAX_NUM_NET_NODES = 1;
int AGFastCoalBuilder :: taxonOutgroup = -1;

AGFastCoalBuilder :: AGFastCoalBuilder(vector<PhylogenyTreeBasic *> &listGeneTreesIn, GenealogicalNetwork &networkInitIn, TaxaMapper &mapperIdTaxaIn) : pnetworkOpt(NULL), listGeneTrees(listGeneTreesIn), networkInit(networkInitIn), mapperTaxaIds(mapperIdTaxaIn),  maxNumMixNodes( DEF_MAX_NUM_NET_NODES ), logprobBestCurr(MAX_NEG_DOUBLE_VAL), fHeuSearch(false), fInitNetOnly(false)
{
    Init();
}

AGFastCoalBuilder :: ~AGFastCoalBuilder()
{
    //
    if( pnetworkOpt != NULL )
    {
        delete pnetworkOpt;
    }
}

void AGFastCoalBuilder :: Init()
{
    // find the maximum num of gene tree lineages
    int maxAlleleNum = 0;
    for(unsigned int i=0; i<listGeneTrees.size(); ++i)
    {
        int nl = listGeneTrees[i]->GetNumLeaves();
        if( nl > maxAlleleNum )
        {
            maxAlleleNum = nl;
        }
    }
    cout << "Maximum number of gene alleles: " << maxAlleleNum << endl;
    ApproxGTPCache2 :: Instance().Init(maxAlleleNum);
}

// search for opt network
double AGFastCoalBuilder :: Search(int numMixNodes)
{
return Search2(numMixNodes);
    //
    SetMaxNumMixNodes(numMixNodes);
    vector<int> setMixTaxa;
    FindMixTaxaImp(setMixTaxa);
    
    //
    if( pnetworkOpt == NULL )
    {
        pnetworkOpt = networkInit.Copy();
    }
    AGFastLenOpt probInitCalc( *pnetworkOpt, listGeneTrees, GetNumThreads() );
    //probInitCalc.SetMultithread(GetNumThreads());
    double probCurBest = probInitCalc.Optimize();
//cout << "Init network prob: " << probCurBest << endl;
//cout << "Init network: ";
//pnetworkOpt->DumpMargTrees(false);
    // now add the admix nodes
    for(int i=0; i<(int)setMixTaxa.size(); ++i)
    {
        probCurBest = MAX_NEG_DOUBLE_VAL;
        
        int tidMix = setMixTaxa[i];
//cout << "Search: addiing this mixing taxon: " << tidMix << endl;
        GenealogicalNetworkNode *pLeaf = pnetworkOpt->GetLeafWithTaxon(tidMix);
        YW_ASSERT_INFO( pLeaf != NULL, "Cannot be null" );
        GenealogicalNetworkBranch *pBrLeaf = pLeaf->GetAnces1();
        YW_ASSERT_INFO(pBrLeaf!=NULL, "Leaf branch cannot be null");
        vector<GenealogicalNetwork *> listNetsOneNewMix;
        pnetworkOpt->RetriveNgbrNetsOneNewMixForEdge( pBrLeaf, listNetsOneNewMix );
//cout << "Number of admixture network candidates: " << listNetsOneNewMix.size() << endl;
//for(unsigned int kkk=0; kkk<listNetsOneNewMix.size(); ++kkk)
//{
//cout << "candidate network " << kkk << ": ";
//listNetsOneNewMix[kkk]->DumpMargTrees(false);
//}
        
        // find the best net to go with
        vector<double> listLogProbs;
        for(int j=0; j<(int)listNetsOneNewMix.size(); ++j)
        {
            double probStep = ScoreForNetwork(listGeneTrees, listNetsOneNewMix[j]);
            //double probStep = ScoreForNetworkMDC(listGeneTrees, listNetsOneNewMix[j]);
            listLogProbs.push_back(probStep);
        }
//cout << "List of prob of candidate networks: ";
//DumpDoubleVec(listLogProbs);
        
        int maxProbPos = -1;
        bool fContSearch = false;
        if( listNetsOneNewMix.size() > 0 )
        {
            int maxProbPosStep = FindMaxValPositionFromList( listLogProbs );
//cout << "Max prob pos: " << maxProbPosStep << endl;
//cout << "Max prob: " << listLogProbs[maxProbPosStep] << ", probCurBest: " << probCurBest << endl;
            double probBestNgbr = listLogProbs[maxProbPosStep];
            //cout << "!Best ngbr prob: " << probBestNgbr << endl;
            
            if( probBestNgbr >= probCurBest + log(MIN_INC_RATIO) )
            {
                maxProbPos = maxProbPosStep;
                
                // now optimize branch length of the selected one
                //if( fFastSearch == false )
                //{
                probCurBest = probBestNgbr;
                //}
                //else
                //{
                //    GenealogicalNetworkLenOpt opt(*listNgbrNetsDist1[maxProbPos], listGeneTrees);
                //    double probBestNgbrLenOpt = opt.Optimize();
                    // improvement found
                    //probCurBest = probBestNgbr;
                //    probCurBest = probBestNgbrLenOpt;
                //}
//cout << "maxProbPosStep: " << maxProbPosStep << ", size of net candidates: " << listNetsOneNewMix.size() << endl;
                delete pnetworkOpt;
                pnetworkOpt = listNetsOneNewMix[maxProbPosStep];
                fContSearch = true;
//cout << "UPDATING THE CURRENT NETWORK TO: ";
// pnetworkOpt->DumpMargTrees(false);
                //pnetworkOpt->Dump();
            }
        }
        
        for(int k=0; k<(int)listNetsOneNewMix.size(); ++k)
        {
            //if( k != maxProbPos || fContSearch == false )
            if( k != maxProbPos  )
            {
                delete listNetsOneNewMix[k];
            }
        }
        if( fContSearch == false )
        {
            break;
        }
    }
    
//cout << "Perform one more round of length optimization...\n";
//exit(1);
//#if 0
    // optimize again
    AGFastLenOpt probInitCalc2( *pnetworkOpt, listGeneTrees, GetNumThreads() );
    //probInitCalc2.SetMultithread(GetNumThreads());
    probInitCalc2.SetLogProbPre(probCurBest);
    probCurBest = probInitCalc2.Optimize();
//#endif
    
    SetCurLogProb(probCurBest);
    
    // now try to improve the network
    //if( fHeuSearch == false )
    //{
        probCurBest = ExploreNgbrs();
    //}
    return probCurBest;
}

// search for opt network ver 2
// add admixture to the network and search for opt
double AGFastCoalBuilder :: Search2(int numMixNodes)
{
#if 0
    // initialize by construct net heuristically
    AGHeuristicBuilder agHeuBuilder(listGeneTrees, taxonOutgroup);
    agHeuBuilder.Build( numMixNodes );
    this->pnetworkOpt = agHeuBuilder.GetNet()->Copy();
    //return agHeuBuilder.GetCurrScore();
    AGFastLenOpt probInitCalc( *pnetworkOpt, listGeneTrees );
//cout << "here1\n";
    //probInitCalc.SetMultithread(GetNumThreads());
//cout << "here2\n";
    double probCurBest = probInitCalc.Optimize();
    cout << "Initial network probability: " << probCurBest << endl;
//cout << "Init network: ";
    
    // if only init net to be built, stop
    if( fInitNetOnly )
    {
        return probCurBest;
    }
#endif
    
    //
    SetMaxNumMixNodes(numMixNodes);
    
//cout << "Search2: numMixNodes=" << numMixNodes << endl;
//#if 0
    vector<int> setMixTaxa;
    FindMixTaxaImp(setMixTaxa);
//cout << "Search2: setMixTaxa: ";
//DumpIntVec(setMixTaxa);
    //
    if( pnetworkOpt == NULL )
    {
        pnetworkOpt = networkInit.Copy();
    }
//cout << "here\n";
//cout << "AGFastCoalBuilder :: Search2: Num of threads: " << GetNumThreads() << endl;
    //AGFastLenOpt probInitCalc( *pnetworkOpt, listGeneTrees, GetNumThreads() );
    AGFastLenOpt probInitCalc( *pnetworkOpt, listGeneTrees );
//cout << "here1\n";
    //probInitCalc.SetMultithread(GetNumThreads());
//cout << "here2\n";
    double probCurBest = probInitCalc.Optimize();
    cout << "Initial network probability: " << probCurBest << endl;
//pnetworkOpt->MapBackUserLabelsByUserLabels( mapperTaxaIds );
//cout << "Init tree: ";
//pnetworkOpt->DumpMargTrees(false);
//exit(1);
    // now add the admix nodes
    for(int i=0; i<(int)setMixTaxa.size(); ++i)
    {
        probCurBest = MAX_NEG_DOUBLE_VAL;
        
        int tidMix = setMixTaxa[i];
//cout << "Search: addiing this mixing taxon: " << tidMix << endl;
        GenealogicalNetworkNode *pLeaf = pnetworkOpt->GetLeafWithTaxon(tidMix);
        YW_ASSERT_INFO( pLeaf != NULL, "Cannot be null" );
        GenealogicalNetworkBranch *pBrLeaf = pLeaf->GetAnces1();
        YW_ASSERT_INFO(pBrLeaf!=NULL, "Leaf branch cannot be null");
        vector<GenealogicalNetwork *> listNetsOneNewMix;
        pnetworkOpt->RetriveNgbrNetsOneNewMixForEdge( pBrLeaf, listNetsOneNewMix );
//cout << "Number of admixture network candidates: " << listNetsOneNewMix.size() << endl;
//for(unsigned int kkk=0; kkk<listNetsOneNewMix.size(); ++kkk)
//{
//cout << "candidate network " << kkk << ": ";
//listNetsOneNewMix[kkk]->DumpMargTrees(false);
//}
        
        // find the best net to go with
        vector<double> listLogProbs;
        for(int j=0; j<(int)listNetsOneNewMix.size(); ++j)
        {
            double probStep = ScoreForNetwork(listGeneTrees, listNetsOneNewMix[j]);
            //double probStep = ScoreForNetworkMDC(listGeneTrees, listNetsOneNewMix[j]);
            listLogProbs.push_back(probStep);
        }
//cout << "List of prob of candidate networks: ";
//DumpDoubleVec(listLogProbs);
        
        int maxProbPos = -1;
        bool fContSearch = false;
        if( listNetsOneNewMix.size() > 0 )
        {
            int maxProbPosStep = FindMaxValPositionFromList( listLogProbs );
//cout << "Max prob pos: " << maxProbPosStep << endl;
//cout << "Max prob: " << listLogProbs[maxProbPosStep] << ", probCurBest: " << probCurBest << endl;
            double probBestNgbr = listLogProbs[maxProbPosStep];
            //cout << "!Best ngbr prob: " << probBestNgbr << endl;
            
            if( probBestNgbr >= probCurBest + log(MIN_INC_RATIO) )
            {
                maxProbPos = maxProbPosStep;
                
                // now optimize branch length of the selected one
                //if( fFastSearch == false )
                //{
                    probCurBest = probBestNgbr;
                //}
                //else
                //{
                //    GenealogicalNetworkLenOpt opt(*listNgbrNetsDist1[maxProbPos], listGeneTrees);
                //    double probBestNgbrLenOpt = opt.Optimize();
                    // improvement found
                    //probCurBest = probBestNgbr;
                //    probCurBest = probBestNgbrLenOpt;
                //}
//cout << "maxProbPosStep: " << maxProbPosStep << ", size of net candidates: " << listNetsOneNewMix.size() << endl;
                delete pnetworkOpt;
                pnetworkOpt = listNetsOneNewMix[maxProbPosStep];
                fContSearch = true;
//cout << "UPDATING THE CURRENT NETWORK TO: ";
// pnetworkOpt->DumpMargTrees(false);
                //pnetworkOpt->Dump();
            }
        }
        
        for(int k=0; k<(int)listNetsOneNewMix.size(); ++k)
        {
            //if( k != maxProbPos || fContSearch == false )
            if( k != maxProbPos  )
            {
                delete listNetsOneNewMix[k];
            }
        }
        if( fContSearch == false )
        {
            break;
        }
    }
    
//cout << "Perform one more round of length optimization...\n";
//exit(1);
//#if 0
    
// hack: don't search for neighbrs
//return probCurBest;
    
//cout << "Opt br length again...\n";
//pnetworkOpt->MapBackUserLabelsByUserLabels( mapperTaxaIds );
//cout << "After adding admixture, network has  probabilty: " << probCurBest << " and with marginal trees: \n";
//pnetworkOpt->DumpMargTrees(false);
    
    // optimize again
    //AGFastLenOpt probInitCalc2( *pnetworkOpt, listGeneTrees, GetNumThreads() );
    AGFastLenOpt probInitCalc2( *pnetworkOpt, listGeneTrees );
    //probInitCalc2.SetMultithread(GetNumThreads());
    probInitCalc2.SetLogProbPre(probCurBest);
    probCurBest = probInitCalc2.Optimize();
//#endif
    
//#endif
    
    
//cout << "**** After branch length optimization, initial network has prob: " << probCurBest << " marginal trees: \n";
//pnetworkOpt->DumpMargTrees(false);
    /*
cout << "re-dump: \n";
pnetworkOpt->DumpMargTrees(false);
cout << "re-re-dump: \n";
pnetworkOpt->DumpMargTrees(false);
cout << "(4)re-re-dump: \n";
pnetworkOpt->DumpMargTrees(false);
cout << "(5)re-re-dump: \n";
pnetworkOpt->DumpMargTrees(false);
exit(1);  */
    
    SetCurLogProb(probCurBest);
    
    // if only init net to be built, stop
    if( fInitNetOnly )
    {
        return probCurBest;
    }
    
    // now try to improve the network
    //if( fHeuSearch == false )
    //{
        probCurBest = ExploreNgbrs();
    //}
    return probCurBest;
}

void AGFastCoalBuilder :: InfMixPops( set<string> &setMixTaxaUser )
{
    //
    setMixTaxaUser.clear();
    vector<int> setMixTaxa;
    FindMixTaxaImp(setMixTaxa);
    for(vector<int> :: iterator it = setMixTaxa.begin(); it != setMixTaxa.end(); ++it)
    {
        string midUser = mapperTaxaIds.GetString( *it );
        setMixTaxaUser.insert(midUser);
    }
}

//***********************************************************************************

void AGFastCoalBuilder :: FindMixTaxaImp( vector<int> &setMixTaxa )
{
    //
    set<int> setTaxaCurr;
    GetAllTaxa(setTaxaCurr);
    
    // now search for mix taxa one by one
    setMixTaxa.clear();
    
    // first add all the user-selected taxa
    for(set<int> :: iterator it = taxaMixInit.begin(); it != taxaMixInit.end(); ++it )
    {
        setMixTaxa.push_back(*it);
        
        // remove from current list
        setTaxaCurr.erase(*it);
    }
    
//cout << "maxNumMixNodes: " << maxNumMixNodes << endl;
    while( setMixTaxa.size() < this->maxNumMixNodes )
    {
        //for(int na=0; na<this->maxNumMixNodes; ++na)
        //{
        int taxonMix = FindOneMixTaxonFrom(setTaxaCurr);
        //cout << "FindOneMixTaxonFrom: taxonMix: " << taxonMix << endl;
        if( taxonMix < 0 )
        {
//cout << "Stop: no more admixed taxon.\n";
            break;
        }
        setMixTaxa.push_back(taxonMix);
        setTaxaCurr.erase(taxonMix);
        
        if(setMixTaxa.size() >= this->maxNumMixNodes)
        {
            break;
        }
        //}
    }
    
//cout << "FindMixTaxaImp: Set of admixed taxa: ";
//DumpIntVec(setMixTaxa);
}

int AGFastCoalBuilder :: FindOneMixTaxonFrom( const set<int> &setTaxaCurr )
{
//cout << "FindOneMixTaxonFrom: setTaxaCurr:";
//DumpIntSet(setTaxaCurr);
    // consider all taxon in the set one by one to see which one fits more a mixture
    if( setTaxaCurr.size() <= 0 )
    {
        // nothing left
//cout << "nothing left\n";
        return -1;
    }
    vector<int> listMultiplicityTrees;
    for(int tr=0; tr<(int)listGeneTrees.size(); ++tr)
    {
        listMultiplicityTrees.push_back(1);
    }
    
    // get the overall MDC
    const int DEF_MDC_LEVEL= 1;
    DeepCoalescence mdcOrigCalc(listGeneTrees, DEF_MDC_LEVEL);
    mdcOrigCalc.SetMultiplictyofTrees(listMultiplicityTrees);
    int mdcOrig = mdcOrigCalc.FindMDC();
//cout << "MCD orig: " << mdcOrig << endl;
    
    //int mdcMin = HAP_MAX_INT;
    //int res = -1;
    vector<pair<int,int> > listMDCLevels;
    vector<int> listTaxa;
    for(set<int> :: iterator it = setTaxaCurr.begin(); it != setTaxaCurr.end(); ++it )
    {
        int taxonCurr = *it;
        
        // disallow og
        if( taxonCurr == taxonOutgroup)
        {
            continue;
        }
        
        set<int> setTaxaStep = setTaxaCurr;
        setTaxaStep.erase(taxonCurr);
        
        vector<PhylogenyTreeBasic *> listPhyTreesSub;
        GetReducedGenetreesForSubsetTaxa( setTaxaStep, listPhyTreesSub );
        
        // convert tree to consecutive zero-based
        TaxaMapper taxonMapper;
        ConvPhyloTreesToZeroBasedId( listPhyTreesSub, &taxonMapper);
//cout << "Curr taxon: " << taxonCurr << ", set of taxa remaining: ";
//DumpIntSet( setTaxaStep );
//for(int tr=0; tr<(int)listPhyTreesSub.size(); ++tr)
//{
//string strNW;
//listPhyTreesSub[tr]->ConsNewickSorted(strNW);
//cout << "Tree: " << strNW << endl;
//}
        
        //
        DeepCoalescence mdc(listPhyTreesSub, DEF_MDC_LEVEL);
        mdc.SetMultiplictyofTrees(listMultiplicityTrees);
        int mdcVal = mdc.FindMDC();
        YW_ASSERT_INFO(listPhyTreesSub.size() > 0, "Cannot be empty list");
        int numLvsStep = listPhyTreesSub[0]->GetNumLeaves();
        pair<int,int> pp(mdcVal, numLvsStep);
        listMDCLevels.push_back(pp);
        listTaxa.push_back(taxonCurr);
        //if( mdcVal < mdcMin)
        //{
        //    mdcMin = mdcVal;
        //    res = taxonCurr;
        //}
        
//cout << "MDC: " << mdcVal;
//cout << "  for curr taxon: " << taxonCurr << endl;
        
//cout << "Curr taxon: " << taxonCurr << ", set of taxa remaining: ";
//DumpIntSet( setTaxaStep );
        for(int tr=0; tr<(int)listPhyTreesSub.size(); ++tr)
        {
//string strNW;
//listPhyTreesSub[tr]->ConsNewickSorted(strNW);
//cout << "Tree: " << strNW << endl;
            delete listPhyTreesSub[tr];
        }
    }
    
    int res = -1;
    double mdcUnitOpt = -1.0;
    for(int i=0; i<(int)listMDCLevels.size(); ++i)
    {
        double mdcUnitStep = ((double)(  mdcOrig-listMDCLevels[i].first ))/listMDCLevels[i].second;
        if(mdcUnitOpt < mdcUnitStep)
        {
            mdcUnitOpt = mdcUnitStep;
            res = listTaxa[i];
        }
    }
    //YW_ASSERT_INFO(res >= 0, "Fail to find initial taxa");
//cout << "** MDC-opt admixure taxon: " << res << endl;
//cout << "mdcUnitOpt: " << mdcUnitOpt << endl;
    return res;
}

double AGFastCoalBuilder :: ScoreForNetwork(const vector<PhylogenyTreeBasic *> &listTrees, GenealogicalNetwork *pNetCurr)
{
    // score
    vector<PhylogenyTreeBasic *> &listTreesUse = const_cast<vector<PhylogenyTreeBasic *> & >(listTrees);
    //AGFastLenOpt opt(*pNetCurr, listTreesUse, GetNumThreads() );
    AGFastLenOpt opt(*pNetCurr, listTreesUse );
    //opt.SetMultithread(GetNumThreads());
    return opt.CalcProb();
}

double AGFastCoalBuilder :: ScoreForNetwork2(const vector<PhylogenyTreeBasic *> &listTrees, GenealogicalNetwork *pNetCurr)
{
    // score
    vector<PhylogenyTreeBasic *> &listTreesUse = const_cast<vector<PhylogenyTreeBasic *> & >(listTrees);
    //AGFastLenOpt opt(*pNetCurr, listTreesUse, GetNumThreads() );
    AGFastLenOpt opt(*pNetCurr, listTreesUse );
    //opt.SetMultithread(GetNumThreads());
    double likeli0 = opt.CalcProb();
    double likeli1 = opt.Optimize();
cout << "Un-optmized likelihood: " << likeli0 << ", optimized likelihood: " << likeli1 << endl;
    return likeli1;
}

double AGFastCoalBuilder :: ScoreForNetworkMDC(const vector<PhylogenyTreeBasic *> &listTrees, GenealogicalNetwork *pNetCurr)
{
    // get marginal trees embedded in the network
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTreesWithFreq;
    pNetCurr->RetriveAllMarginalTrees( mapMargTreesWithFreq );
    
    vector<MarginalTree *> listMargTrees;
    for(map<GenealogicalNetworkMTreeCode,GenealogicalNetworkMTreeInfo> :: iterator it=mapMargTreesWithFreq.begin(); it!=mapMargTreesWithFreq.end(); ++it)
    {
        listMargTrees.push_back(it->second.GetMargTree());
    }
    
    // for each tree, take the best (lowest) MDC
    int res = 0;
    for(int i=0; i<(int)listTrees.size(); ++i)
    {
        int costStep = HAP_MAX_INT;
        for(int j=0; j<(int)listMargTrees.size(); ++j)
        {
            int cstep = ScoreMDCGTandST(listTrees[i], listMargTrees[j]);
            if( cstep < costStep )
            {
                costStep = cstep;
            }
        }
        res += costStep;
    }
    
    // free trees
    for(auto x : listMargTrees )
    {
        delete x;
    }
    
    //tbd;
    return -1.0*res;
}

// support multithreading
static void UtilExploreNgbrNetsMT( int tid, AGFastCoalBuilder *pCoalBuilder, int numTotNgbrs, vector<bool> *pvecProcessed, vector<GenealogicalNetwork *> *plistNgbrNetsDist1, vector<GenealogicalNetwork *> *plistNgbrNetsDist1New, vector<int> *plistNgbrNetPos, vector<PhylogenyTreeBasic *> *plistGeneTrees, vector<double> *plistLogProbNgbrNets, vector<set<GenealogicalNetworkBranch *> > *plistNgbrNetsDist1CritBrs )
{
    int indNet = 0;
    static std::mutex mut;
    while(indNet < numTotNgbrs)
    {
        // find a node that is not being processed
        mut.lock();
        while( indNet < numTotNgbrs && (*pvecProcessed)[indNet] == true )
        {
            ++indNet;
        }
        if( indNet < numTotNgbrs )
        {
            (*pvecProcessed)[indNet] = true;
        }
        mut.unlock();
        if( indNet >= numTotNgbrs )
        {
            break;
        }
        
        // real work
        //YW_ASSERT_INFO( i<listNgbrNetPos.size(), "Overflow444" );
        int netIndex = (*plistNgbrNetPos)[indNet];
#if 0
        cout << "&&&&&&&&&&&&&&&&&&Processing ngbr network " << netIndex << endl;
        //listNgbrNetsDist1[netIndex]->DumpMargTrees(false);
        //listNgbrNetsDist1[netIndex]->Dump();
#endif

        // mark current processed
        pCoalBuilder->MarkNetProcessed((*plistNgbrNetsDist1)[netIndex]);
        
//cout << "Creating length optimization...\n";
        // don't use multithreading for individual net
        AGFastLenOpt *popt = new AGFastLenOpt(*(*plistNgbrNetsDist1)[netIndex], *plistGeneTrees, tid );
        //opt.SetMultithread(GetNumThreads());
        //double logprobBest = opt.Optimize();
        
        double probStep = 0.0;
        if( fFastSearch == false )
        {
            popt->SetCritcalBrs((*plistNgbrNetsDist1CritBrs)[netIndex]);
//cout << "Number of critical branches: " << listNgbrNetsDist1CritBrsNew[netIndex].size() << endl;
//cout << "Now optimize branch length...\n";
            probStep = popt->Optimize();
//cout << "(1a) Network ngbr " << i << ": optimized likelihood: " << logprobBest << endl;
        }
        else
        {
//cout << "Now computing prob without opt branch length...\n";
            // YW: Feb 18, 2016: don't optimize branch length after a change of topology
            probStep = popt->CalcProb();
//cout << "(1) Network ngbr " << i << ": un-optimized likelihood: " << logprobChange << endl;
        }
        
        // set prob
        (*plistLogProbNgbrNets)[indNet] = probStep;

        delete popt;
    }
}


double AGFastCoalBuilder :: ExploreNgbrs()
{
    YW_ASSERT_INFO(pnetworkOpt != NULL, "Opt network is not initialized");
    //GenealogicalNetFilter netFilter(pnetworkOpt);
    
    // mark current processed
    MarkNetProcessed(pnetworkOpt);
    
    // clear out previous stored tree prob. YW: TBD
    //AGGeneTreeProbDepot::Instance().Clear();
    
    //double probCurBest = MAX_NEG_DOUBLE_VAL;
    double probCurBest = GetCurrBestLogprob();
    while(true)
    {
        // mark current processed
        //MarkNetProcessed(pnetworkOpt);
        
#if 0
        cout << "--Processing current network (w/ prob: " << probCurBest << "): \n";
//cout << "   num of reticulate nodes: " << pnetworkOpt->GetNumMixNodes() << endl;
//pnetworkOpt->MapBackUserLabelsByUserLabels( mapperTaxaIds );
cout << "Current network: \n";
pnetworkOpt->DumpMargTrees(false);
pnetworkOpt->DumpAdmixNodeIds();
//pnetworkOpt->Dump();
#endif
        //
        vector<GenealogicalNetwork *> listNgbrNetsDist1;
        vector<set<GenealogicalNetworkBranch *> > listNgbrNetsDist1CritBrs;
        int og = AGFastCoalBuilder :: GetOutgroup();
        //int og = -1;
//cout << "og: " << og << endl;
        //pnetworkOpt->RetriveNgbrNetsOneEvt( listNgbrNetsDist1, og );
        pnetworkOpt->RetriveNgbrNetsOneEvt2( listNgbrNetsDist1, og, &listNgbrNetsDist1CritBrs );
        //pnetworkOpt->RetriveNgbrNetsOneEvt2( listNgbrNetsDist1, og );
//cout << "listNgbrNetsDist1: size: " << listNgbrNetsDist1.size() << endl;
//cout << "Found ngbr networks: " << listNgbrNetsDist1.size() << ", critical branch sets: " << listNgbrNetsDist1CritBrs.size() << endl;
        
        // first find new networks
        vector<int> listNgbrNetPos;
        vector<GenealogicalNetwork *> listNgbrNetsDist1New;
        //vector<set<GenealogicalNetworkBranch *> > listNgbrNetsDist1CritBrsNew;
        for(int i=0; i<(int)listNgbrNetsDist1.size(); ++i)
        {
            YW_ASSERT_INFO( listNgbrNetsDist1[i] != NULL, "BAD"  );
//cout << "^^^ raw ngbr network " << i << endl;
//listNgbrNetsDist1[i]->DumpMargTrees(false);
            //listNgbrNetsDist1[i]->Dump();
            
            // skip any net violating outgroup
            if( listNgbrNetsDist1[i]->CheckOutgroup(taxonOutgroup) == false )
            {
                continue;
            }
            
            // make sure it has no trivial cycles
            if( listNgbrNetsDist1[i]->CheckTrivialCycles() )
            {
                //
                continue;
            }
            
            // if this net has been processed, skip
            if( IsNetProcBefore(listNgbrNetsDist1[i]) )
            {
                AGProcessedNetsDepot::Instance().SkipNet(listNgbrNetsDist1[i]);
                continue;
            }
            
            //bool fFoundBefore = IsNetProcBefore( listNgbrNetsDist1[i] );
            bool fOGOK = IsNetworkOGGood( *listNgbrNetsDist1[i] );
            
            if( fOGOK )
            {
                // make sure this net has the desired number of network nodes
                if( listNgbrNetsDist1[i]->GetNumMixNodes() != this->maxNumMixNodes )
                {
                    continue;
                }
                
//cout << "^^^ keeping raw ngbr network " << i << endl;
//listNgbrNetsDist1[i]->DumpMargTrees(false);
                listNgbrNetsDist1New.push_back( listNgbrNetsDist1[i] );
                listNgbrNetPos.push_back(i);
                //listNgbrNetsDist1CritBrsNew.push_back(listNgbrNetsDist1CritBrs[i]);
                
                // mark how many nets we have processed
                //infoSearch.ProcNet( listNgbrNetsDist1[i] );
            }
            else
            {
                //infoSearch.SkipNet( listNgbrNetsDist1[i] );
                // free it
            }
        }
        
cout << "Number of neighboring nets to evaluate: " << listNgbrNetsDist1New.size() << endl;
        // find ngbr networks (no adding new mix nodes)
        vector<double> listLogProbNgbrNets;
        
        if( GetNumThreads() > 1 )
        {
            // multi-threading
            // first initialize prob vector
            listLogProbNgbrNets.resize( listNgbrNetsDist1New.size() );
            
            // multithreading
            int numThreadsUse = GetNumThreads();
            if( (int)listNgbrNetsDist1New.size() < numThreadsUse )
            {
                numThreadsUse = listNgbrNetsDist1New.size();
            }
//cout << "MT: num of threads to use: " << numThreadsUse << endl;
            vector<bool> vecProcessed(listNgbrNetsDist1New.size());
            for(int i=0; i<(int)listNgbrNetsDist1New.size(); ++i)
            {
                vecProcessed[i] = false;
            }
            // vector for a node being processed
            vector<thread *> listPtrThreads;
            for(int t=0; t<numThreadsUse; ++t)
            {
                thread *pthr = new thread(UtilExploreNgbrNetsMT, t, this, listNgbrNetsDist1New.size(), &vecProcessed, &listNgbrNetsDist1, &listNgbrNetsDist1New, &listNgbrNetPos, &listGeneTrees, &listLogProbNgbrNets, &listNgbrNetsDist1CritBrs );
                listPtrThreads.push_back(pthr);
            }
            
            // wait to join
            for(unsigned int j=0; j<listPtrThreads.size(); ++j)
            {
              listPtrThreads[j]->join();
            }
            for(unsigned int j=0; j<listPtrThreads.size(); ++j)
            {
              delete listPtrThreads[j];
            }
            listPtrThreads.clear();
            
        }
        else
        {
            for(int i=0; i<(int)listNgbrNetsDist1New.size(); ++i)
            {
                YW_ASSERT_INFO( i<listNgbrNetPos.size(), "Overflow444" );
                int netIndex = listNgbrNetPos[i];
#if 0
                cout << "&&&&&&&&&&&&&&&&&&Processing ngbr network " << netIndex << endl;
                listNgbrNetsDist1[netIndex]->DumpMargTrees(false);
                listNgbrNetsDist1[netIndex]->DumpAdmixNodeIds();
                //listNgbrNetsDist1[netIndex]->Dump();
#endif
                bool fSearch = true;
                
                if( fHeuSearch == true )
                {
                    //    fSearch = netFilter.ScoreNetwork( listNgbrNetsDist1[netIndex] );
                }
                
                
                if( fSearch == true )
                {
                    // mark current processed
                    MarkNetProcessed(listNgbrNetsDist1[netIndex]);
                    
                    //cout << "Creating length optimization...\n";
                    AGFastLenOpt opt(*listNgbrNetsDist1[netIndex], listGeneTrees );
                    //opt.SetMultithread(GetNumThreads());
                    //double logprobBest = opt.Optimize();
                    
                    if( fFastSearch == false )
                    {
                        opt.SetCritcalBrs(listNgbrNetsDist1CritBrs[netIndex]);
                        //cout << "Number of critical branches: " << listNgbrNetsDist1CritBrsNew[netIndex].size() << endl;
                        //cout << "Now optimize branch length...\n";
                        double logprobBest = opt.Optimize();
                        listLogProbNgbrNets.push_back( logprobBest );
//cout << "(1a) Network ngbr " << i << ": optimized likelihood: " << logprobBest << endl;
                    }
                    else
                    {
                        //cout << "Now computing prob without opt branch length...\n";
                        // YW: Feb 18, 2016: don't optimize branch length after a change of topology
                        double logprobChange = opt.CalcProb();
                        //listLogProbNgbrNets.push_back( logprobBest );
                        listLogProbNgbrNets.push_back(logprobChange);
//cout << "(1) Network ngbr " << i << ": un-optimized likelihood: " << logprobChange << endl;
                    }
                    
                    
                    // mark it
                    //bool fNoWorse = (logprobChange >= probCurBest);
                    //netFilter.MarkNetwork( listNgbrNetsDist1[netIndex], fNoWorse );
                }
            }
        }
//cout << "**** Number of neighbr nets: " << listLogProbNgbrNets.size() << endl;
        // if nothing left, done
        int maxProbPos = -1;
        bool fContSearch = false;
        if( listLogProbNgbrNets.size() > 0 )
        {
            int maxProbPosStep = FindMaxValPositionFromList( listLogProbNgbrNets );
            double probBestNgbr = listLogProbNgbrNets[maxProbPosStep];
//cout << "!Best ngbr prob: " << probBestNgbr << ", probCurBest: " << probCurBest << endl;
            
            if( probBestNgbr >= probCurBest + log(MIN_INC_RATIO) )
            {
                maxProbPos = listNgbrNetPos[maxProbPosStep];
                
                // now optimize branch length of the selected one
                if( fFastSearch == false )
                {
                    probCurBest = probBestNgbr;
//cout << "-- Likelihood improved to: " << probCurBest << endl;
                }
                else
                {
                    //
                    //AGFastLenOpt opt(*listNgbrNetsDist1[maxProbPos], listGeneTrees, GetNumThreads() );
                    AGFastLenOpt opt(*listNgbrNetsDist1[maxProbPos], listGeneTrees );
                    //opt.SetMultithread(GetNumThreads());
                    double probBestNgbrLenOpt = opt.Optimize();
                    
                    // improvement found
                    //probCurBest = probBestNgbr;
                    probCurBest = probBestNgbrLenOpt;
                }
                
                delete pnetworkOpt;
                pnetworkOpt = listNgbrNetsDist1[maxProbPos];
                fContSearch = true;
                //cout << "UPDATING THE CURRENT NETWORK TO: ";
                //pnetworkOpt->Dump();
                
                //netFilter.AdoptNet(pnetworkOpt);
            }
        }
        //YW_ASSERT_INFO(maxProbPos >=0, "WRONG");
        int numNetsDeleted = 0;
        for(int i=0; i<(int)listNgbrNetsDist1.size(); ++i)
        {
            if( i != maxProbPos || fContSearch == false )
            {
                delete listNgbrNetsDist1[i];
                ++numNetsDeleted;
            }
        }
//cout << "Number of ngbr nets deleted: " << numNetsDeleted << endl;
//cout << "@@@@@@@@@@@@@@@ Current best log-likelihood: " << probCurBest << endl;
        //cout << "Best network: ";
        //pnetworkOpt->DumpMargTrees(false);
        if( fContSearch == false )
        {
            break;
        }
    }
    
    return probCurBest;
}

//***********************************************************************************

GenealogicalNetwork * AGFastCoalBuilder :: GetBestNet()
{
    return pnetworkOpt;
}

void AGFastCoalBuilder :: SetOutgroup(int r)
{
    //cout << "Interval outgroup id: " << r << endl;
    taxonOutgroup = r;
}
int AGFastCoalBuilder :: GetOutgroup()
{
    return taxonOutgroup;
}

void AGFastCoalBuilder :: GetAllTaxa(set<int> &setTaxa) const
{
    YW_ASSERT_INFO(listGeneTrees.size() > 0, "Must have some trees");
    listGeneTrees[0]->GetLeafIntLabels( setTaxa );
}

void AGFastCoalBuilder :: GetReducedGenetreesForSubsetTaxa( const set<int> &setTaxaSub, vector<PhylogenyTreeBasic *> &listGeneTreesSub ) const
{
    //
    listGeneTreesSub.clear();
    for(int tr=0; tr<(int)listGeneTrees.size(); ++tr)
    {
        PhylogenyTreeBasic *ptreesub = ConsPhyTreeSubsetTaxa( listGeneTrees[tr], setTaxaSub );
        listGeneTreesSub.push_back(ptreesub);
    }
}

bool AGFastCoalBuilder :: IsNetworkOGGood( const GenealogicalNetwork &netTest )
{
    // if outgroup is set, then each embedded tree should be consistent with the outgroup: yes each
    if( taxonOutgroup < 0 )
    {
        return true;
    }
    //cout << "IsNetworkOGGood: network is ";
    //netTest.DumpMargTrees(false);
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTreesWithFreq;
    netTest.RetriveAllMarginalTrees( mapMargTreesWithFreq );
    bool res = true;
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: iterator it = mapMargTreesWithFreq.begin(); it != mapMargTreesWithFreq.end(); ++it )
    {
        MarginalTree *ptree = it->second.GetMargTree();
        if( ptree->IsOutgroup( taxonOutgroup ) == false )
        {
            res = false;
            break;
        }
    }
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> :: iterator it = mapMargTreesWithFreq.begin(); it != mapMargTreesWithFreq.end(); ++it )
    {
        MarginalTree *ptree = it->second.GetMargTree();
        delete ptree;
    }
    
    return res;
}

int AGFastCoalBuilder :: ScoreMDCGTandST( PhylogenyTreeBasic *pGeneTree, MarginalTree *pST )
{
    // get all clades of ST
    vector<set<int> > leafNodeLabels ;
    pST->ConsDecedentLeavesInfoLabels( leafNodeLabels );
    
    // cons string-based labels
    vector<set<string> > leafNodeLabelsUse;
    for(int i=0;i<(int)leafNodeLabels.size(); ++i)
    {
        //
        set<string> ss;
        for(set<int>::iterator it=leafNodeLabels[i].begin(); it!=leafNodeLabels[i].end(); ++it)
        {
            string str = std::to_string( *it );
            ss.insert(str);
        }
        leafNodeLabelsUse.push_back(ss);
    }
    
    // now get the
    int res = 0;
    
    //
    set<set<TreeNode *> > setClades;
    pGeneTree->GetAllCladeNodess(setClades);
    for(int i=0;i<(int)leafNodeLabelsUse.size(); ++i)
    {
#if 0
cout << "Taxa: ";
for(set<string>:: iterator it=leafNodeLabelsUse[i].begin(); it!=leafNodeLabelsUse[i].end();++it)
{
cout << *it << " ";
}
cout << endl;
#endif
        set<TreeNode *> setLvNodes;
        pGeneTree->GetLeavesWithLabels( leafNodeLabelsUse[i], setLvNodes );
        
        //
        set< set<TreeNode *> > setSubtreeClades;
        PhylogenyTreeBasic :: GroupLeavesToSubtrees( setLvNodes,  setClades, setSubtreeClades );
        res += setSubtreeClades.size()-1;
//cout << "Curr MDC score: " << setSubtreeClades.size()-1 << endl;
    }
    string strNW;
    pGeneTree->ConsNewickSorted(strNW);
//cout << "For gene tree: " << strNW << ", marginalt ree: " << pST->GetNewickSorted(false) << ", MDC = " << res << endl;
    return res;
}

bool AGFastCoalBuilder :: IsNetProcBefore(GenealogicalNetwork *pNet) const
{
return false;
    return AGProcessedNetsDepot::Instance().IsNetProcessed(pNet);
}

void AGFastCoalBuilder :: MarkNetProcessed(GenealogicalNetwork *pNet) const
{
return;
    AGProcessedNetsDepot::Instance().MarkNetProcessed(pNet);
}
