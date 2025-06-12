//
//  AGHeuristicSearch.cpp
//  
//
//  Created by Yufeng Wu on 2/4/25.
//

#include "AGHeuristicSearch.hpp"
#include "GenealogicalNetwork.h"
#include "PhylogenyTreeBasic.h"
#include "GenealogicalNetworkHeuSearch.h"
#include "Utils4.h"
#include "UtilsNumerical.h"
#include "AGFastCoalBuilder.hpp"

//***********************************************************************************
// Evaluate how a gene tree fits AG


AGHeuristicSearchEval :: AGHeuristicSearchEval(GenealogicalNetwork &agCurrIn) : agCurr(agCurrIn)
{
    ProcessAG();
}

double AGHeuristicSearchEval :: EvalTree( PhylogenyTreeBasic *pTree )
{
    // first, collect allele counts for all taxa
    map<int,int> mapAlleleCnts;
    vector<int> listLabels;
    pTree->GetAllLeafIntLabeles(listLabels);
    int numLeaves = listLabels.size();
    for(unsigned int i=0; i<listLabels.size();++i)
    {
        ++mapAlleleCnts[ listLabels[i] ];
    }
//cout << "Allele counts: ";
//for(auto it = mapAlleleCnts.begin(); it != mapAlleleCnts.end(); ++it)
//{
//cout << "allele " << it->first << ": #" << it->second  << "   ";
//}
//cout << endl;
    
    // deal with each clade
    double res = 0.0;
    
    set<set<TreeNode *> > setClades;
    pTree->GetAllCladeNodess(setClades);
    for(auto it = setClades.begin(); it != setClades.end(); ++it)
    {
        // ignore all leaves
        if( it->size() <= 1 || it->size() == numLeaves )
        {
            continue;
        }
        
        // calculate a profile based on frequency
        map<int,int> mapAlleleCntClade;
        set<int> clade;
        vector<int> cladeVec;
        for(auto x : *it )
        {
            int lbl = x->GetIntLabel();
            ++mapAlleleCntClade[lbl];
            clade.insert(lbl);
            cladeVec.push_back(lbl);
        }
        
        // ignore intra-taxa coalescent
        if( clade.size() <= 1 )
        {
            continue;
        }
        
        //
        double cladeFreq = 1.0;
        for(auto it2 = mapAlleleCntClade.begin(); it2 != mapAlleleCntClade.end(); ++it2)
        {
            int cnt = mapAlleleCnts[it2->first];
            YW_ASSERT_INFO(cnt >=1 && cnt >= it2->second, "WRONG 2222");
            cladeFreq *= (1.0*it2->second)/cnt;
        }
        
        // check all branches
        double scoreStep = 0.0;
        for(auto it3 = mapAGInfo.begin(); it3 != mapAGInfo.end(); ++it3)
        {
            if( it3->second.find(clade) != it3->second.end() )
            {
                scoreStep += it3->second[clade] * cladeFreq;
//cout << "----Find a hit with clade score in AG: " << it3->second[clade] << endl;
            }
        }
        
        res += scoreStep;
        
//cout << "scoreStep: " << scoreStep << ", cladeFreq: " << cladeFreq << "  CaldeVec: ";
//DumpIntVec(cladeVec);
    }
    
    return res;
}

//***********************************************************************************
// evaluate the AG
void AGHeuristicSearchEval :: ProcessAG()
{
    //
    mapAGInfo.clear();
    // traverse AG
    vector<GenealogicalNetworkBranch *> setAllBranches;
    agCurr.GetAllBranchesBottomUp( setAllBranches );
    for(unsigned int i=0; i<setAllBranches.size(); ++i)
    {
        GenealogicalNetworkBranch *br = setAllBranches[i];
        GenealogicalNetworkNode *pn = br->GetDestNode();
        // if leaf
        if( pn->IsLeaf())
        {
            int tid = pn->GetTaxonId();
            set<int> ss;
            ss.insert(tid);
            mapAGInfo[br][ss] = 1.0;
        }
        else if( pn->IsMixNode() == false )
        {
            // internal node, coalesce
            GenealogicalNetworkBranch *pbrChild1 = pn->GetLeftDesc();
            GenealogicalNetworkBranch *pbrChild2 = pn->GetRightDesc();
            YW_ASSERT_INFO(pbrChild1!=NULL, "pbrChild1 is NULL");
            YW_ASSERT_INFO(pbrChild2!=NULL, "pbrChild2 is NULL");
            // coal these two
            for( auto it1 = mapAGInfo[pbrChild1].begin(); it1 != mapAGInfo[pbrChild1].end(); ++it1 )
            {
                for( auto it2 = mapAGInfo[pbrChild2].begin(); it2 != mapAGInfo[pbrChild2].end(); ++it2 )
                {
                    set<int> ss = it1->first;
                    UnionSets(ss, it2->first);
                    
                    set<int> sint;
                    JoinSets(it1->first, it2->first, sint);
                    
                    // if disjoint, then multiply; otherwise, add up
                    double score = it1->second * it2->second;
                    if( sint.size() > 0 )
                    {
                        score = it1->second + it2->second;
                    }
                    auto it3 = mapAGInfo[br].find(ss);
                    if( it3 != mapAGInfo[br].end() )
                    {
                        mapAGInfo[br][ss] += score;
                    }
                    else
                    {
                        mapAGInfo[br][ss] = score;
                    }
                }
            }
        }
        else
        {
            // mix node
            double mr = pn->GetMixRatio();
            bool fLeft = pn->IsLeftAncesBranch(br);
            if( fLeft == false )
            {
                mr = 1.0 - mr;  // the right branch mr = 1.0 - left mr
            }
            //
            GenealogicalNetworkBranch *pbrBelow = pn->GetLeftDesc();
            if( pbrBelow == NULL )
            {
                pbrBelow = pn->GetRightDesc();
            }
            YW_ASSERT_INFO(pbrBelow != NULL, "Canot be NULL678");
            YW_ASSERT_INFO( mapAGInfo.find(pbrBelow) != mapAGInfo.end() && mapAGInfo[pbrBelow].size() > 0, "WRONG678" );
            for( auto it1 = mapAGInfo[pbrBelow].begin(); it1 != mapAGInfo[pbrBelow].end(); ++it1 )
            {
                mapAGInfo[br][it1->first] = it1->second * mr;
            }
        }
    }
#if 0
    // dump
    cout << "******* List of branch clades: \n";
    for(auto it1 = mapAGInfo.begin(); it1 != mapAGInfo.end(); ++it1)
    {
        cout << "Branch: ";
        it1->first->Dump();
        cout << "   set of branch clades: \n";
        for(auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
            cout << "[" << it2->second << "]: ";
            DumpIntSet(it2->first);
        }
    }
#endif
}


//***********************************************************************************
// Evaluate how a gene tree fits AG based on trees contained in AG


AGHeuristicSearchTreeBaseEval :: AGHeuristicSearchTreeBaseEval(GenealogicalNetwork &agCurrIn): agCurr(agCurrIn)
{
    ProcessAG();
}

AGHeuristicSearchTreeBaseEval :: ~AGHeuristicSearchTreeBaseEval()
{
}

double AGHeuristicSearchTreeBaseEval :: EvalTree( PhylogenyTreeBasic *pTree )
{
//string strTree;
//cout << "AGHeuristicSearchTreeBaseEval :: EvalTree: tree:";
//pTree->ConsNewickSorted(strTree);
//cout << strTree << endl;
    // first, collect allele counts for all taxa
    map<int,int> mapAlleleCnts;
    vector<int> listLabels;
    pTree->GetAllLeafIntLabeles(listLabels);
    int numLeaves = listLabels.size();
    for(unsigned int i=0; i<listLabels.size();++i)
    {
        ++mapAlleleCnts[ listLabels[i] ];
    }
//cout << "Allele counts: ";
//for(auto it = mapAlleleCnts.begin(); it != mapAlleleCnts.end(); ++it)
//{
//cout << "allele " << it->first << ": #" << it->second  << "   ";
//}
//cout << endl;
    
    // clade info in pTree
    vector<pair<set<int>, double> > listGTreeInfo;
    
    set<set<TreeNode *> > setClades;
    pTree->GetAllCladeNodess(setClades);
    for(auto it = setClades.begin(); it != setClades.end(); ++it)
    {
        // ignore all leaves
        if( it->size() <= 1 || it->size() == numLeaves )
        {
            continue;
        }
        
        // calculate a profile based on frequency
        map<int,int> mapAlleleCntClade;
        set<int> clade;
        vector<int> cladeVec;
        for(auto x : *it )
        {
            int lbl = x->GetIntLabel();
            ++mapAlleleCntClade[lbl];
            clade.insert(lbl);
            cladeVec.push_back(lbl);
        }
        
        // ignore intra-taxa coalescent
        if( clade.size() <= 1 )
        {
            continue;
        }
        
        //
        double cladeFreq = 1.0;
        for(auto it2 = mapAlleleCntClade.begin(); it2 != mapAlleleCntClade.end(); ++it2)
        {
            int cnt = mapAlleleCnts[it2->first];
            YW_ASSERT_INFO(cnt >=1 && cnt >= it2->second, "WRONG 2222");
            cladeFreq *= (1.0*it2->second)/cnt;
        }
        
        listGTreeInfo.push_back(std::make_pair(clade, cladeFreq));
//cout << "^cladeFreq:" << cladeFreq << " for clade: ";
//DumpIntSet(clade);
    }
    
    // for all found GT clade
    // now check all contained trees
    double scoreMax = 0.0;
    for(unsigned int i=0; i<listContainedAGTreeInfo.size(); ++i)
    {
        double scoreStep = 0.0;
        for(auto xInfo : listGTreeInfo)
        {
            if( listContainedAGTreeInfo[i].find(xInfo.first) != listContainedAGTreeInfo[i].end() )
            {
                scoreStep += xInfo.second;
//cout << "----Find a hit with clade score in AG: " << xInfo.second << endl;
            }
        }
        
        if( scoreMax < scoreStep )
        {
            scoreMax = scoreStep;
        }
        
        //cout << "scoreStep: " << scoreStep << "  CaldeVec: ";
        //DumpIntVec(cladeVec);
    }
    
    return scoreMax;
}
    
// evaluate the AG
void AGHeuristicSearchTreeBaseEval :: ProcessAG()
{
    listContainedAGTreeInfo.clear();
    
    // first extract marginal trees
    vector<GenealogicalNetworkMTreeInfo> listMargTreeInfo;
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTrees;
    this->agCurr.RetriveAllMarginalTrees(mapMargTrees);
    for(auto it = mapMargTrees.begin(); it != mapMargTrees.end(); ++it )
    {
        listMargTreeInfo.push_back(it->second);
    }
//cout << "# of marginal trees: " << listMargTreeInfo.size() << endl;
    
    // collect info
    for(unsigned int i=0; i<listMargTreeInfo.size(); ++i)
    {
//cout << "Marg tree: " <<
        string strNW = listMargTreeInfo[i].GetMargTree()->GetNewickSorted(false);
        PhylogenyTreeBasic tr;
        tr.ConsOnNewick(strNW);
        set<set<int> > setClades;
        tr.GetAllClades(setClades);
        //listMargTreeInfo[i].GetMargTree()->BuildDescendantInfo();
        //listMargTreeInfo[i].GetMargTree()->FindAllSplits( setClades);
        //listMargTreeInfo[i].GetMargTree()->ConsDecedentLeavesInfoLabels(setClades);
//cout << "Number of clades: " << setClades.size() << endl;
        set<set<int> > setCladesRem;
        for(auto cl: setClades)
        {
//cout << "Clade: ";
//DumpIntSet(setClades[i]);
            if( cl.size() > 1 )
            {
                setCladesRem.insert(cl);
            }
        }
        listContainedAGTreeInfo.push_back(setCladesRem);
//cout << "--Init: add " << i << ", setCladesRem: ";
//for(auto x: setCladesRem)
//{
//DumpIntSet(x);
//}
    }
    
    // free marg trees
    for(unsigned int i=0; i<listMargTreeInfo.size(); ++i)
    {
        delete listMargTreeInfo[i].GetMargTree();
    }
}
    

//***********************************************************************************
// construct network from trees

AGHeuristicBuilder :: AGHeuristicBuilder(vector<PhylogenyTreeBasic *> &listGeneTreePtrsIn, int ogIn) : listGeneTreePtrs(listGeneTreePtrsIn), og(ogIn), pnetCurr(NULL), scoreCurrNet(MAX_NEG_DOUBLE_VAL)
{
}

AGHeuristicBuilder :: ~AGHeuristicBuilder()
{
    delete pnetCurr;
}

void AGHeuristicBuilder :: Build(int numMixEvts)
{
    // init first
    Init();
    
    // now add reticulations
    double scoreCurr = MAX_NEG_DOUBLE_VAL;
    for(int i=0; i<numMixEvts; ++i)
    {
        scoreCurr = FindBestNetAddOneMix();
//cout << "***** After adding one reticulation, score: " << scoreCurr << endl;
    }
    // assign score
    this->scoreCurrNet = scoreCurr;
    
//cout << "***** After adding reticulation, score: " << this->scoreCurrNet << endl;
//cout << "   num of reticulate nodes: " << this->pnetCurr->GetNumMixNodes() << endl;
//this->pnetCurr->DumpMargTrees(false);
    
    // now search for optimal ngbr nets
    this->scoreCurrNet = FindOptNetLocalSearch();
//cout << "***** After local search, score: " << this->scoreCurrNet << endl;
//cout << "   num of reticulate nodes: " << this->pnetCurr->GetNumMixNodes() << endl;
//this->pnetCurr->DumpMargTrees(false);
}


// init
void AGHeuristicBuilder :: Init()
{
    // construct from NJ
    string strNW;
    
    // build a simple tree from trees heuristicly
    GNHeuSearchClusterSupport clusterSupport;
    clusterSupport.SetNumTrees(listGeneTreePtrs.size());
    clusterSupport.SetOutgroup(og);
    GNHeuSearchEvaluator gnEvalutor;
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        GNHeuSearchGTreeInfo ginfo( *listGeneTreePtrs[i] );
        gnEvalutor.Evaluate( *listGeneTreePtrs[i], ginfo, clusterSupport );
    }
//cout << "******* Found support: ";
//clusterSupport.Dump();
    strNW = clusterSupport.BuildNJTree();
//cout << "Constructed NJ tree: ";
//cout << strNW << endl;
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        GNHeuSearchGTreeInfo ginfo( *listGeneTreePtrs[i] );
        gnEvalutor.Evaluate2( *listGeneTreePtrs[i], ginfo, clusterSupport );
    }
    string strNW2 = clusterSupport.ReWeightTree(strNW);
    int numTaxa = clusterSupport.GetNumTaxa();
    // init net with this single tree
    pnetCurr = new GenealogicalNetwork;
    pnetCurr->InitWithTree(strNW2, numTaxa );
    
    double sc = ScoreNet(*this->pnetCurr);
    this->scoreCurrNet = sc;
//cout << "After init, score of curr net: " << sc << endl;
}
    
// score the curr net
double AGHeuristicBuilder :: ScoreNet(GenealogicalNetwork &net)
{
    //
    double res = 0.0;
    //AGHeuristicSearchEval agHeu(net);
    AGHeuristicSearchTreeBaseEval agHeu(net);
    for(unsigned int i=0; i<listGeneTreePtrs.size(); ++i)
    {
        double sc = agHeu.EvalTree(listGeneTreePtrs[i]);
        res += sc;
//string strNW;
//listGeneTreePtrs[i]->ConsNewickSorted(strNW);
//cout << "Tree " << i << ":" << strNW << ", score: " << sc << endl;
    }
    return res;
}

// search for one best net by making one tree branch to be reticulate and then regraft this edge to somewhere else
// for now, only consider tree edges (not those edges already reticulate)
double AGHeuristicBuilder :: FindBestNetAddOneMix()
{
    YW_ASSERT_INFO(this->pnetCurr != NULL, "Net: not initialized");
    //
    vector<GenealogicalNetworkBranch *> setAllBranches;
    this->pnetCurr->GetAllBranchesBottomUp( setAllBranches );
    double scoreBest = GetCurrScore();
    vector<GenealogicalNetwork *> listNetsOneNewMix;
    for(unsigned int i=0; i<setAllBranches.size(); ++i)
    {
        if( setAllBranches[i]->IsMixing() )
        {
            continue;
        }
        //
        GenealogicalNetworkBranch *pbr = setAllBranches[i];
        vector<GenealogicalNetwork *> listNetsOneNewMixStep;
        this->pnetCurr->RetriveNgbrNetsOneNewMixForEdge( pbr, listNetsOneNewMixStep );
        for(unsigned int kkk=0; kkk<listNetsOneNewMixStep.size(); ++kkk)
        {
            listNetsOneNewMix.push_back(listNetsOneNewMixStep[kkk]);
        }
    }
     
    // find the best net to go with
    vector<double> listScores;
    for(int j=0; j<(int)listNetsOneNewMix.size(); ++j)
    {
        double score = ScoreNet( *listNetsOneNewMix[j]);
        //double probStep = ScoreForNetworkMDC(listGeneTrees, listNetsOneNewMix[j]);
        listScores.push_back(score);
    }
//cout << "---List of scores of candidate networks: ";
//DumpDoubleVec(listScores);
    
    int maxScorePos = -1;
    if( listNetsOneNewMix.size() > 0 )
    {
        int maxScoreStep = FindMaxValPositionFromList( listScores );
//cout << "Max score pos: " << maxScoreStep << endl;
//cout << "Max score: " << listScores[maxScoreStep] << ", scoreBest: " << scoreBest << endl;
        double scoreNgbr = listScores[maxScoreStep];
        
        const double MIN_INC_SCORE = 0.000001;
        if( scoreNgbr >= scoreBest + MIN_INC_SCORE )
        {
            maxScorePos = maxScoreStep;
            scoreBest = scoreNgbr;

//cout << "maxScorePos: " << maxScorePos << ", size of net candidates: " << listNetsOneNewMix.size() << endl;
            delete this->pnetCurr;
            this->pnetCurr = listNetsOneNewMix[maxScoreStep];
        }
    }
//cout << "Now cleanup...\n";
    //
    for(int k=0; k<(int)listNetsOneNewMix.size(); ++k)
    {
        if( k != maxScorePos  )
        {
            delete listNetsOneNewMix[k];
        }
    }

    return scoreBest;
}

// Find best neighboring net by local search
double AGHeuristicBuilder :: FindOptNetLocalSearch()
{
    //
    YW_ASSERT_INFO(this->pnetCurr != NULL, "Current network is not initialized");
    
    //double probCurBest = MAX_NEG_DOUBLE_VAL;
    double scoreCurBest = GetCurrScore();
    while(true)
    {
        //#if 0
        //cout << "--Processing current network (w/ score: " << scoreCurBest << "): \n";
        //cout << "   num of reticulate nodes: " << pnetworkOpt->GetNumMixNodes() << endl;
        //pnetworkOpt->DumpMargTrees(false);
        //pnetworkOpt->Dump();
        //#endif
        //
        vector<GenealogicalNetwork *> listNgbrNetsDist1;
        vector<set<GenealogicalNetworkBranch *> > listNgbrNetsDist1CritBrs;
        int og = AGFastCoalBuilder :: GetOutgroup();
        //int og = -1;
        //cout << "og: " << og << endl;
        this->pnetCurr->RetriveNgbrNetsOneEvt2( listNgbrNetsDist1, og, &listNgbrNetsDist1CritBrs );
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
            
            // make sure it has no trivial cycles
            if( listNgbrNetsDist1[i]->CheckTrivialCycles() )
            {
                //
                continue;
            }
            
            //bool fFoundBefore = IsNetProcBefore( listNgbrNetsDist1[i] );
            //bool fOGOK = IsNetworkOGGood( *listNgbrNetsDist1[i] );
            
            bool fGOOK = true;
            if( fGOOK )
            {
                // make sure this net has the desired number of network nodes
                //if( listNgbrNetsDist1[i]->GetNumMixNodes() != this->maxNumMixNodes )
                //{
                //    continue;
                //}
                
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
        
        //cout << "Number of neighboring nets to evaluate: " << listNgbrNetsDist1New.size() << endl;
        // find ngbr networks (no adding new mix nodes)
        vector<double> listScoreNgbrNets;
        
        for(int i=0; i<(int)listNgbrNetsDist1New.size(); ++i)
        {
            YW_ASSERT_INFO( i<listNgbrNetPos.size(), "Overflow444" );
            int netIndex = listNgbrNetPos[i];
#if 0
            cout << "&&&&&&&&&&&&&&&&&&Processing ngbr network " << netIndex << endl;
            listNgbrNetsDist1[netIndex]->DumpMargTrees(false);
            listNgbrNetsDist1[netIndex]->Dump();
#endif
            //cout << "Creating length optimization...\n";
            double scoreStep = ScoreNet(*listNgbrNetsDist1[netIndex]);
            listScoreNgbrNets.push_back( scoreStep );
        }
        //cout << "**** Number of neighbr nets: " << listLogProbNgbrNets.size() << endl;
        // if nothing left, done
        int maxPos = -1;
        bool fContSearch = false;
        if( listScoreNgbrNets.size() > 0 )
        {
            int maxPosStep = FindMaxValPositionFromList( listScoreNgbrNets );
            double scoreBestNgbr = listScoreNgbrNets[maxPosStep];
            //cout << "!Best ngbr prob: " << probBestNgbr << ", probCurBest: " << probCurBest << endl;
            
            const double MIN_INC_RATIO = 1.00000001;
            if( scoreBestNgbr >= scoreCurBest + log(MIN_INC_RATIO) )
            {
                maxPos = listNgbrNetPos[maxPosStep];
                
                scoreCurBest = scoreBestNgbr;
                cout << "-- Score improved to: " << scoreCurBest << endl;
                
                delete this->pnetCurr;
                this->pnetCurr = listNgbrNetsDist1[maxPos];
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
            if( i != maxPos || fContSearch == false )
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
    
    return scoreCurBest;
}
