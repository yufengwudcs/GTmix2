//
//  GenealogicalNetworkHeuSearch.cpp
//  
//
//  Created by Yufeng Wu on 1/14/16.
//
//

#include "GenealogicalNetworkHeuSearch.h"
//#include "GenealogicalNetworkSearch.h"
#include "Utils3.h"
#include "TreeBuilder.h"
#include "GenealogicalNetwork.h"
#include "MarginalTree.h"
#include "Utils4.h"

//***********************************************************************************
// Tetsing routine

void GNHeuSearchInitNet( const vector<PhylogenyTreeBasic *> &listGenetreePtrs, GenealogicalNetwork &netInit )
{
    GNHeuSearchClusterSupport clusterSupport;
    
    GNHeuSearchEvaluator gnEvalutor;
    for(int i=0; i<(int)listGenetreePtrs.size(); ++i)
    {
        GNHeuSearchGTreeInfo ginfo( *listGenetreePtrs[i] );
        gnEvalutor.Evaluate( *listGenetreePtrs[i], ginfo, clusterSupport );
    }
    
//cout << "******* Found support: ";
//clusterSupport.Dump();
    
//cout << "Constructed NJ tree: ";
    string strNW = clusterSupport.BuildNJTree();
//cout << strNW << endl;
    
    // test network
    int numTaxa = clusterSupport.GetNumTaxa();
//cout << "numTaxa: " << numTaxa << endl;
    //GenealogicalNetwork net1;
    netInit.InitWithTree(strNW, numTaxa );
//cout << "Constructed network: ";
//net1.Dump();
}


//***********************************************************************************
// Basic info about gene tree


GNHeuSearchGTreeInfo :: GNHeuSearchGTreeInfo( PhylogenyTreeBasic &gtree )
{
    Init(gtree);
}
void GNHeuSearchGTreeInfo :: Init( PhylogenyTreeBasic &gtree )
{
    // collect num of alleles for each taxon
    vector< TreeNode *> listLeafNodes;
    gtree.GetAllLeafNodes(listLeafNodes);
    for(int i=0; i<(int)listLeafNodes.size(); ++i )
    {
        int lbl = listLeafNodes[i]->GetIntLabel();
        if( mapTaxonAlleleCounts.find(lbl) == mapTaxonAlleleCounts.end() )
        {
            mapTaxonAlleleCounts.insert( map<int,int> :: value_type(lbl, 0) );
        }
        ++mapTaxonAlleleCounts[lbl];
    }
    
    // now record for all nodes
    vector< TreeNode *> listAllNodes;
    gtree.GetAllNodes(listAllNodes);
    for(int i=0;i<(int)listAllNodes.size(); ++i)
    {
        set<TreeNode *> setDescendents;
        listAllNodes[i]->GetAllLeavesUnder(setDescendents);
        multiset<int> ss;
        for( set<TreeNode *> :: iterator it = setDescendents.begin(); it != setDescendents.end(); ++it )
        {
            int lbl = (*it)->GetIntLabel();
            ss.insert(lbl);
        }
        mapNodeLeafIds.insert( map<TreeNode *, multiset<int> > :: value_type(listAllNodes[i], ss)) ;
    }
    // collect taxa and leaves
    setTaxa.clear();
    gtree.GetRoot()->GetAllDescendIntLbls( setTaxa );
    
    for(set<int> :: iterator it = setTaxa.begin(); it != setTaxa.end(); ++it)
    {
        char buf[1000];
        snprintf(buf, 1000, "%d", *it);
        string strBuf = buf;
        set<string> setLabels;
        setLabels.insert(strBuf);
        set<TreeNode *> setLvNodes;
        gtree.GetLeavesWithLabels( setLabels, setLvNodes );
    
        mapTaxonToLeaves.insert( map<int, set<TreeNode *> > :: value_type(*it, setLvNodes) );
    }
    
    // collect distance info
    for(int i=0; i<(int)listLeafNodes.size(); ++i)
    {
        for(int j=i+1; j<(int)listLeafNodes.size(); ++j)
        {
            TreeNode *pMRCA = listLeafNodes[i]->GetMRCA(listLeafNodes[j]);
            int ne1 = listLeafNodes[i]->GetNumEdgesToAncestor( pMRCA );
            int ne2 = listLeafNodes[j]->GetNumEdgesToAncestor( pMRCA );
            int ne = ne1+ne2;
            pair<TreeNode *, TreeNode *> pp1(listLeafNodes[i], listLeafNodes[j]), pp2(listLeafNodes[j], listLeafNodes[i]);
            mapEdgeNumBetweenPairLeaves.insert( map< pair<TreeNode *, TreeNode *>, int > :: value_type(  pp1, ne ) );
            mapEdgeNumBetweenPairLeaves.insert( map< pair<TreeNode *, TreeNode *>, int > :: value_type(  pp2, ne ) );
        }
    }
}
int GNHeuSearchGTreeInfo :: GetNumOfAlleles(int taxon) const
{
    YW_ASSERT_INFO( mapTaxonAlleleCounts.find(taxon) != mapTaxonAlleleCounts.end(), "Fail to find" );
    return mapTaxonAlleleCounts.find(taxon)->second;
}
int GNHeuSearchGTreeInfo :: GetNumOfAllelesFor(const set<int> &taxa) const
{
    int res = 0;
    for(set<int> :: iterator it = taxa.begin(); it !=taxa.end(); ++it)
    {
        res += GetNumOfAlleles(*it);
    }
    return res;
}
int GNHeuSearchGTreeInfo :: GetNumOfLeavesUnder(TreeNode *pncurr) const
{
    GNHeuSearchGTreeInfo *pthis = const_cast<GNHeuSearchGTreeInfo *>(this);
    YW_ASSERT_INFO( mapNodeLeafIds.find(pncurr) != mapNodeLeafIds.end(), "Fail to find" );
    return pthis->mapNodeLeafIds[pncurr].size();
}
set<TreeNode *> & GNHeuSearchGTreeInfo :: GetSetLeavesFor(int taxon)
{
    //
    YW_ASSERT_INFO( mapTaxonToLeaves.find(taxon) != mapTaxonToLeaves.end(), "Fail to find" );
    return mapTaxonToLeaves[taxon];
}

int GNHeuSearchGTreeInfo :: GetNumEdgesBtwLeaves( TreeNode *pLeaf1, TreeNode *pLeaf2 )
{
    //
    pair<TreeNode *,TreeNode *> pp(pLeaf1, pLeaf2);
    YW_ASSERT_INFO( mapEdgeNumBetweenPairLeaves.find(pp) != mapEdgeNumBetweenPairLeaves.end(), "Fail to find");
    return mapEdgeNumBetweenPairLeaves[pp];
}

//***********************************************************************************
// How much support a cluster receives


void GNHeuSearchClusterSupport :: AddSupport( const set<int> &setTaxa, double val )
{
    //
    if( mapClusterSupport.find(setTaxa) == mapClusterSupport.end() )
    {
        //
        mapClusterSupport.insert( map<set<int>, double> :: value_type( setTaxa, 0.0 )  );
    }
    mapClusterSupport[setTaxa] += val;
}

void GNHeuSearchClusterSupport :: AddSupport2( const set<int> &setTaxa, double val )
{
    mapClusterSupport2[setTaxa] += val;
}

void GNHeuSearchClusterSupport :: AddDistEst(int taxon1, int taxon2, double distEst)
{
    //
    pair<int,int> pp(taxon1,taxon2);
    if( mapTaxaPairDist.find( pp ) == mapTaxaPairDist.end() )
    {
        //
        mapTaxaPairDist.insert( map< pair<int,int>, double > :: value_type(pp, 0.0) );
        mapTaxaPairDistNumEst.insert( map< pair<int,int>, int > :: value_type(pp, 0) );
    }
    mapTaxaPairDist[pp] += distEst;
    ++mapTaxaPairDistNumEst[pp];
}

void GNHeuSearchClusterSupport :: Dump() const
{
    GNHeuSearchClusterSupport *pthis=const_cast<GNHeuSearchClusterSupport *>(this);
    //
    for( map<set<int>,double> :: const_iterator it = mapClusterSupport.begin(); it != mapClusterSupport.end(); ++it )
    {
        cout << "[" << it->second << "]: "  << "Cluster:";
        DumpIntSet(it->first);
    }
    
    cout << "Distance estimate:\n";
    for( map< pair<int,int>, double > :: const_iterator it= mapTaxaPairDist.begin(); it != mapTaxaPairDist.end(); ++it )
    {
        cout << "<num-of-est:" <<  pthis->mapTaxaPairDistNumEst[it->first] << "> [" << it->first.first << ", " << it->first.second << "]: dist=" << it->second << endl;
    }
}

string GNHeuSearchClusterSupport :: BuildNJTree() const
{
    // build NJ tree based on pairwise distance
    GNHeuSearchClusterSupport *pthis=const_cast<GNHeuSearchClusterSupport *>(this);
    PhyloDistance phyDist;
    for( map< pair<int,int>, double > :: const_iterator it= mapTaxaPairDist.begin(); it != mapTaxaPairDist.end(); ++it )
    {
        double dist = it->second/( pthis->mapTaxaPairDistNumEst[ it->first ] );
        phyDist.SetDistance( it->first.first, it->first.second, dist );
    }
//cout << "GNHeuSearchClusterSupport :: BuildNJTree(): phyDist: ";
//phyDist.Dump();
    DistanceTreeBuilder distNJ(phyDist);
    //distNJ.SetOutgroup( GenealogicalNetworkSearch::GetOutgroup() );
    distNJ.SetOutgroup( GetOutgroup() );
    return distNJ.NJ();
}

string GNHeuSearchClusterSupport :: ReWeightTree(const string &strTree) const
{
    // weight tree branch length by the frequency of clades
    PhylogenyTreeBasic phTree;
    phTree.ConsOnNewick(strTree);
//cout << "ReWeightTree: tree is: ";
//string strNW;
//phTree.ConsNewick(strNW);
//cout << strNW << endl;
    MarginalTree mtree;
    ReadinMarginalTreesNewickWLenString(strTree, phTree.GetNumLeaves(), mtree );
//cout << "mtree: " << mtree.GetNewick() << endl;
    
    vector<set<int> > leafNodeLabels;
    mtree.ConsDecedentLeavesInfoLabels( leafNodeLabels );
    
    //
    for(unsigned int i=0; i<leafNodeLabels.size(); ++i)
    {
        set<int> setIntLbs = leafNodeLabels[i];
        
        auto it = mapClusterSupport2.find(setIntLbs);
        const double MIN_LEN = 0.001;
        double len = MIN_LEN;
        if( it != mapClusterSupport2.end() )
        {
            len = it->second/numTrees;
        }
        if( len < MIN_LEN )
        {
            len = MIN_LEN;
        }
        mtree.SetBranchLen(i, len);
    }

    string strNW1 = mtree.GetNewick();
//cout << "AFTER ReWeightTree: tree is: ";
//cout << strNW1 << endl;
    return strNW1;
}

/*
string GNHeuSearchClusterSupport :: BuildTree2() const
{
    // get max wt
    YW_ASSERT_INFO( mapClusterSupport2.size() > 0, "Fatal error: cluster support2 is not initialized" );
    double maxWt = 0.0;
    for(auto it = mapClusterSupport2.begin(); it != mapClusterSupport2.end(); ++it)
    {
        if( it->second > maxWt )
        {
            maxWt = it->second;
        }
    }
    
    // build tree by prioritizing clade frequencies
    cout << "GNHeuSearchClusterSupport :: BuildTree2 :: Clade frequncies: \n";
    for(auto it = mapClusterSupport2.begin(); it != mapClusterSupport2.end(); ++it)
    {
        cout << "[" << (double)it->second/maxWt << "]: ";
        DumpIntSet(it->first);
    }
    // tbd;
    string strTree;
    
    return strTree;
}
 */

int GNHeuSearchClusterSupport :: GetNumTaxa() const
{
    //
    set<int> ss;
    for( map< pair<int,int>, double > :: const_iterator it= mapTaxaPairDist.begin(); it != mapTaxaPairDist.end(); ++it )
    {
        ss.insert( it->first.first );
        ss.insert( it->first.second );
    }
    return ss.size();
}

//***********************************************************************************
// Evaluate the tree


void  GNHeuSearchEvaluator ::  Evaluate( PhylogenyTreeBasic &treeGene, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport )
{
    // examine each internal node (coalescent) to obtain support for clusters
    PhylogenyTreeIterator itor(treeGene);
    itor.Init();
    
    while( itor.IsDone() == false )
    {
        TreeNode *pncurr = itor.GetCurrNode();
        
        if( pncurr->IsLeaf() == false )
        {
            // test it
            EvaluateNode( pncurr, treeGeneInfo, clusterSupport );
        }
        
        itor.Next();
    }
    
    // estimate pairwise distance
    set<int> &setTaxa = treeGeneInfo.GetAllTaxa();
    for( set<int> :: iterator it1 = setTaxa.begin(); it1 != setTaxa.end(); ++it1 )
    {
        //
        set<int> :: iterator it2 = it1;
        ++it2;
        for(; it2 != setTaxa.end(); ++it2)
        {
            EvaluateTaxaDist( *it1, *it2, treeGeneInfo, clusterSupport );
        }
    }
}

void  GNHeuSearchEvaluator ::  Evaluate2( PhylogenyTreeBasic &treeGene, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport )
{
    /*
cout << "Evaluate2: evaluting tree: ";
string strNW;
treeGene.ConsNewick(strNW);
cout << strNW << endl;
    // get allele frequency
    map<int,double> mapFreqs;
    PhylogenyTreeIterator itor2(treeGene);
    itor2.Init();
    
    while( itor2.IsDone() == false )
    {
        TreeNode *pncurr = itor2.GetCurrNode();
        
        if( pncurr->IsLeaf() )
        {
            mapFreqs[pncurr->GetIntLabel()] += 1.0/treeGene.GetNumLeaves();
        }
        
        itor2.Next();
    }
cout << "Allele frequency: \n";
for(auto x : mapFreqs)
{
cout << "Allele: " << x.first << ": " << x.second << endl;
}
 */
    
    // examine each node to obtain support for clusters
    PhylogenyTreeIterator itor(treeGene);
    itor.Init();
    
    set<set<int> > setClades;
    
    while( itor.IsDone() == false )
    {
        TreeNode *pncurr = itor.GetCurrNode();
        
        // extract the subsets
        //EvaluateNode2( pncurr, treeGeneInfo, clusterSupport, mapFreqs );
        set<int> setIntLbs;
        pncurr->GetAllDescendIntLbls( setIntLbs );
        setClades.insert(setIntLbs);
        
        itor.Next();
    }
    
    // add up
    for(auto &x : setClades)
    {
        clusterSupport.AddSupport2( x, 1.0);
    }
    
}

void GNHeuSearchEvaluator :: EvaluateNode( TreeNode *pncurr, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport )
{
    // process the node to find support of clusters
    set<int> setIntLbs;
    pncurr->GetAllDescendIntLbls( setIntLbs );
    
    // add a record by one
    int numTotAlleles = treeGeneInfo.GetNumOfAllelesFor( setIntLbs );
    int numAllelesUnder = treeGeneInfo.GetNumOfLeavesUnder(  pncurr);
    double val = (1.0*numAllelesUnder)/numTotAlleles;
    clusterSupport.AddSupport( setIntLbs, val);
}

/*
void GNHeuSearchEvaluator :: EvaluateNode2( TreeNode *pncurr, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport, const map<int,double> &mapFreqs )
{
    // process the node to find support of clusters
    set<int> setIntLbs;
    pncurr->GetAllDescendIntLbls( setIntLbs );
    
    // add support2 but skip clades that are purely within one population
    if( setIntLbs.size() > 1 || pncurr->IsLeaf() )
    {
        double val = 1.0;
        // weight by 1 over product of frequcnies
        for(auto x: setIntLbs)
        {
            auto it = mapFreqs.find(x);
            YW_ASSERT_INFO(it != mapFreqs.end(), "Fail to find4445");
            val /= it->second;
        }
        
        clusterSupport.AddSupport2( setIntLbs, val);
cout << "EvaluateNode2: val:" << val << " for clade: ";
DumpIntSet(setIntLbs);
    }
}
*/

void GNHeuSearchEvaluator :: EvaluateTaxaDist( int taxon1, int taxon2, GNHeuSearchGTreeInfo &treeGeneInfo, GNHeuSearchClusterSupport &clusterSupport  )
{
    //
    set<TreeNode *> &setLeaves1 = treeGeneInfo.GetSetLeavesFor(taxon1);
    set<TreeNode *> &setLeaves2 = treeGeneInfo.GetSetLeavesFor(taxon2);
    // take average
    int totDist = 0;
    for( set<TreeNode *> :: iterator it1 = setLeaves1.begin(); it1 != setLeaves1.end(); ++it1 )
    {
        for( set<TreeNode *> :: iterator it2 = setLeaves2.begin(); it2 != setLeaves2.end(); ++it2)
        {
            int ne = treeGeneInfo.GetNumEdgesBtwLeaves( *it1, *it2 );
            totDist += ne;
        }
    }
    int numPairs = setLeaves1.size() * setLeaves2.size();
    YW_ASSERT_INFO(numPairs > 0, "Can not vanish");
    double dist = ((double)totDist)/numPairs;
    clusterSupport.AddDistEst( taxon1, taxon2, dist );
}

//***********************************************************************************
// Search network based on subtree cluster scoring
// i.e. prefer networks with marginal tree that provides more support from subnet structures


GNHeuInfSubtreeCluster :: GNHeuInfSubtreeCluster(GenealogicalNetwork &netInitIn, vector<PhylogenyTreeBasic *> &listGenetreePtrsIn, TaxaMapper *pMapperTaxonIdIn, map< set<int>,double > &mapClusterSupportIn) : netInit(netInitIn), listGenetreePtrs(listGenetreePtrsIn), pMapperTaxonId(pMapperTaxonIdIn), mapClusterSupport(mapClusterSupportIn), pnetworkOpt(NULL)
{
    //
}

GNHeuInfSubtreeCluster::~GNHeuInfSubtreeCluster()
{
    if( pnetworkOpt != NULL )
    {
        delete pnetworkOpt;
    }
}

void GNHeuInfSubtreeCluster :: Infer()
{
    // starting from the init net and perform neighbor search
    // initialize if needed
    if( pnetworkOpt == NULL )
    {
        pnetworkOpt = netInit.Copy();
        // TBD??????
        pnetworkOpt->MapBackUserLabels(*pMapperTaxonId);
    }
    //
    //double scoreCurBest = ScoreNetwork(pnetworkOpt);
    
    // for now allow one network node so explore one time
    double scoreCurBest = SearchWithNewMix();
    
//cout << "INITIAL network score: " << scoreCurBest << endl;
    
    while(true)
    {
        //
        vector<GenealogicalNetwork *> listNgbrNetsDist1;
        pnetworkOpt->RetriveNgbrNetsOneEvt( listNgbrNetsDist1 );
        
        // find ngbr networks (no adding new mix nodes)
        vector<double> listScoresNgbrNets;
        for(int i=0; i<(int)listNgbrNetsDist1.size(); ++i)
        {
            double scoreStep = ScoreNetwork( listNgbrNetsDist1[i] );
            listScoresNgbrNets.push_back( scoreStep );
//cout << "Network ngbr " << i << ": score: " << scoreStep << endl;
        }
        
        // if nothing left, done
        int maxPos = -1;
        bool fContSearch = false;
        maxPos = FindMaxValPositionFromList( listScoresNgbrNets );
        double resStepBest = listScoresNgbrNets[maxPos];
        
        if( resStepBest > scoreCurBest)
        {
            // improvement found
            scoreCurBest = resStepBest;
            delete pnetworkOpt;
            pnetworkOpt = listNgbrNetsDist1[maxPos];
            fContSearch = true;
            //cout << "UPDATING THE CURRENT NETWORK TO: ";
            //pnetworkOpt->Dump();
        }
        
        for(int i=0; i<(int)listNgbrNetsDist1.size(); ++i)
        {
            if( i != maxPos || fContSearch == false )
            {
                delete listNgbrNetsDist1[i];
            }
        }
        //cout << "@@@@@@@@@@@@@@@ Current best log-likelihood: " << probCurBest << endl;
        
        if( fContSearch == false )
        {
            break;
        }
    }
//cout << "Finally best score found:" << scoreCurBest << endl;
}

double GNHeuInfSubtreeCluster :: SearchWithNewMix()
{
    // search from current best network and adding one new mixing node
    YW_ASSERT_INFO( pnetworkOpt!=NULL, "Wrong: must be initialized" );
    
    double scoreBestCurrFromPrior = ScoreNetwork(pnetworkOpt);
    double scoreCurBest = scoreBestCurrFromPrior;
    //cout << "SearchWithNewMix:probCurBest (at entry)" << probCurBest << endl;
    //cout << "CURRENT NETWORK: ";
    //pnetworkOpt->Dump();
    //
    vector<GenealogicalNetwork *> listNgbrNetsDist1;
    YW_ASSERT_INFO(pnetworkOpt != NULL, "Must have already been initialized");
    pnetworkOpt->RetriveNgbrNetsOneNewMix( listNgbrNetsDist1 );
    
    // find ngbr networks (no adding new mix nodes)
    vector<double> listScoreNgbrNets;
    for(int i=0; i<(int)listNgbrNetsDist1.size(); ++i)
    {
        //cout << "############################Processing network: ";
        //listNgbrNetsDist1[i]->Dump();
        // here, don't optimize over branch; rather, simply search for the single best network
        //double logprobBest = opt.CalcProb();
        double scoreStep = ScoreNetwork(listNgbrNetsDist1[i]);        // optimize over branch length
        listScoreNgbrNets.push_back( scoreStep );
        //cout << "Network ngbr " << i << ": optimized likelihood: " << logprobBest << endl;
    }
    
    int maxPos = -1;
    if( listScoreNgbrNets.size() > 0 )
    {
        maxPos = FindMaxValPositionFromList( listScoreNgbrNets );
        double scoreBestNgbr = listScoreNgbrNets[maxPos];
        if( scoreBestNgbr > scoreCurBest )
        {
            // improvement found
            scoreCurBest = scoreBestNgbr;
            delete pnetworkOpt;
            pnetworkOpt = listNgbrNetsDist1[maxPos];
            //cout << "UPDATING THE CURRENT NETWORK TO: ";
            //pnetworkOpt->Dump();
        }
    }
    
    for(int i=0; i<(int)listNgbrNetsDist1.size(); ++i)
    {
        if( i != maxPos )
        {
            delete listNgbrNetsDist1[i];
        }
    }
    //cout << "@@@@@@@@@@@@@@@ Current best log-likelihood: " << probCurBest << endl;
    
    
    return scoreCurBest;
}


double GNHeuInfSubtreeCluster :: ScoreNetwork( GenealogicalNetwork *pNetToEval )
{
    // score the current network: extract all marginal trees
    map<string,string> mapTaxonIdToUserId;
    pNetToEval->CreateMapTaxonIdToUserId( mapTaxonIdToUserId );
//cout << "network: ";
//pNetToEval->Dump();
//cout << "mapTaxonIdToUserId: ";
//for(map<string,string> :: iterator itt = mapTaxonIdToUserId.begin(); itt != mapTaxonIdToUserId.end(); ++itt)
//{
//cout << itt->first << " to " << itt->second << endl;
//}
    // get all the marg trees
    map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> mapMargTreesWithFreq;
    pNetToEval->RetriveAllMarginalTrees( mapMargTreesWithFreq );
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo>::const_iterator it= mapMargTreesWithFreq.begin(); it!=mapMargTreesWithFreq.end(); ++it)
    {
//cout <<"BEFORE remapping: tree: ";
//cout << it->second.GetMargTree()->GetNewick() << endl;
        // should remap id first
        RemapLeafIntLabelsTaxaMap(*(it->second.GetMargTree()), mapTaxonIdToUserId);
//cout <<"After remapping: tree: ";
//cout << it->second.GetMargTree()->GetNewick() << endl;
    }
    
    double res = 0.0;
    for( map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo>::const_iterator it= mapMargTreesWithFreq.begin(); it!=mapMargTreesWithFreq.end(); ++it)
    {
        // should remap id first
        double resstep = ScoreMargTree( it->second.GetMargTree() );
        res += resstep* it->second.GetFreq();
        it->second.FreeTree();
    }
    return res;
}

// Score an individual marg tree
double GNHeuInfSubtreeCluster :: ScoreMargTree(MarginalTree *ptree)
{
//cout << "Scoring margintree: " << ptree->GetNewick() << endl;
    // format: [freq] newick-format-tree
    vector< set<int> > listSplits;
    ptree->FindAllSplits(listSplits);
    
    double res = 0.0;
    for(int i=0; i<(int)listSplits.size(); ++i)
    {
        double ressplit = ScoreSplit( listSplits[i] );
//cout << "ScoreMargTree: score=" << ressplit << ": split: ";
//DumpIntSet( listSplits[i]);
        res+= ressplit;
    }

    return res;
}

double GNHeuInfSubtreeCluster :: ScoreSplit(const set<int> &split)
{
    //
    double resstep = 0.0;
    if( mapClusterSupport.find( split ) != mapClusterSupport.end() )
    {
        resstep = mapClusterSupport[ split ];
    }
    return resstep;
}

