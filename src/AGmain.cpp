//
//  
//  
//
//  Modified by Yufeng Wu on 3/22/24.
//
//

#include "AGGeneTreeQuickCoal.hpp"
#include "GenealogicalNetwork.h"
#include "PhylogenyTreeBasic.h"
#include "AGGeneTreeProb.hpp"
#include "AGFastLenOpt.hpp"
#include "AGFastCoalBuilder.hpp"
#include "ApproxGeneTreeProb2.h"
#include "Utils3.h"
#include <string>
#include <sstream>
#include "GeneTreeProcUtils.hpp"
#include "GenealogicalNetworkClades.hpp"
#include "GenealogicalNetworkHeuSearch.h"
#include "AGHeuristicSearch.hpp"
//#include "GML2ExtNewick.h"

// definitions
static char fileNetworkGMLDefault[100]="optimal-network.gml";

// settings
static char *fileNetworkGMLOut = fileNetworkGMLDefault;
static char *filenameGenetrees = NULL;
//static char *filenameHaps = NULL;
static char *filenamePop = NULL;
//static int szSubsetTaxa = -1;
static int szCompositeTaxa = -1;
//static int szMaxSubtreeIdent = -1;
//static int szMaxSubtreeSz = -1;
//static int fSummaryMode = false;
static int numThreads = 1;
static int numAdmixNodes = 1;
static int maxNumAdmixNode = 1;
static double costAddOneAdmix = 50.0;
static bool fSlowSearch = false;
static string taxonOutgroup = "";
static set<string> setMixTaxa;      // user-speficied mix taxa (for initial graph construction only for now)
static int numTreesToPick = 10000;   // no more than 10,000 trees for inference, for the sake of efficiency
static bool fHeuSearch = false;     // for large data
static bool fInitNetOnly = false;   // only build initial net; no optimization
static string netFileIn;

int GetNumThreads()
{
    return numThreads;
}


// Local functions
static bool CheckArguments(int argc, char **argv)
{
    if( argc <= 1  )
    {
        return false;
    }
    
    // Check argument one by one
    int argpos = 1;
    while( argpos < argc)
    {
        //
        if( argv[ argpos ][ 0 ] != '-' )
        {
            // this is the gene tree file
            filenameGenetrees = argv[argpos];
            argpos++;
        }
        else if( argv[argpos][1] == 'P' && argpos+1< argc )
        {
            // setup filename
            //argpos ++;
            //filenameHaps = argv[argpos];
            argpos++;
            filenamePop = argv[argpos];
            argpos++;
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 'o' )
        {
            // setup filename
            argpos ++;
            fileNetworkGMLOut = argv[argpos];
            argpos++;
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 't' )
        {
            // setup filename
            argpos ++;
            int numThreadsProb = numThreads;
            sscanf( argv[argpos], "%d", &numThreadsProb);
            cout << "Set number of threads to: " << numThreadsProb << endl;
            numThreads = numThreadsProb;
            //GenealogicalNetworkSearch::SetNumThreads(numThreads);
            //GenealogicalNetworkProb::SetNumThreads(numThreads);
            ApproxGTPCache2 :: Instance().InitNumThreads(numThreads);
            argpos++;
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 'm' )
        {
            // setup admixture upper bound
            argpos ++;
            sscanf( argv[argpos], "%d", &maxNumAdmixNode);
            cout << "Set maximum number of admixture nodes to: " << maxNumAdmixNode << endl;
            argpos++;
            float costHere;
            sscanf( argv[argpos], "%f", &costHere);
            costAddOneAdmix = costHere;
            cout << "Set admixture cost to: " << costAddOneAdmix << endl;
            argpos++;
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 'n' )
        {
            argpos ++;
            sscanf( argv[argpos], "%d", &numAdmixNodes);
            cout << "Fix number of admixture nodes to: " << numAdmixNodes << endl;
            if( maxNumAdmixNode < numAdmixNodes )
            {
                maxNumAdmixNode = numAdmixNodes;
            }
            argpos++;
        }
        else if( argv[argpos][1] == 'T' )
        {
            // setup filename
            argpos ++;
            sscanf( argv[argpos], "%d", &numTreesToPick);
            //cout << "Set the number of trees to use: " << numTreesToPick << endl;
            argpos++;
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 'r' )
        {
            // setup filename
            argpos ++;
            //char buf[100];
            //sscanf( argv[argpos], "%", &taxonOutgroup);
            taxonOutgroup = argv[argpos];
            cout << "Set the outgroup to: " << taxonOutgroup << endl;
            argpos++;
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 'i' )
        {
            // setup filename
            argpos ++;
            //char buf[100];
            //sscanf( argv[argpos], "%", &taxonOutgroup);
            netFileIn = argv[argpos];
            cout << "Read network file (in TreeMix format) from: " << netFileIn << endl;
            argpos++;
        }
        else if( argv[argpos][1] == 'q' )
        {
            // only build initial network
            fInitNetOnly = true;
            cout << "Only build initial network; no optimization " << endl;
            argpos++;
        }
        else if( argv[argpos][1] == 'a' )
        {
            argpos++;
            string ogs = argv[argpos];
            cout << "Adding initial outgroup: " << ogs << endl;
            setMixTaxa.insert(ogs);
            argpos++;
        }
        
        else
        {
            return false;
        }
    }
    
    
    return true;
}

static void ReadPopulationFile( const char *filenamePopIn, vector<pair<string,vector<int> > > &listPopRowsIn )
{
    // read population file: format: each row refers to one population: <population name> <num haps> <list of haplotype rows>
    //cout << "pass 1\n";
    listPopRowsIn.clear();
    ifstream inFileGT( filenamePopIn );
    if(!inFileGT)
    {
        cout << "Can not open population file: "<< filenamePopIn <<endl;
        exit(1);
    }
    while( inFileGT.eof() == false  )
    {
        // read one line
        const int SZ_BUF=10240;
        char buf[SZ_BUF];
        inFileGT.getline(buf, SZ_BUF);
        string strBuf(buf);
        if( strBuf.size() == 0 )
        {
            break;
        }
        std::stringstream ss(strBuf);
        //
        string strPopName;
        int nind;
        ss >> strPopName >> nind;
//cout << "Population: " << strPopName << ", num of haplotypes: " << nind;
        vector<int> listRows;
        for(int i=0; i<nind; ++i)
        {
            int row;
            ss >> row;
            listRows.push_back(row);
            //listRows.push_back(row-1);
//cout << " " << row;
        }
//cout << endl;
        pair<string, vector<int> > pp(strPopName, listRows);
        listPopRowsIn.push_back( pp );
    }

    inFileGT.close();
}

static void InitProbs()
{
return;
    /*
    const int maxCfgNum = 100;
    ApproxLineageConfigSore :: SetMaxLinCfgNum(maxCfgNum);
     */
}

static void TestNetInference( const char *filenameGenetrees, TaxaMapper *pmapperTaxaIdsUpperLevel = NULL )
{
#if 0
if( pmapperTaxaIdsUpperLevel != NULL )
{
cout << "TestNetInference: pmapperTaxaIdsUpperLevel: ";
pmapperTaxaIdsUpperLevel->Dump();
}
#endif
    YW_ASSERT_INFO( filenamePop != NULL, "Population file is not provided" );
    
    // init
    InitProbs();
    
    // read in haplotypes if provided
    //MultiBinaryMatrices listHapMats;
    vector<pair<string,vector<int> > > listPopRows;
    //if( filenameHaps != NULL )
    //{
    //    listHapMats.ReadFromFile( filenameHaps );
    //cout << "Number of loci: " << listHapMats.GetNumMats() << endl;
    YW_ASSERT_INFO( filenamePop != NULL, "Population file is not provided" );
    ReadPopulationFile(filenamePop, listPopRows);
    //}
    // construct population mapping
    TaxaMapper mapperPops;
    for(int i=0; i<(int)listPopRows.size(); ++i)
    {
        mapperPops.AddTaxaString(listPopRows[i].first);
    }
    
    //
    // now read in the set of gene trees
	vector<PhylogenyTreeBasic *> listGeneTreePtrs;
    //cout << "pass 1\n";
	ifstream inFileGT( filenameGenetrees );
	if(!inFileGT)
	{
		cout << "Can not open gene tree: "<< filenameGenetrees <<endl;
		exit(1);
	}
	TaxaMapper mapperTaxaIds;
	// read in phylogenetic tree
	if( ReadInGeneTrees( inFileGT, listGeneTreePtrs, &mapperTaxaIds ) == false )
	{
		cout << "CAUTION: checking the input file. It seems some part of input may be wrong\n";
	}
	inFileGT.close();
	YW_ASSERT_INFO( listGeneTreePtrs.size() > 0, "Must have at least one tree" );

#if 0
// test code
AdjustLabelsTaxaMapper(listGeneTreePtrs, &mapperTaxaIds);
    // dec leaf label by one
    map<int,int> mapDecBy1;
    for(int i=0; i<listGeneTreePtrs[0]->GetNumLeaves(); ++i)
    {
        mapDecBy1[i+1] = i;
    }
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        ChangeLeafIntLabelOfTree( *listGeneTreePtrs[i], mapDecBy1, true );
    }
#endif
    
#if 0
cout << "Trees after conversion..\n";
for(int tr=0; tr<(int)listGeneTreePtrs.size(); ++tr)
{
string strNewick1;
listGeneTreePtrs[tr]->ConsNewick(strNewick1);
cout << "Constructed one gene Tree = " << strNewick1 << endl;
}
//TestGNClades( listPopRows, listGeneTreePtrs );
#endif
    
    // format conversion of gene trees
    AdjustLabelsForPopInfo(listGeneTreePtrs, listPopRows, &mapperTaxaIds);
    
    int numTreesOrig = listGeneTreePtrs.size();
    SetNumTreesToPick(numTreesToPick);
    PickTrees( listGeneTreePtrs );
cout << "Number of gene trees to use: " << listGeneTreePtrs.size() << endl;
    if( (int)listGeneTreePtrs.size() < numTreesOrig )
    {
        cout << "Warning: not all input gene trees are used for inference.\n";
    }
#if 0
	for(int tr=0; tr<(int)listGeneTreePtrs.size(); ++tr)
	{
		string strNewick1;
		listGeneTreePtrs[tr]->ConsNewick(strNewick1);
cout << "Constructed one gene Tree = " << strNewick1 << endl;
	}
cout << "TestNetInference: Taxamapper: ";
mapperTaxaIds.Dump();
//exit(1);
#endif
    
    // setup outgroup if needed
    int og = -1;
    if( taxonOutgroup.length() > 0 )
    {
        og = mapperPops.GetId( taxonOutgroup );
        //int og = mapperTaxaIds.GetId( taxonOutgroup );
        //GenealogicalNetworkSearch::SetOutgroup(og);
        AGFastCoalBuilder::SetOutgroup(og);
    }
    
    // setup initial mix taxa if needed
    set<int> setMixTaxaInit;
    for(set<string> :: iterator it = setMixTaxa.begin(); it != setMixTaxa.end(); ++it)
    {
        int tid = mapperPops.GetId( *it );
        if( tid >= 0 )
        {
            setMixTaxaInit.insert(tid);
        }
        else
        {
            cout << "Warning: the designated admixed taxon " << *it << " is invalid!\n";
        }
    }
    if( setMixTaxaInit.size() > numAdmixNodes )
    {
        cout << "You specified too many initial admixed taxa. \n";
        exit(1);
    }
    
	//vector<int> listMultiplicity;
	//for(int tt=0; tt<(int)listGeneTreePtrs.size(); ++tt)
	//{
	//	listMultiplicity.push_back(1);
	//}
    
    // construct initial network; for now, just 
	YW_ASSERT_INFO( listGeneTreePtrs.size() > 0, "Tree list is empty" );
    
    // test code
    //AGHeuristicBuilder agHeuBuilder(listGeneTreePtrs, og);
    //agHeuBuilder.Build( numAdmixNodes );
    //exit(1);
    
    
    string strNW;
    int numTaxa = 0;
    
    //if( filenameHaps == NULL )
    //{
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
//cout << "After rewighting tree: ";
//cout << strNW2;
//exit(1);
    // test network
    numTaxa = clusterSupport.GetNumTaxa();
    //}
    //else
    //{
        // use haplotypes to initiialze the network
    //    numTaxa = listPopRows.size();
        //
    //    ScaffoldTreeHapBuillder treeBuilder( listHapMats, listPopRows  );
    //    strNW = treeBuilder.Build();
    //}
    YW_ASSERT_INFO(numTaxa > 0, "Number of population is not set");
    cout << "Scaffold tree: " << mapperPops.ConvIdStringWithOrigTaxa( strNW2 ) << endl;
    
//cout << "numTaxa: " << numTaxa << endl;

#if 0
cout << "Trees in initial network:\n";
net1.DumpMargTrees();
//cout << "Subtrees(init): \n";
//const int szSubsetTaxa=3;
//map<set<int>, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> > mapSubMargTreesWithFreq;
//net1.RetriveAllSubMarginalTrees(szSubsetTaxa, mapSubMargTreesWithFreq );
#endif
    
    long tmLast = GetCurrentTimeTick();
    
    GenealogicalNetwork net1;
    //net1.InitWithTree(strNW, numTaxa );
    net1.InitWithTree(strNW2, numTaxa );
//cout << "Constructed network: ";
//net1.Dump();
    
    // even no network file is given, we still want to initialize cacche, so just do this here
    AGFastCoalBuilder optNet2(listGeneTreePtrs, net1, mapperPops);
    optNet2.SetInitNetOnly(fInitNetOnly);
    optNet2.SetInitMixTaxa(setMixTaxaInit);
    
    if( netFileIn.length() == 0 )
    {
        double logprobBest = optNet2.Search( numAdmixNodes );
        optNet2.GetBestNet()->MapBackUserLabels( mapperPops );
        cout << "****** Highest log-probabiliyt of optimized network (searching over network space): " << logprobBest << endl;
        //cout << "The optimized network: ";
        //optNet2.GetBestNet()->Dump();
        
        // test code
        //AGHeuristicSearchEval agHeu(*optNet2.GetBestNet());
        //for(unsigned int i=0; i<listGeneTreePtrs.size(); ++i)
        //{
        //    double sc = agHeu.EvalTree(listGeneTreePtrs[i]);
        //    string strNW;
        //    listGeneTreePtrs[i]->ConsNewickSorted(strNW);
        //    cout << "Tree " << i << ":" << strNW << ", score: " << sc << endl;
        //}
//exit(1);
        
        
        // if network label has been changed prior to invoking this function, change it back now
        if( pmapperTaxaIdsUpperLevel != NULL )
        {
            //cout << "Before mapping back: trees=\n";
            //optNet.GetBestNet()->DumpMargTrees();
            
            optNet2.GetBestNet()->MapBackUserLabelsByUserLabels( *pmapperTaxaIdsUpperLevel );
        }
        
        
        long tmCur = GetElapseTime(tmLast);
        cout << "Time needed to find the optimal network: " << tmCur << endl;
        //cout << "Number of networks processed: " << optNet.GetNumProcNets() << endl;
        //cout << "Number of networks skipped: " << optNet.GetNumSkippedNets() << endl;
        //const char fileNetworkGML[100]="optimal-network.gml";
        optNet2.GetBestNet()->OutputGML(fileNetworkGMLOut);
        cout << "Optimal network: output to file: " << fileNetworkGMLOut << endl;
        
        // output the list of marginal trees
        cout << "List of marginal trees in the optimal network:" << endl;
        optNet2.GetBestNet()->DumpMargTrees();
        //cout << "Subtrees: \n";
        //map<set<int>, map<GenealogicalNetworkMTreeCode, GenealogicalNetworkMTreeInfo> > mapSubMargTreesWithFreqOpt;
        //optNet.GetBestNet()->RetriveAllSubMarginalTrees(szSubsetTaxa, mapSubMargTreesWithFreqOpt );
        
        // find all the admixture nodes
        optNet2.GetBestNet()->DumpAdmixNodes();
    }
    else
    {
        // YW: for now, assume taxa name are 1, 2, 3 and so on...
        // create a taxa mapper, map taxa x to x-1
        /*TaxaMapper mapperNetIn;
        for(int i=0; i<numTaxa; ++i)
        {
            mapperNetIn.AddTaxaString(std::to_string(i+1));
        }
cout << "mapperNetIn: ";
mapperNetIn.Dump(); */

//cout << "mapperPops: ";
//mapperPops.Dump();
        
        // read in the given network and optimize over it
        GenealogicalNetwork net1;
        //net1.ReadFromFileTreeMixFormat(netFileIn, &mapperTaxaIds);
        net1.ReadFromFileTreeMixFormat(netFileIn, &mapperPops);
//cout << "now creating optBr...\n";
        AGFastLenOpt optBr(net1, listGeneTreePtrs);
        const int MAX_OPT_ROUND = 3;
        optBr.SetMaxNumOptRounds(MAX_OPT_ROUND);
//cout << "now optimizing branch length...\n";
        //const int NUM_THREADS = 5;
        //optBr.SetMultithread(GetNumThreads());
        double logprobBest = optBr.Optimize();
cout << "logprobBest:" << logprobBest << endl;
        net1.MapBackUserLabels( mapperPops );
        cout << "****** Highest log-probabiliyt of optimized network (searching over network space): " << logprobBest << endl;
        //cout << "The optimized network: ";
        //net1.GetBestNet()->Dump();
        // if network label has been changed prior to invoking this function, change it back now
        if( pmapperTaxaIdsUpperLevel != NULL )
        {
            //cout << "Before mapping back: trees=\n";
            //optNet.GetBestNet()->DumpMargTrees();
            net1.MapBackUserLabelsByUserLabels( *pmapperTaxaIdsUpperLevel );
        }
        long tmCur = GetElapseTime(tmLast);
        cout << "Time needed to optimize branch length the network: " << tmCur << endl;
        //cout << "Number of networks processed: " << optNet.GetNumProcNets() << endl;
        //cout << "Number of networks skipped: " << optNet.GetNumSkippedNets() << endl;
        //const char fileNetworkGML[100]="optimal-network.gml";
        net1.OutputGML(fileNetworkGMLOut);
        cout << "Optimal network: output to file: " << fileNetworkGMLOut << endl;
        
        // output the list of marginal trees
        cout << "List of marginal trees in the optimal network:" << endl;
        net1.DumpMargTrees();
    }
    
    // dump out some stats
    AGGeneTreeProbDepot :: Instance().DumpStats();
    
    //SubtreeCoalProbDepot2::Instance().DumpStats();
    
    //SubtreeCoalACStore::Instance().DumpStats();
    
    //GTSubtreeDepot ::Instance().DumpStats();
    
    //STSubtreeDepot::Instance().DumpStats();
    
    AGProcessedNetsDepot::Instance().DumpStats();
    
    //ApproxGTPStats::Instance().DumpStats();
    
    ApproxGTPCache :: Instance().DumpStats();
    ApproxGTPCache2 :: Instance().DumpStats();
    
#if 0
    vector<GenealogicalNetworkNode *> listMixNodes;
    optNet.GetBestNet()->GetMixNodes( listMixNodes );
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
#endif
    
#if 0
    // // now output the extended newick format of the network
    bool fHVLines = true;
    string networkGMLNameTemp = fileNetworkGMLOut;
    networkGMLNameTemp += ".tmp";
    refineNetwork(fileNetworkGMLOut, networkGMLNameTemp, fHVLines);
    string networkGMLNameTemp2 = fileNetworkGMLOut;
    networkGMLNameTemp2 += ".tmp2";
    refineNetwork(networkGMLNameTemp,networkGMLNameTemp2,false);
    //const string networkGMLNameEWK ="HybridizationNetwork.ewk";
    string strEWK = convertGML2ExtNewick(networkGMLNameTemp2,false);
    //string strEWK = convertGML2ExtNewick(fileNetworkGMLOut,false);
    cout << "**************************************************************\n";
    cout << "The constructed network (in the extended Newick format): \n" << strEWK << endl;
#endif
    
    // cleanup
	for(int tr=0; tr<(int)listGeneTreePtrs.size(); ++tr)
	{
		delete listGeneTreePtrs[tr];
	}
	listGeneTreePtrs.clear();
}


static void Usage()
{
	cout << "Usage: ./gtmix2 <OPTIONS> -P <Population id file>  <GeneTreeFile> " << endl;
    cout << "List of options\n";
    cout << "\t -r <outgroup>:\t\t\t\t specify outgroup (default: no outgroup)\n";
    cout << "\t -n <number of admixture events>:\t specify the number of admixture events (default is 1 admixture)\n";
    cout << "\t -o <output admixture network name>:\t specify the name of the output admixture network file (in GML format). Default name is optimal-network.gml\n";
    cout << "\t -T <num of gene trees for inference>:\t set the number of gene trees used in inference (when the number of given input gene trees is large). Default is 500 trees.\n";
    cout << "\t -t <num of threads>:\t set the number of threads to use. Default is 1 (i.e., no parallization).\n";
	exit(1);
}

//const char *GTMIX_VER_INFO ="*** GTmix2 ver. 0.1.0.1, May 13, 2024 ***";
//const char *GTMIX_VER_INFO ="*** GTmix2 ver. 0.1.0.2, July 11, 2024 ***";
//const char *GTMIX_VER_INFO ="*** GTmix2 ver. 0.1.1.5, September 12, 2024 ***";
//const char *GTMIX_VER_INFO ="*** GTmix2 ver. 0.1.2.9, November 11, 2024 ***";
//const char *GTMIX_VER_INFO ="*** GTmix2 ver. 0.1.3.0, April 18, 2025 ***";
//const char *GTMIX_VER_INFO ="*** GTmix2 ver. 0.1.3.1, May 8, 2025 ***";
const char *GTMIX_VER_INFO ="*** GTmix2 ver. 1.0.0.0, June 12, 2025 ***";

//***********************************************************************

int main(int argc, char **argv)
{
    
    //TestGN3();
//return 1;
    //TestGN();
    //TestGN2();
    
    //YW_ASSERT_INFO( argc == 2, "Usage: must supply a file containing the list of gene trees" );
    
    cout << GTMIX_VER_INFO << endl << endl;
    
	// first verify usage
	if( CheckArguments(  argc, argv  ) == false)
	{
		Usage();
	}
    
    YW_ASSERT_INFO(filenameGenetrees != NULL, "Must provide gene trees.");
    
    long tstart1 = GetCurrentTimeTick();
    
    //if( fSummaryMode == true )
    //{
    //InfSummaryNet(filenameGenetrees);
    //}
    //else
    TestNetInference( filenameGenetrees );

    
    cout << "Elapsed time = " << GetElapseTime( tstart1 ) << " seconds." << endl;
    
    return 0;
}

