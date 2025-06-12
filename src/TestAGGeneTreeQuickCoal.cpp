// build test code:
// g++ TestAGGeneTreeQuickCoal.cpp GenealogicalNetwork.cpp Utils.cpp Utils2.cpp Utils3.cpp UnWeightedGraph.cpp PhylogenyTreeBasic.cpp MarginalTree.cpp Utils4.cpp AGGeneTreeQuickCoal.cpp  -std=c++11
#include "AGGeneTreeQuickCoal.hpp"
#include "GenealogicalNetwork.h"
#include "PhylogenyTreeBasic.h"
#include "AGFastLenOpt.hpp"
#include "AGFastCoalBuilder.hpp"
#include "Utils3.h"
using namespace std;

int GetNumThreads()
{
    int numThreads = 1;
    return numThreads;
}

void Test3()
{
    // create a network
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "((0,1),(2,3))";
    const int NUM_TAXA = 4;
    const int TAXON = 1;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == 2 )
        {
            pt2 = listLeaves[i];
        }
        
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    // find the parent of pt2
    vector<GenealogicalNetworkNode *> listAnces2;
    pt2->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt1->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // create a map of lineage count
    map<int,int> mapLinCnts;
    const int NUM_ALLELES = 2;
    for(int pop=0; pop<NUM_TAXA; ++pop)
    {
        mapLinCnts[pop] = NUM_ALLELES;
    }
    
    // creae a phylogenetic tree
    //string strTree = "((0,1),((0,3),(2,2)))";
    //string strTree = "((0,1),(1,(2,3)))";
    string strTree = "((0,0),(1,(2,3)))";
    PhylogenyTreeBasic phTree;
    phTree.ConsOnNewick(strTree);
    string strt2;
    phTree.ConsNewick(strt2);
    cout << "Constructed tree: " << strt2 << endl;
    
    // create a taxa mapper
    TaxaMapper taxMap;
    taxMap.AddTaxaStringWithId(0, "0");
    taxMap.AddTaxaStringWithId(1, "1");
    taxMap.AddTaxaStringWithId(2, "2");
    taxMap.AddTaxaStringWithId(3, "3");
    
    // run quick coal
    AGGeneTreeQuickCoal2 coalFinder( net );
    
    cout << "********************* Computing log-likelihood....\n";
    double loglikeli = coalFinder.CalcGeneTreeProbHeu(phTree);
    cout << "loglikeli: " << loglikeli << endl;
}

// test on input tree to see if works on larger trees
void Test4()
{
    cout << "******** TEST 4\n";
    // create a network with more taxa
    // create a network (top network on my board)
    const int NUM_TAXA = 6;
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "(((0,1),(2,3)),(4,5))";
    const int TAXON1 = 1;
    const int TAXON2 = 2;
    const int TAXON3 = 3;
    const int TAXON4 = 4;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    GenealogicalNetworkNode *pt3 = NULL;
    GenealogicalNetworkNode *pt4 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON1 )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON2 )
        {
            pt2 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON3 )
        {
            pt3 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON4 )
        {
            pt4 = listLeaves[i];
        }
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    YW_ASSERT_INFO(pt3!=NULL, "Taxon 3: didn't find");
    YW_ASSERT_INFO(pt4!=NULL, "Taxon 4: didn't find");
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt2->GetAnces1(), pt1->GetAnces1() );
    // find the parent of pt3
    vector<GenealogicalNetworkNode *> listAnces2;
    pt3->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches(  pt4->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // creae a phylogenetic tree with 30 leaves
    //string strTree = "((((((0,0),0),(1,1)),(5,0)),((((0,0),4),1),(((2,2),(2,3)),((2,3),(2,(3,3)))))),((((4,4),4),(4,((2,5),5))),((5,5),5)))";
    string strTree = "((((((0,1),0),(0,1)),(5,0)),((((0,4),0),1),(((2,2),(2,3)),((2,3),(2,(3,4)))))),((((3,4),4),(5,((2,5),5))),((4,5),5)))";
    PhylogenyTreeBasic phTree;
    phTree.ConsOnNewick(strTree);
    string strt2;
    phTree.ConsNewick(strt2);
    cout << "Constructed tree: " << strt2 << endl;
    
    // create a taxa mapper
    TaxaMapper taxMap;
    taxMap.AddTaxaStringWithId(0, "0");
    taxMap.AddTaxaStringWithId(1, "1");
    taxMap.AddTaxaStringWithId(2, "2");
    taxMap.AddTaxaStringWithId(3, "3");
    taxMap.AddTaxaStringWithId(4, "4");
    taxMap.AddTaxaStringWithId(5, "5");
    
    // create a map of lineage count
    map<int,int> mapLinCnts;
    mapLinCnts[0] = 6;
    mapLinCnts[1] = 3;
    mapLinCnts[2] = 6;
    mapLinCnts[3] = 4;
    mapLinCnts[4] = 5;
    mapLinCnts[5] = 6;
    
    // run quick coal
    AGGeneTreeQuickCoal2 coalFinder( net );
    
    cout << "********************* Computing log-likelihood....\n";
    double loglikeli = coalFinder.CalcGeneTreeProbHeu(phTree);
    cout << "loglikeli: " << loglikeli << endl;
}

// test on input tree to see if works on larger trees
void Test4a()
{
    cout << "******** TEST 4a\n";
    // create a network with more taxa
    // create a network (top network on my board)
    const int NUM_TAXA = 6;
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "(((0,1),(2,3)),(4,5))";
    const int TAXON1 = 1;
    const int TAXON2 = 2;
    const int TAXON3 = 3;
    const int TAXON4 = 4;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    GenealogicalNetworkNode *pt3 = NULL;
    GenealogicalNetworkNode *pt4 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON1 )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON2 )
        {
            pt2 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON3 )
        {
            pt3 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON4 )
        {
            pt4 = listLeaves[i];
        }
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    YW_ASSERT_INFO(pt3!=NULL, "Taxon 3: didn't find");
    YW_ASSERT_INFO(pt4!=NULL, "Taxon 4: didn't find");
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt2->GetAnces1(), pt1->GetAnces1() );
    // find the parent of pt3
    vector<GenealogicalNetworkNode *> listAnces2;
    pt3->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches(  pt4->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // creae a phylogenetic tree with 30 leaves
    //string strTree = "((((((0,0),0),(1,1)),(5,0)),((((0,0),4),1),(((2,2),(2,3)),((2,3),(2,(3,3)))))),((((4,4),4),(4,((2,5),5))),((5,5),5)))";
    string strTree = "(((((((0,0),1),0),(0,(1,1))),((4,5),0)),((((0,(5,4)),0),1),((((1,2),(2,3)),((3,2),3)),(((2,2),3),(2,(3,4)))))),(((((5,3),4),4),(5,((2,(3,5)),5))),((4,5),5)))";
    PhylogenyTreeBasic phTree;
    phTree.ConsOnNewick(strTree);
    string strt2;
    phTree.ConsNewick(strt2);
    cout << "Constructed tree: " << strt2 << endl;
    
    // create a taxa mapper
    TaxaMapper taxMap;
    taxMap.AddTaxaStringWithId(0, "0");
    taxMap.AddTaxaStringWithId(1, "1");
    taxMap.AddTaxaStringWithId(2, "2");
    taxMap.AddTaxaStringWithId(3, "3");
    taxMap.AddTaxaStringWithId(4, "4");
    taxMap.AddTaxaStringWithId(5, "5");
    
    // create a map of lineage count
    map<int,int> mapLinCnts;
    mapLinCnts[0] = 6;
    mapLinCnts[1] = 3;
    mapLinCnts[2] = 6;
    mapLinCnts[3] = 4;
    mapLinCnts[4] = 5;
    mapLinCnts[5] = 6;
    
    // run quick coal
    AGGeneTreeQuickCoal2 coalFinder( net );
    
    cout << "********************* Computing log-likelihood....\n";
    long tmCur = GetCurrentTimeTick();
    double loglikeli = 0.0;
    for(int k=0; k<1000; ++k)
    {
        loglikeli = coalFinder.CalcGeneTreeProbHeu(phTree);
    }
    cout << "Elapsed time: " << GetElapseTime(tmCur) << endl;
    cout << "loglikeli: " << loglikeli << endl;
}


void Test5()
{
    // create a network
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "((0,1),(2,3))";
    const int NUM_TAXA = 4;
    const int TAXON = 1;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == 2 )
        {
            pt2 = listLeaves[i];
        }
        
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    // find the parent of pt2
    vector<GenealogicalNetworkNode *> listAnces2;
    pt2->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt1->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // create a map of lineage count
    map<int,int> mapLinCnts;
    const int NUM_ALLELES = 2;
    for(int pop=0; pop<NUM_TAXA; ++pop)
    {
        mapLinCnts[pop] = NUM_ALLELES;
    }
    
    // creae a phylogenetic tree
    //string strTree = "((0,1),((0,3),(2,2)))";
    //string strTree = "((0,1),(1,(2,3)))";
    string strTree1 = "((0,0),(1,(2,3)))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);

    string strTree2 = "((0,1),(1,(2,3)))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    
    string strTree3 = "((0,2),(2,(1,3)))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
    
    vector<PhylogenyTreeBasic *> listTrees;
    listTrees.push_back( &phTree1 );
    listTrees.push_back( &phTree2 );
    listTrees.push_back( &phTree3 );
    
    // optimize branch length
    //AGGeneTreeQuickCoal2 coalFinder( net );
    AGFastLenOpt optBr(net, listTrees);
    double prob = optBr.Optimize();
    
    cout << "********************* Optimized likelihood after branch length optimization: " << prob << endl;
}

// test on input tree to see if works on larger trees
void Test6()
{
    cout << "******** TEST 6\n";
    // create a network with more taxa
    // create a network (top network on my board)
    const int NUM_TAXA = 6;
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "(((0,1),(2,3)),(4,5))";
    const int TAXON1 = 1;
    const int TAXON2 = 2;
    const int TAXON3 = 3;
    const int TAXON4 = 4;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    GenealogicalNetworkNode *pt3 = NULL;
    GenealogicalNetworkNode *pt4 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON1 )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON2 )
        {
            pt2 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON3 )
        {
            pt3 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON4 )
        {
            pt4 = listLeaves[i];
        }
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    YW_ASSERT_INFO(pt3!=NULL, "Taxon 3: didn't find");
    YW_ASSERT_INFO(pt4!=NULL, "Taxon 4: didn't find");
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt2->GetAnces1(), pt1->GetAnces1() );
    // find the parent of pt3
    vector<GenealogicalNetworkNode *> listAnces2;
    pt3->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches(  pt4->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // creae a phylogenetic tree with 30 leaves
    //string strTree = "((((((0,0),0),(1,1)),(5,0)),((((0,0),4),1),(((2,2),(2,3)),((2,3),(2,(3,3)))))),((((4,4),4),(4,((2,5),5))),((5,5),5)))";
    string strTree1 = "((((((0,1),0),(0,1)),(5,0)),((((0,4),0),1),(((2,2),(2,3)),((2,3),(2,(3,4)))))),((((3,4),4),(5,((2,5),5))),((4,5),5)))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);
    string strTree2 = "((((((0,2),1),(1,1)),(2,0)),((((0,4),0),1),(((2,3),(2,3)),((2,3),(3,(3,4)))))),((((3,4),4),(5,((2,5),5))),((4,4),5)))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    string strTree3 = "((((((0,0),0),(0,1)),(1,0)),((((0,3),0),1),(((2,2),(2,3)),((3,3),(2,(2,4)))))),((((4,4),4),(5,((3,5),5))),((5,5),5)))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
    
    // add two more trees
    string strTree4 = "((((((0,0),0),(0,1)),(1,0)),((((0,3),0),1),(((2,2),(2,3)),((3,3),(2,(2,4)))))),((((4,4),4),(5,((4,5),5))),((3,5),5)))";
    PhylogenyTreeBasic phTree4;
    phTree4.ConsOnNewick(strTree4);
    string strTree5 = "((((((0,1),0),(0,1)),(1,1)),((((1,3),0),1),(((2,3),(2,3)),((3,3),(2,(2,4)))))),((((4,4),4),(5,((4,5),5))),((3,5),5)))";
    PhylogenyTreeBasic phTree5;
    phTree5.ConsOnNewick(strTree5);
    
    vector<PhylogenyTreeBasic *> listTrees;
    listTrees.push_back( &phTree1 );
    listTrees.push_back( &phTree2 );
    listTrees.push_back( &phTree3 );
    listTrees.push_back( &phTree4 );
    listTrees.push_back( &phTree5 );
    
    // optimize branch length
    //AGGeneTreeQuickCoal2 coalFinder( net );
    AGFastLenOpt optBr(net, listTrees);
    const int NUM_THREADS = 5;
    optBr.SetMultithread(NUM_THREADS);
    double prob = optBr.Optimize();
    
    cout << "********************* Optimized likelihood after branch length optimization: " << prob << endl;
}

void Test7()
{
    // create a network
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "((0,1),(2,3))";
    const int NUM_TAXA = 4;
    const int TAXON = 1;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == 2 )
        {
            pt2 = listLeaves[i];
        }
        
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    // find the parent of pt2
    vector<GenealogicalNetworkNode *> listAnces2;
    pt2->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt1->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // create a map of lineage count
    map<int,int> mapLinCnts;
    const int NUM_ALLELES = 2;
    for(int pop=0; pop<NUM_TAXA; ++pop)
    {
        mapLinCnts[pop] = NUM_ALLELES;
    }
    
    // creae a phylogenetic tree
    //string strTree = "((0,1),((0,3),(2,2)))";
    //string strTree = "((0,1),(1,(2,3)))";
    string strTree1 = "((0,0),(1,(2,3)))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);

    string strTree2 = "((0,1),(1,(2,3)))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    
    string strTree3 = "((0,2),(2,(1,3)))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
    
    vector<PhylogenyTreeBasic *> listTrees;
    listTrees.push_back( &phTree1 );
    listTrees.push_back( &phTree2 );
    listTrees.push_back( &phTree3 );
    
    // create a taxa mapper
    TaxaMapper taxMap;
    taxMap.AddTaxaStringWithId(0, "0");
    taxMap.AddTaxaStringWithId(1, "1");
    taxMap.AddTaxaStringWithId(2, "2");
    taxMap.AddTaxaStringWithId(3, "3");
    //taxMap.AddTaxaStringWithId(4, "4");
    //taxMap.AddTaxaStringWithId(5, "5");
    
    long tmCur = GetCurrentTimeTick();
    
    // optimize branch length
    //AGGeneTreeQuickCoal2 coalFinder( net );
    //AGFastLenOpt optBr(net, listTrees);
    AGFastCoalBuilder opt(listTrees, net, taxMap );
    double prob = opt.Search();
    
    cout << "********************* Optimized likelihood after network search: " << prob << endl;
    cout << "The optimal network: \n";
    opt.GetBestNet()->DumpMargTrees(false);
    
    cout << "Elapsed time: " << GetElapseTime(tmCur) << endl;
}

// test on input tree to see if works on larger trees
void Test8()
{
    cout << "******** TEST 8\n";
    // create a network with more taxa
    // create a network (top network on my board)
    const int NUM_TAXA = 6;
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "(((0,1),(2,3)),(4,5))";
    const int TAXON1 = 1;
    const int TAXON2 = 2;
    const int TAXON3 = 3;
    const int TAXON4 = 4;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    GenealogicalNetworkNode *pt3 = NULL;
    GenealogicalNetworkNode *pt4 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON1 )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON2 )
        {
            pt2 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON3 )
        {
            pt3 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == TAXON4 )
        {
            pt4 = listLeaves[i];
        }
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    YW_ASSERT_INFO(pt3!=NULL, "Taxon 3: didn't find");
    YW_ASSERT_INFO(pt4!=NULL, "Taxon 4: didn't find");
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt2->GetAnces1(), pt1->GetAnces1() );
    // find the parent of pt3
    vector<GenealogicalNetworkNode *> listAnces2;
    pt3->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches(  pt4->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // creae a phylogenetic tree with 30 leaves
    //string strTree = "((((((0,0),0),(1,1)),(5,0)),((((0,0),4),1),(((2,2),(2,3)),((2,3),(2,(3,3)))))),((((4,4),4),(4,((2,5),5))),((5,5),5)))";
    string strTree1 = "((((((0,1),0),(0,1)),(5,0)),((((0,4),0),1),(((2,2),(2,3)),((2,3),(2,(3,4)))))),((((3,4),4),(5,((2,5),5))),((4,5),5)))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);
    string strTree2 = "((((((0,2),1),(1,1)),(2,0)),((((0,4),0),1),(((2,3),(2,3)),((2,3),(3,(3,4)))))),((((3,4),4),(5,((2,5),5))),((4,4),5)))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    string strTree3 = "((((((0,0),0),(0,1)),(1,0)),((((0,3),0),1),(((2,2),(2,3)),((3,3),(2,(2,4)))))),((((4,4),4),(5,((3,5),5))),((5,5),5)))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
    
    vector<PhylogenyTreeBasic *> listTrees;
    listTrees.push_back( &phTree1 );
    listTrees.push_back( &phTree2 );
    listTrees.push_back( &phTree3 );
    
    // create a taxa mapper
    TaxaMapper taxMap;
    taxMap.AddTaxaStringWithId(0, "0");
    taxMap.AddTaxaStringWithId(1, "1");
    taxMap.AddTaxaStringWithId(2, "2");
    taxMap.AddTaxaStringWithId(3, "3");
    taxMap.AddTaxaStringWithId(4, "4");
    taxMap.AddTaxaStringWithId(5, "5");
    
    long tmCur = GetCurrentTimeTick();
    
    // optimize branch length
    //AGGeneTreeQuickCoal2 coalFinder( net );
    //AGFastLenOpt optBr(net, listTrees);
    AGFastCoalBuilder opt(listTrees, net, taxMap );
    double prob = opt.Search();
    
    cout << "********************* Optimized likelihood after network search: " << prob << endl;
    cout << "The optimal network: \n";
    opt.GetBestNet()->DumpMargTrees(false);
    
    cout << "Elapsed time: " << GetElapseTime(tmCur) << endl;
}

void Test9()
{
    // create a network
    GenealogicalNetwork net;
    // init with a tree
    string strNW = "((0,1),(2,3))";
    const int NUM_TAXA = 4;
    const int TAXON = 1;
    net.InitWithTree(strNW, NUM_TAXA );
    //#if 0
    // get leaf with taxon 1
    vector<GenealogicalNetworkNode *> listLeaves;
    net.GetLeaves(listLeaves);
    GenealogicalNetworkNode *pt1 = NULL;
    GenealogicalNetworkNode *pt2 = NULL;
    for(unsigned int i=0; i<listLeaves.size(); ++i)
    {
        if( listLeaves[i]->GetTaxonId() == TAXON )
        {
            pt1 = listLeaves[i];
        }
        if( listLeaves[i]->GetTaxonId() == 2 )
        {
            pt2 = listLeaves[i];
        }
        
    }
    YW_ASSERT_INFO(pt1!=NULL, "Taxon 1: didn't find");
    YW_ASSERT_INFO(pt2!=NULL, "Taxon 2: didn't find");
    // find the parent of pt2
    vector<GenealogicalNetworkNode *> listAnces2;
    pt2->GetParents(listAnces2);
    YW_ASSERT_INFO( listAnces2.size() == 1, "WRONG" );
    GenealogicalNetworkNode *pt2p = listAnces2[0];
    GenealogicalNetworkBranch *pb2 = pt2p->GetAnces1();
    YW_ASSERT_INFO( pb2 != NULL, "WRONG2" );
    // add admixture from pb2 to pt1p
    net.AddMixBtwBranches( pt1->GetAnces1(), pb2 );
    //#endif
    
    //
    cout << "Network: ";
    //net.Dump();
    net.DumpMargTrees(false);
    
    // create a map of lineage count
    map<int,int> mapLinCnts;
    const int NUM_ALLELES = 2;
    for(int pop=0; pop<NUM_TAXA; ++pop)
    {
        mapLinCnts[pop] = NUM_ALLELES;
    }
    
    // creae a phylogenetic tree
    //string strTree = "((0,1),((0,3),(2,2)))";
    //string strTree = "((0,1),(1,(2,3)))";
#if 0
    string strTree1 = "((0,0),((1,2),3))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);
    string strTree2 = "((0,3),((0,2),1))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    string strTree3 = "(((0,1),0),(2,3))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
#endif
#if 0
    string strTree1 = "((0,0),(((1,1),(2,2)),(3,3)))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);
    string strTree2 = "((0,0),(((((1,1),2),2),3),3))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    //string strTree3 = "((((0,1),0),1),((2,2),(3,3)))";
    string strTree3 = "((((0,2),0),1),((1,2),(3,3)))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
#endif
    
#if 0
    string strTree1 = "(((0,0),(0,0)),((((1,1),(1,1)),((2,2),(2,2))),((3,3),(3,3))))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);
    string strTree2 = "((((0,0),0),0),(((((((1,1),(1,1)),2),2),(2,2)),(3,3)),(3,3)))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    //string strTree3 = "((((0,1),0),1),((2,2),(3,3)))";
    string strTree3 = "(((((0,0),(2,2)),(0,0)),(1,1)),(((1,1),(2,2)),((3,3),(3,3))))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
#endif
    string strTree1 = "(((0,0),(0,0)),((((1,1),(1,1)),((2,2),(2,2))),((3,3),(3,3))))";
    PhylogenyTreeBasic phTree1;
    phTree1.ConsOnNewick(strTree1);
    string strTree2 = "((((0,0),0),0),((((((1,1),(1,1)),2),(2,(2,2))),(3,3)),(3,3)))";
    PhylogenyTreeBasic phTree2;
    phTree2.ConsOnNewick(strTree2);
    //string strTree3 = "((((0,1),0),1),((2,2),(3,3)))";
    string strTree3 = "(((((0,0),2),(0,0)),(1,1)),(((1,(2,1)),(2,2)),((3,3),(3,3))))";
    PhylogenyTreeBasic phTree3;
    phTree3.ConsOnNewick(strTree3);
    
    // run quick coal
    AGGeneTreeQuickCoal2 coalFinder( net );
    
    cout << "********************* Computing log-likelihood....\n";
    double loglikeli = coalFinder.CalcGeneTreeProbHeu(phTree1);
    cout << "loglikeli for tree1: " << loglikeli << endl;
    loglikeli = coalFinder.CalcGeneTreeProbHeu(phTree2);
    cout << "loglikeli for tree2: " << loglikeli << endl;
    loglikeli = coalFinder.CalcGeneTreeProbHeu(phTree3);
    cout << "loglikeli for tree3: " << loglikeli << endl;
}

int main()
{
    //Test3();
    //Test4();
    //Test4a();
    //Test5();
    //Test6();
    //Test7();
    Test8();
    //Test9();
}
