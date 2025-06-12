//
//  GeneTreeProcUtils.cpp
//  
//
//  Created by Yufeng Wu on 4/24/19.
//
//

#include "GeneTreeProcUtils.hpp"
#include "PhylogenyTreeBasic.h"
#include <string>
#include <sstream>
#include <iterator>

static int numTreesToPick = 500;

void SetNumTreesToPick(int n) { numTreesToPick = n; }

//*****************************************************************************

void ConsPopHapToPopId( const vector<pair<string, vector<int> > > &listPopInfo, map<int,int> &mapPopRowToPopId  )
{
    //
    for(int i=0;i<(int)listPopInfo.size(); ++i)
    {
        for(int j=0;j<(int)listPopInfo[i].second.size(); ++j)
        {
            mapPopRowToPopId[ listPopInfo[i].second[j] ] = i;
            
//cout << "ConsPopHapToPopId: pop hap: " << listPopInfo[i].second[j] << ", to pop " << i << endl;
        }
    }
}


bool ReadInGeneTrees( ifstream &inFile, vector<PhylogenyTreeBasic *> &treePtrList, TaxaMapper *pTMapper )
{
    // read in RentPlus trees. Format: ignore lines starting with "//"; each tree can have an optional position then start with Newick format; empty lines: ignored
    int nLvs = -1;
    
    // read marginal trees in newick format
    // here there is no preamble, one line per tree
    while( inFile.eof() == false )
    {
        // ensure the first char is '('; otherwise stop
        // read in a file
        char buf[102400];
        inFile.getline(buf, sizeof(buf) );
        string strbuf=buf;
//cout << "Read one line: " << strbuf << endl;
        if(strbuf.length()==0)
        {
            continue;
        }
        if(strbuf.length()>=2 && strbuf[0]=='/' && strbuf[1]=='/')
        {
            continue;
        }
        string treeNewick;
        char ch=strbuf[0];
        if( ch != '(' )
        {
            // look for the second term
            istringstream iss(strbuf);
            std::vector<std::string> strTmps((std::istream_iterator<std::string>(iss)),
                                             std::istream_iterator<std::string>());
            YW_ASSERT_INFO(strTmps.size()>=2, "File format wrong: tree file missing");
            treeNewick = strTmps[1];
        }
        else
        {
            treeNewick = strbuf;
        }
        if( treeNewick.size() == 0 )
        {
            break;
        }
        if( treeNewick[0] != '(')
        {
            YW_ASSERT_INFO(false, "File format wrong: tree file format wrong");
        }
        
        //cout << "newick tree = " << treeNewick << endl;
        
        //#if 0
        // update numleaves
        multiset<string> setLabels;
        NewickUtils :: RetrieveLabelSet( treeNewick, setLabels);
#if 0
        for(multiset<string> :: iterator it22 = setLabels.begin(); it22 != setLabels.end(); ++it22)
        {
            cout << "Label found: " << *it22 << endl;
        }
#endif
        nLvs = setLabels.size();
        //#endif
        //
        PhylogenyTreeBasic *pphTree = new PhylogenyTreeBasic;
        //if( fDup == false )
        //{
        pphTree->ConsOnNewick(treeNewick, -1, false, pTMapper);
        //cout << "Done phylogenetic tree construction...\n";
        //pphTree->OutputGML("tmp.gml");
        //}
        //else
        //{
        //	phTree.ConsOnNewickDupLabels(treeNewick, pTMapper);
        //}
        
        if( pTMapper != NULL)
        {
            pTMapper->SetInitialized(true);
        }
        //string strTr;
        //pphTree->ConsNewick(strTr);
        //cout << "After reconstruction: strTr = " << strTr << endl;
        //see if zero is in, if not, must have 1 and decrease by 1
        set<int> lvids;
        pphTree->GetLeaveIds(lvids);
        //cout << "lvids : ";
        //DumpIntSet( lvids );
        int idInternal = lvids.size();
        YW_ASSERT_INFO( lvids.find(0) != lvids.end(), "Must adjust leaf label first (to start with 0)" );
        
        //	YW_ASSERT_INFO( lvids.find(1) != lvids.end(), "Wrong" );
        
        // decrease by one
        PhylogenyTreeIterator itorTree(*pphTree);
        itorTree.Init();
        //pphTree->InitPostorderWalk();
        while(itorTree.IsDone() == false)
        {
            //				TreeNode *pn = pphTree->NextPostorderWalk( ) ;
            TreeNode *pn = itorTree.GetCurrNode() ;
            itorTree.Next();
            if( pn == NULL)
            {
                break;      // done with all nodes
            }
            if( pn->IsLeaf() == false )
            {
                pn->SetID( idInternal++ );
            }
        }
        
        // mark the change
        //	fNoChange = false;
        
        vector<int> nidsList, nparsList;
        pphTree->GetNodeParInfo(nidsList, nparsList);
        //phTree.GetNodeParInfoNew(nidsList, nparsList);
        //phTree.GetNodeParInfo(nidsList, nparsList);
        //if( nLvs <= 0 )
        //{
        //string strTrNW;
        //pphTree->ConsNewick(strTrNW);
        //cout << "strTrNW: " << strTrNW << endl;
        treePtrList.push_back( pphTree );
        //cout << "Newick format of this marginal tree: ";
        //cout << tree.GetNewick() << endl;
    }
    return true;
}

void AdjustLabelsForPopInfo( vector<PhylogenyTreeBasic *> &listGeneTreePtrs, const vector<pair<string, vector<int> > > &listPopInfo, TaxaMapper *pTMapperTree )
{
    // create map
    map<int,int> mapPopRowToPopId ;
    ConsPopHapToPopId( listPopInfo, mapPopRowToPopId );
    set<int> setLvLabels;
    YW_ASSERT_INFO(listGeneTreePtrs.size() > 0, "Must have at least one tree");
    listGeneTreePtrs[0]->GetLeafIntLabels(setLvLabels);
    map<int,int> mapOldLblsToNewLbls;
    for(set<int> :: iterator it = setLvLabels.begin(); it != setLvLabels.end(); ++it)
    {
        int lblOld = *it;
        string lblStrOrig=pTMapperTree->GetString(lblOld);
        int lblUser = std::stoi( lblStrOrig );
//cout << "lblOld: " << lblOld << ", lblUser: " << lblUser << endl;
        YW_ASSERT_INFO( mapPopRowToPopId.find(lblUser)!=mapPopRowToPopId.end(), "Fail to find" );
        int lblNew = mapPopRowToPopId[lblUser];
        mapOldLblsToNewLbls[lblOld] = lblNew;
    }
    
    
    // map the raw row index to population index (starting from 0, 1, ....
    for(int tr=0; tr<(int)listGeneTreePtrs.size(); ++tr)
    {
        //
        set<int> setLvLabels;
        listGeneTreePtrs[tr]->GetLeafIntLabels(setLvLabels);
//cout << "Tree " << tr << ": before changing, ";
//string strTreePrev;
//listGeneTreePtrs[tr]->ConsNewickSorted(strTreePrev);
//cout << strTreePrev << endl;
        
        ChangeLeafIntLabelOfTree(*listGeneTreePtrs[tr], mapOldLblsToNewLbls, true);
        
//cout << "After changing: tree: ";
//string strTreeAfter;
//listGeneTreePtrs[tr]->ConsNewickSorted(strTreeAfter);
//cout << strTreeAfter << endl;
    }
    
    
    
//exit(1);
}

void AdjustLabelsTaxaMapper( vector<PhylogenyTreeBasic *> &listGeneTreePtrs, TaxaMapper *pTMapperTree )
{
    // create map
    set<int> setLvLabels;
    YW_ASSERT_INFO(listGeneTreePtrs.size() > 0, "Must have at least one tree");
    listGeneTreePtrs[0]->GetLeafIntLabels(setLvLabels);
    map<int,int> mapOldLblsToNewLbls;
    for(set<int> :: iterator it = setLvLabels.begin(); it != setLvLabels.end(); ++it)
    {
        int lblOld = *it;
        string lblStrOrig=pTMapperTree->GetString(lblOld);
        int lblUser = std::stoi( lblStrOrig );
        //cout << "lblOld: " << lblOld << ", lblUser: " << lblUser << endl;
        mapOldLblsToNewLbls[lblOld] = lblUser;
    }
    
    
    // map the raw row index to population index (starting from 0, 1, ....
    for(int tr=0; tr<(int)listGeneTreePtrs.size(); ++tr)
    {
        //
        set<int> setLvLabels;
        listGeneTreePtrs[tr]->GetLeafIntLabels(setLvLabels);
        //cout << "Tree " << tr << ": before changing, ";
        //string strTreePrev;
        //listGeneTreePtrs[tr]->ConsNewickSorted(strTreePrev);
        //cout << strTreePrev << endl;
        
        ChangeLeafIntLabelOfTree(*listGeneTreePtrs[tr], mapOldLblsToNewLbls, true);
        
        //cout << "After changing: tree: ";
        //string strTreeAfter;
        //listGeneTreePtrs[tr]->ConsNewickSorted(strTreeAfter);
        //cout << strTreeAfter << endl;
    }
    
    
    
    //exit(1);
}



//***********************************************************************
// Yet another way to pick tree: score trees based on their cluster frequency

static double ScoreTreeFromClusFreq(PhylogenyTreeBasic *pTree, map<set<int>, int> &mapClusFreqAll, int totWt )
{
    double res=1.0;
    map<TreeNode *, set<int> > mapClusOfNode;
    PhylogenyTreeIterator itor(*pTree);
    itor.Init();
    while(itor.IsDone()==false)
    {
        //
        TreeNode *pn=itor.GetCurrNode();
        if(pn->IsLeaf())
        {
            set<int> ss;
            ss.insert(pn->GetIntLabel());
            mapClusOfNode[pn]=ss;
        }
        else
        {
            // union of the two children and save
            int nc=pn->GetChildrenNum();
            set<int> ssCombo;
            for(int i=0;i<nc;++i)
            {
                map<TreeNode *, set<int> > :: iterator it1 = mapClusOfNode.find( pn->GetChild(i) );
                YW_ASSERT_INFO(it1 != mapClusOfNode.end(), "Fail to find");
                UnionSets(ssCombo, it1->second);
            }
            mapClusOfNode[pn]=ssCombo;
            
            YW_ASSERT_INFO( mapClusFreqAll.find(ssCombo)!=mapClusFreqAll.end(), "Fail to find" );
            int numOccClus=mapClusFreqAll[ssCombo];
            res *= ((double)numOccClus)/totWt;
        }
        
        itor.Next();
    }
    
    return res;
}

static void AnalyzeClusterOfTree( PhylogenyTreeBasic *pTree, map<set<int>, int> &mapClusFreq )
{
    // deal only with internal node
    map<TreeNode *, set<int> > mapClusOfNode;
    PhylogenyTreeIterator itor(*pTree);
    itor.Init();
    while(itor.IsDone()==false)
    {
        //
        TreeNode *pn=itor.GetCurrNode();
        if(pn->IsLeaf())
        {
            set<int> ss;
            ss.insert(pn->GetIntLabel());
            mapClusOfNode[pn]=ss;
        }
        else
        {
            // union of the two children and save
            int nc=pn->GetChildrenNum();
            set<int> ssCombo;
            for(int i=0;i<nc;++i)
            {
                map<TreeNode *, set<int> > :: iterator it1 = mapClusOfNode.find( pn->GetChild(i) );
                YW_ASSERT_INFO(it1 != mapClusOfNode.end(), "Fail to find");
                UnionSets(ssCombo, it1->second);
            }
            mapClusOfNode[pn]=ssCombo;
            
            if(mapClusFreq.find(ssCombo)==mapClusFreq.end())
            {
                mapClusFreq[ssCombo]=0;
            }
            ++mapClusFreq[ssCombo];
        }
        
        itor.Next();
    }
    //string strNW;
    //pTree->ConsNewickSorted(strNW);
    //cout << "Tree: " << strNW << ", clusters frequencey: \n";
    //for(map<set<int>,int> :: iterator it = mapClusFreq.begin(); it != mapClusFreq.end(); ++it )
    //{
    //cout << "[" << it->second << "], clsuter ";
    //DumpIntSet(it->first);
    //}
}

static void PickTreesImp3( const vector<PhylogenyTreeBasic *> &listGeneTreePtrs, int numTreesToChoseHere, vector<PhylogenyTreeBasic *> &listChosenTrees)
{
    //
    //int numTotTree = listGeneTreePtrs.size();
    map<string, int> mapTopoFreq;
    map<PhylogenyTreeBasic *, string> mapTreeTopo;
    map<string, set<PhylogenyTreeBasic *> > mapTreeTopoTrees;
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        string strNW;
        listGeneTreePtrs[i]->ConsNewickSorted( strNW );
        mapTreeTopo[ listGeneTreePtrs[i] ] = strNW;
        if(mapTopoFreq.find(strNW) == mapTopoFreq.end() )
        {
            mapTopoFreq[strNW] = 0;
        }
        ++mapTopoFreq[strNW];
        mapTreeTopoTrees[strNW].insert( listGeneTreePtrs[i] );
    }
    //cout << "List of freqencies: \n";
    //for(map<string,int> :: iterator it=mapTopoFreq.begin(); it!=mapTopoFreq.end(); ++it)
    //{
    //cout << "[" << it->second << "]: " << it->first << endl;
    //}
    
    // accumulate clusters
    map<set<int>,int> mapClusFreqAll;
    int totWt=0;
    for( map<string, set<PhylogenyTreeBasic *> > :: iterator it = mapTreeTopoTrees.begin(); it != mapTreeTopoTrees.end(); ++it  )
    {
        YW_ASSERT_INFO(it->second.size()>0, "Empty");
        PhylogenyTreeBasic *ptree = *( it->second.begin() );
        map<set<int>,int> mapClusFreqStep;
        AnalyzeClusterOfTree(ptree, mapClusFreqStep);
        
        int numTreeThis = mapTopoFreq[it->first];
        
        for( map<set<int>,int>::iterator it2= mapClusFreqStep.begin();it2!=mapClusFreqStep.end(); ++it2 )
        {
            if(mapClusFreqAll.find(it2->first)==mapClusFreqAll.end())
            {
                mapClusFreqAll[it2->first]=0;
            }
            double wtStep=it2->second * numTreeThis;
            mapClusFreqAll[it2->first] += wtStep;
            totWt+=wtStep;
        }
    }
    //cout << "Overall, clusters frequencey: \n";
    //for(map<set<int>,int> :: iterator it = mapClusFreqAll.begin(); it != mapClusFreqAll.end(); ++it )
    //{
    //cout << "[" << it->second << "], clsuter ";
    //DumpIntSet(it->first);
    //}
    
    // now score the trees based on cluster frq
    vector< pair<double, void * > > listTreeWithScore;
    for( map<string, set<PhylogenyTreeBasic *> > :: iterator it = mapTreeTopoTrees.begin(); it != mapTreeTopoTrees.end(); ++it  )
    {
        PhylogenyTreeBasic *ptree = *( it->second.begin() );
        double score = ScoreTreeFromClusFreq(ptree, mapClusFreqAll, totWt);
        for(set<PhylogenyTreeBasic*> :: iterator it2 = it->second.begin(); it2 !=it->second.end(); ++it2)
        {
            pair<double, void *> pp( score, (void *)(*it2) );
            listTreeWithScore.push_back(pp);
        }
    }
    
    // sort
    SortPairsByNumsDouble(listTreeWithScore);
    
    // pick from the back
    int numOutTree=0;
    listChosenTrees.clear();
    bool fPick=true;
    for(int i=(int)listTreeWithScore.size()-1; i>=0; --i)
    {
        ++numOutTree;
        //string strNW=mapTreeTopo[ (PhylogenyTreeBasic *)listTreeWithScore[i].second ];
        //cout << strNW << endl;
        if( fPick )
        {
            listChosenTrees.push_back( (PhylogenyTreeBasic *)listTreeWithScore[i].second );
        }
        else
        {
            delete (PhylogenyTreeBasic *)listTreeWithScore[i].second;
        }
        if(numOutTree >=numTreesToChoseHere )
        {
            fPick=false;
        }
    }
}

//***********************************************************************
// Pick trees by narrowing down: idea is trying to remove noisy trees step by step

static void PickTreesImpStepwise( const vector<PhylogenyTreeBasic *> &listGeneTreePtrs, vector<PhylogenyTreeBasic *> &listChosenTrees)
{
    // drop half of trees each time until reaching the number of needed trees
    vector<PhylogenyTreeBasic *> listTreesToChooseCurr = listGeneTreePtrs;
    while( true )
    {
        //
        int numTreesAvail = (int)listTreesToChooseCurr.size();
        int numTreesToChoseHere = numTreesAvail/2;
        if( numTreesAvail <= 2*numTreesToPick )
        {
            numTreesToChoseHere = numTreesToPick;
        }
        vector<PhylogenyTreeBasic *> listTreesToChooseCurrNext;
        PickTreesImp3(listTreesToChooseCurr, numTreesToChoseHere, listTreesToChooseCurrNext);
        if( (int)listTreesToChooseCurrNext.size() <= numTreesToPick )
        {
            listChosenTrees = listTreesToChooseCurrNext;
            break;
        }
        else
        {
            listTreesToChooseCurr = listTreesToChooseCurrNext;
        }
    }
    
}



//***********************************************************************


void PickTrees(vector<PhylogenyTreeBasic *> &listGeneTreePtrsToProc)
{
    // construct initial network; for now, just
    YW_ASSERT_INFO( listGeneTreePtrsToProc.size() > 0, "Tree list is empty" );
    
    
    // perform selection
    //PickTreesImp( listGeneTreePtrs, mapNontrivSites );
    vector<PhylogenyTreeBasic *> listTreesPicked;
    //PickTreesImp3( listGeneTreePtrs, numTreesToPick, listTreesPicked );
    PickTreesImp3( listGeneTreePtrsToProc, numTreesToPick,  listTreesPicked );
    
    listGeneTreePtrsToProc = listTreesPicked;
}




