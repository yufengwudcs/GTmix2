#include "DeepCoalescence.h"
#include "Utils3.h"

extern void GetLeavesSpeciesFromGT( PhylogenyTreeBasic &treeGene, set<int> &species );

///////////////////////////////////////////////////////////////////////////////////
static int maxSrcEntryNumStored = 20;

static void SetMaxSrcEntryNumStored(int numSrcEntries)
{
	maxSrcEntryNumStored = numSrcEntries;
}

static int GetMaxSrcEntryNumStored()
{
	return maxSrcEntryNumStored;
}


// useful functions
void CollectSpecies(TreeNode *pn, multiset<int> &listSpecies)
{
	// 
	YW_ASSERT_INFO( pn != NULL, "Can not be empty" );
	listSpecies.clear();
	vector<string> listLeafLabels;
	pn->GetAllLeafLabeles(listLeafLabels);
	// convert to int
	for(int i=0; i<(int)listLeafLabels.size(); ++i)
	{
		int id;
		sscanf( listLeafLabels[i].c_str(), "%d", &id );
		listSpecies.insert(id);
	}
}

void GetAllClusterForTree(PhylogenyTreeBasic &tree, set< set<int> > &clusters)
{
	PhylogenyTreeIterator titor( tree );
	titor.Init();
	while( titor.IsDone() == false )
	{
		// get current node
		TreeNode *pnCurr = titor.GetCurrNode();
		multiset<int> setLabels;
		CollectSpecies(pnCurr, setLabels);
		// add record
		//mapNodesUnder.insert( map<TreeNode *, multiset<int> > :: value_type( pnCurr, setLabels ) ) ;


		//YW_ASSERT_INFO( mapNodesUnder.find( pnCurr ) != mapNodesUnder.end(), "Fail" );
		//multiset<int> setLabels = mapNodesUnder[pnCurr];
		// 
		set<int> setDistinctLabels;
		ConvMSetToSet(setLabels, setDistinctLabels);
		// 
		clusters.insert(setDistinctLabels);

		titor.Next();
	}
}

///////////////////////////////////////////////////////////////////////////////////

MDCTableEntry :: MDCTableEntry() : fDirty(false)
{
	// do nothing
}

// access
//void MDCTableEntry :: CheckSol(int val)
//{
//	TBD????
//}

void MDCTableEntry :: Dump() 
{
	for( map<int, pair<int, set<MDCEntrySrc> > > :: iterator itt = mapNearpOptInfo.begin(); itt != mapNearpOptInfo.end(); ++itt )
	{
		cout << "Near optial value: " << itt->first << ", new optimal tree num: " << itt->second.first << ", src config: \n";
		for(set<MDCEntrySrc> :: iterator itt2  = itt->second.second.begin(); itt2 != itt->second.second.end();  ++itt2)
		{
			cout << "src level1 = " << itt2->first.second << ", level2 = " << itt2->second.second << ", subset1 = ";
			DumpIntSet( itt2->first.first );
		}
	}

#if 0
	cout << "Near optimal value: ";
	DumpIntVec( valsNearOpt );
	cout << "New optimal tree num: ";
	DumpIntVec( counstNearOpt );
	cout << "Src configs: \n";
	for( int i=0; i<(int)listSrcClusters.size(); ++i  )
	{
		for(set<MDCEntrySrc> :: iterator itt  = listSrcClusters[i].begin(); itt != listSrcClusters[i].end();  ++itt)
		{
			cout << "src level1 = " << itt->first.second << ", level2 = " << itt->second.second << ", subset1 = ";
			DumpIntSet( itt->first.first );
		}
	}
#endif
}


void MDCTableEntry :: AddItem(int levelLimit, const MDC_Cluster &clusterSrc1,int ss1, const MDC_Cluster &clusterSrc2, int ss2, int val, int cnt)
{ 
//cout << "AddItem: val = " << val << ", cnt = " << cnt << endl;
//Dump();

	// see if this value is in 
	MDCEntrySrc mdceNew;
	mdceNew.first.first = clusterSrc1;
	mdceNew.second.first = clusterSrc2;
	mdceNew.first.second = ss1;
	mdceNew.second.second = ss2;
	if( clusterSrc1 > clusterSrc2)
	{
		mdceNew.first.first = clusterSrc2;
		mdceNew.second.first = clusterSrc1;
		mdceNew.first.second = ss2;
		mdceNew.second.second = ss1;
	}

	// are we having a new entry?
	if( mapNearpOptInfo.find( val ) == mapNearpOptInfo.end() )
	{
		// if so, create a new entry to it
		set<MDCEntrySrc> mdcSingle;
		mdcSingle.insert( mdceNew );
		pair<int, set<MDCEntrySrc> > pp(cnt, mdcSingle);
		mapNearpOptInfo.insert( map<int, pair<int, set<MDCEntrySrc> > > :: value_type(val, pp) );
	}
	else
	{
		// append to the existing entry
		mapNearpOptInfo[val].first += cnt; 
		if( (int) mapNearpOptInfo[val].second.size() < GetMaxSrcEntryNumStored() )
		{
			mapNearpOptInfo[val].second.insert( mdceNew );
		}
	}

	// if we have too many, get rid the last one
	if( (int)mapNearpOptInfo.size() > levelLimit )
	{
		int keyLarger = mapNearpOptInfo.rbegin()->first;
		mapNearpOptInfo.erase( keyLarger );
	}
	YW_ASSERT_INFO( (int)mapNearpOptInfo.size() <= levelLimit, "Fatal error2" );



#if 0
	//YW_ASSERT_INFO(rank<(int)counstNearOpt.size(), "Overflow"); counstNearOpt[rank] += cnt;  
	// create a sorted list by the values
	vector<int> listValsNew;
	vector<int> listNumsNew;
	vector< set<MDCEntrySrc> > listSrcClustersNew;
	bool fAdded = false;
	MDCEntrySrc mdceNew;
	mdceNew.first.first = clusterSrc1;
	mdceNew.second.first = clusterSrc2;
	mdceNew.first.second = ss1;
	mdceNew.second.second = ss2;
	//pair<MDC_Cluster, MDC_Cluster> ppNew(clusterSrc1, clusterSrc2);
	if( clusterSrc1 > clusterSrc2)
	{
		mdceNew.first.first = clusterSrc2;
		mdceNew.second.first = clusterSrc1;
		mdceNew.first.second = ss2;
		mdceNew.second.second = ss1;
	}

	for( int i=0; i<(int)valsNearOpt.size(); ++i )
	{
		// 
		if( valsNearOpt[i] == val )
		{
			listValsNew.push_back( val);
			listNumsNew.push_back( counstNearOpt[i]+ cnt);
			set<MDCEntrySrc>  vecPPs = listSrcClusters[i];
			if( (int) vecPPs.size() < GetMaxSrcEntryNumStored() )
			{
				vecPPs.insert( mdceNew );
			}
			listSrcClustersNew.push_back( vecPPs );
			fAdded = true;
		}
		else 
		{
			if( valsNearOpt[i] > val && fAdded == false )
			{
				listValsNew.push_back( val);
				listNumsNew.push_back( cnt);
				set<MDCEntrySrc>  vecPPs;
				vecPPs.insert( mdceNew );
				listSrcClustersNew.push_back( vecPPs );
				fAdded = true;
			}

			// add the current val anyway
			listValsNew.push_back( valsNearOpt[i] );
			listNumsNew.push_back( counstNearOpt[i] );
			listSrcClustersNew.push_back( listSrcClusters[i] );

		}
	}
	// add it if not yet
	if( fAdded == false)
	{
		listValsNew.push_back( val);
		listNumsNew.push_back( cnt);
		set<MDCEntrySrc>  vecPPs;
		vecPPs.insert( mdceNew );
		listSrcClustersNew.push_back(vecPPs);
		fAdded = true;
	}

	// if over the limit, pop back
	if( (int)listValsNew.size() > levelLimit )
	{
		listValsNew.pop_back();
		listNumsNew.pop_back();
		listSrcClustersNew.pop_back();
	}
	YW_ASSERT_INFO( (int)listNumsNew.size() <= levelLimit && listNumsNew.size() == listValsNew.size(), "Fatal error2" );

	valsNearOpt = listValsNew;
	counstNearOpt = listNumsNew;
	listSrcClusters = listSrcClustersNew;
#endif
	// finally mark entry to be dirty
	fDirty = true;

//cout << "After AddItem\n";
//Dump();
}

int MDCTableEntry :: GetVal(int i) const 
{ 
	if( fDirty == true )
	{
		const_cast<MDCTableEntry *>(this)->Init();
	}
	YW_ASSERT_INFO( i<(int)vecOrderedNearOptVals.size(), "FATAL ERROR" );
	return vecOrderedNearOptVals[i]; 
}

int MDCTableEntry :: GetCount(int i) const 
{ 
	int val = GetVal(i);
	return const_cast<MDCTableEntry *>(this)->mapNearpOptInfo[val].first;

	//return counstNearOpt[i]; 
}

int MDCTableEntry :: GetSrcOptionsNum(int i)  const
{ 
	int val = GetVal(i);
	return const_cast<MDCTableEntry *>(this)->mapNearpOptInfo[val].second.size();
//	return listSrcClusters[i].size();  
}


void MDCTableEntry :: GetSrcOptionsAt(int i, set<MDCEntrySrc> &subsetOptions)
{  
	if( i >=  (int)vecOrderedNearOptVals.size() )
	{
		// no more to explore
		return;
	}

	int val = GetVal(i);
	subsetOptions = mapNearpOptInfo[val].second;
	//subsetOptions = listSrcClusters[i];
}

void MDCTableEntry:: Init()
{
	// when record is dirty, re-initialize it
	vecOrderedNearOptVals.clear();
	for( map<int, pair<int, set<MDCEntrySrc> > > :: iterator itt = mapNearpOptInfo.begin(); itt != mapNearpOptInfo.end(); ++itt )
	{
		vecOrderedNearOptVals.push_back( itt->first );
	}
	fDirty = false;
//cout << "vecOrderedNearOptVals: ";
//DumpIntVec( vecOrderedNearOptVals );
}


///////////////////////////////////////////////////////////////////////////////////
// this is the main DP table

MDCTable :: MDCTable(int nl ) : numLevels(nl)
{
}

void MDCTable :: Dump() 
{
	for( map< MDC_Cluster, MDCTableEntry > :: iterator it = mdcTable.begin(); it != mdcTable.end(); ++it )
	{	
		// 
		cout << "Cluster: ";
		DumpIntSet( it->first );
		it->second.Dump();
	}
}


void MDCTable :: AddEntry( const MDC_Cluster &cluster, const MDC_Cluster &clusterSrc1, int ss1, const MDC_Cluster &clusterSrc2, int ss2, int val, int num )
{
//cout << "AddEntry: val = " << val << ", num = " << num << " for cluster: ";
//DumpIntSet( cluster);
	// create an entry if not there yet
	if( mdcTable.find( cluster ) == mdcTable.end() )
	{
		MDCTableEntry mdcEntry;
		mdcTable.insert( map< MDC_Cluster, MDCTableEntry > :: value_type( cluster, mdcEntry ) );
	}
	// add it: be careful, need to ensure we donot have too many levels
	mdcTable[cluster].AddItem( numLevels, clusterSrc1, ss1, clusterSrc2, ss2, val, num );
}

int MDCTable :: GetValAt(const MDC_Cluster &cluster, int rank ) const
{
    if( mdcTable.find( cluster ) == mdcTable.end() )
    {
        //cout << "Fail to find in MDC table: ";
        //DumpIntSet( cluster );
        return -1;
    }
	//YW_ASSERT_INFO( mdcTable.find( cluster ) != mdcTable.end(), "GetValAt: Cluster: not found" );
	if(rank >=  const_cast<MDCTable *>(this)->mdcTable[cluster].GetNumItems() )
	{
		return -1;
	}
	return const_cast<MDCTable *>(this)->mdcTable[cluster].GetVal(rank);
}

int MDCTable :: GetNumAt(const MDC_Cluster &cluster, int rank )
{
	YW_ASSERT_INFO( mdcTable.find( cluster ) != mdcTable.end(), "GetNumAt: Cluster: not found" );
	if(rank >= mdcTable[cluster].GetNumItems() )
	{
		return 0;
	}
	return mdcTable[cluster].GetCount(rank);
}


void MDCTable :: RetrieveNearOptSTrees( const MDC_Cluster &cluster, int mdcLevel, set<string> &listNewickTrees, int maxNeeded )
{
//cout << "RetrieveNearOptSTrees: mdcLevel = " << mdcLevel << ", cluster = ";
//DumpIntSet( cluster );
	YW_ASSERT_INFO( mdcTable.find( cluster ) != mdcTable.end(), "Cluster: not found" );
	listNewickTrees.clear();
	// retrive the optimal solutions here. We try all subsets
	// use recurisve call
	YW_ASSERT_INFO(cluster.size() >= 1, "Can not have empty cluster");
	if( cluster.size() == 1)
	{
		// only return a single item
		char buf[100];
		snprintf(buf, 100, "%d", *(cluster.begin() ) );
		string strName = buf;
		listNewickTrees.insert( strName );
	}
	else
	{
		if( mdcLevel >= (int)mdcTable[cluster].GetNumItems() )
		{
			//cout << "Cluster: ";
			//DumpIntSet( cluster );
			//cout << "mdcLevel = " << mdcLevel << ", mdcTable[cluster].GetNumItems() = " << mdcTable[cluster].GetNumItems() << endl;
			// can not find it
			return;
		}

		//YW_ASSERT_INFO( mdcLevel <(int)mdcTable[cluster].GetNumItems(), "Error in level" );
		// now consider all entries in the item
		//for(int tt = 0; tt<mdcTable[cluster].GetNumItems(); ++tt )
		//{
		set<MDCEntrySrc>  subsetOptions;
		mdcTable[cluster].GetSrcOptionsAt(mdcLevel, subsetOptions);

		// treat each items
		for( set<MDCEntrySrc> :: iterator itt3 =  subsetOptions.begin(); itt3 != subsetOptions.end(); ++itt3 )
//			for(int kk=0; kk<mdcTable[cluster].GetSrcOptionsNum(tt); ++kk)
		{
			MDC_Cluster srcss1 = itt3->first.first, srcss2 = itt3->second.first;
			int level1 = itt3->first.second, level2 = itt3->second.second;
			//mdcTable[cluster].GetSrcOptionsAt(tt, kk, srcss1, srcss2);
			// look down
			set<string> listSS1, listSS2;
			RetrieveNearOptSTrees( srcss1, level1, listSS1, maxNeeded  );
			RetrieveNearOptSTrees( srcss2, level2, listSS2, maxNeeded  );
			// combine it
			for( set<string> :: iterator ittg1 = listSS1.begin(); ittg1 != listSS1.end(); ++ittg1 )
			{
				for( set<string> :: iterator ittg2 = listSS2.begin(); ittg2 != listSS2.end(); ++ittg2 )
				{
					// 
					if( (int)listNewickTrees.size() <  maxNeeded )
					{
						string strTree = "(";
						strTree += *ittg1;
						strTree += ",";
						strTree += *ittg2;
						strTree += ")";

						// save it
						//YW_ASSERT_INFO(listNewickTrees.find( strTree) == listNewickTrees.end(), "Can not have duplicate");
						listNewickTrees.insert( strTree );
					}
					else
					{
						break;
					}
				}
			}
		}
		//}
	}
//cout << "listNewickTrees = \n";
//for(set<string> :: iterator itt = listNewickTrees.begin(); itt != listNewickTrees.end(); ++itt )
//{
//cout << *itt << endl;
//}
}



///////////////////////////////////////////////////////////////////////////////////
// helper to figure out the tree specific

DeepCoalescenceHelper :: DeepCoalescenceHelper(const vector<PhylogenyTreeBasic *> &listGTP)
{
	listGTreePtrs = listGTP;
	Init();
}

void DeepCoalescenceHelper:: Init()
{
	// collect info for all leaf labels
	for(int tr=0; tr<(int) listGTreePtrs.size(); ++tr )
	{
		PhylogenyTreeIterator titor( *listGTreePtrs[tr] );
		titor.Init();
		while( titor.IsDone() == false )
		{
			// get current node
			TreeNode *pnCurr = titor.GetCurrNode();
			// if this is leaf, nothing
			//if( pnCurr->IsLeaf() == false )
			//{
			// 
			multiset<int> setLabels;
			CollectSpecies(pnCurr, setLabels);
			// add record
			mapNodesUnder.insert( map<TreeNode *, multiset<int> > :: value_type( pnCurr, setLabels ) ) ;
			//}

			titor.Next();
		}
	}
}

int DeepCoalescenceHelper :: GetExtraLineageNum( const MDC_Cluster &cluster)
{
//cout << "For cluster: ";
//DumpIntSet( cluster );
	int res = 0;
	for( int tr=0; tr<(int) listGTreePtrs.size(); ++tr  )
	{
		int numEL = GetExtraLineageNumFor(tr, cluster);
//cout << "Tree " << tr << ": extra lineage " << numEL << endl;
		res += numEL;
	}
//cout << "listGTreeMultiplicty: ";
//DumpIntVec( listGTreeMultiplicty );
	return res;
}

int DeepCoalescenceHelper :: GetExtraLineageNumFor( int tr, const MDC_Cluster &cluster)
{
//cout << "GetExtraLineageNumFor: tr = " << tr << ", cluster = ";
//DumpIntSet( cluster );
	// compute the extra coalecen for a given cluster of a given gene tree
	// approach: traverse the tree. Find all internal nodes (clusters)
	// with taxa all contained inside clusters
	// then go over these nodes once again, and count the number of such clusters whose 
	// paretn is not in the list
	set< TreeNode *> setNodesContained;
	PhylogenyTreeIterator titor( *listGTreePtrs[tr] );
	titor.Init();
	while( titor.IsDone() == false )
	{
		// get current node
		TreeNode *pnCurr = titor.GetCurrNode();
		// if this is leaf, nothing

		// 
		YW_ASSERT_INFO( mapNodesUnder.find( pnCurr ) != mapNodesUnder.end(), "FATAL ERROR" );
		multiset<int> setLabels = mapNodesUnder[pnCurr];
		set<int> sscur;
		ConvMSetToSet(setLabels, sscur);
		if( IsSetContainer(cluster, sscur) == true )
		{
			// add it
			YW_ASSERT_INFO( pnCurr != NULL, "FATAL ERROR" );
			setNodesContained.insert( pnCurr );
		}
		titor.Next();
	}

	// 
	int res = 0;
	for( set< TreeNode *> :: iterator it = setNodesContained.begin(); it != setNodesContained.end(); ++it)
	{
		TreeNode *ppar = (*it)->GetParent();
		if( setNodesContained.find( ppar) == setNodesContained.end() )
		{
			res++;
		}
	}
//cout << "Extra lineage num = " << res-1 << endl;
	YW_ASSERT_INFO(res >= 1, "FATAL ERROR");
	YW_ASSERT_INFO( listGTreePtrs.size() == listGTreeMultiplicty.size(), "Tree and multiplicty size: mismatch");
	return (res-1)* listGTreeMultiplicty[tr];
}



//void DeepCoalescenceHelper :: CollectSpecies(TreeNode *pn, multiset<int> &listSpecies)
//{
//	// 
//	YW_ASSERT_INFO( pn != NULL, "Can not be empty" );
//	listSpecies.clear();
//	vector<string> listLeafLabels;
//	pn->GetAllLeafLabeles(listLeafLabels);
//	// convert to int
//	for(int i=0; i<(int)listLeafLabels.size(); ++i)
//	{
//		int id;
//		sscanf( listLeafLabels[i].c_str(), "%d", &id );
//		listSpecies.insert(id);
//	}
//}

void DeepCoalescenceHelper :: GetAllClusterFor(int tr, set< set<int> > &clusters)
{
	GetAllClusterForTree( *listGTreePtrs[tr], clusters );
	//PhylogenyTreeIterator titor( *listGTreePtrs[tr] );
	//titor.Init();
	//while( titor.IsDone() == false )
	//{
	//	// get current node
	//	TreeNode *pnCurr = titor.GetCurrNode();
	//	YW_ASSERT_INFO( mapNodesUnder.find( pnCurr ) != mapNodesUnder.end(), "Fail" );
	//	multiset<int> setLabels = mapNodesUnder[pnCurr];
	//	// 
	//	set<int> setDistinctLabels;
	//	ConvMSetToSet(setLabels, setDistinctLabels);
	//	// 
	//	clusters.insert(setDistinctLabels);
	//	titor.Next();
	//}
}


///////////////////////////////////////////////////////////////////////////////////
// main interface class
int DEF_MAX_NUM_NEW_CLUSTER_NUM = 5000;

// need a list of gene trees
DeepCoalescence :: DeepCoalescence(  const vector<PhylogenyTreeBasic *> &listGTP, int nl )
: listGTreePtrs(listGTP),  numLevels(nl), mdcHelper(listGTP), mdcTable( nl ), fIncMoreClades(true), maxNumExtraClades(DEF_MAX_NUM_NEW_CLUSTER_NUM)
{
	set<int> species;
	for(int i=0;i<(int)listGTreePtrs.size(); ++i)
	{
		set<int> speciesStep;
		GetLeavesSpeciesFromGT( *listGTreePtrs[i], speciesStep );

		//vector<string> listLeafLabels;
		//listGTreePtrs[i]->GetAllLeafLabeles(listLeafLabels);
		// fill in the set
		//for(int i=0; i<(int)listLeafLabels.size(); ++i)
		//{
		//	// convert to integer
		//	int id = -1;
		//	sscanf( listLeafLabels[i].c_str(), "%d", &id );
		//	speciesStep.insert(id);
		//}


		//listGTreePtrs[i]->GetLeavesSpeciesGT(speciesStep);
//cout << "SpeicesStep: ";
//DumpIntSet( speciesStep );
		if( species.size() == 0 )
		{
			species = speciesStep;
		}
        if( species != speciesStep )
        {
            string strNW;
            listGTreePtrs[i]->ConsNewickSorted(strNW);
            cout << "Tree " << i << ": this tree seems odd: " << strNW << endl;
        }
		YW_ASSERT_INFO( species == speciesStep, "Fatal error: gene trees species inconsistent" );
	}
	// convert to vector format
	//PopulateVecBySet(setSpecies, species);
	this->numSpecies = species.size();

	// setup default multiplicity
	for(int i=0; i<(int)listGTP.size(); ++i)
	{
		listGTreeMultiplicty.push_back(1);
	}
	mdcHelper.SetMultiplictyofTrees(listGTreeMultiplicty);
}

int DeepCoalescence :: FindMDC()
{
	mdcTable.Init();

	// initiailize the base cases
	for(int i=0; i<numSpecies; ++i)
	{
		// 
		MDC_Cluster singletonCluster;
		singletonCluster.insert( i );
		int numDCCur = mdcHelper.GetExtraLineageNum(singletonCluster);
//cout << "Singleton " << i << ": extra lineage = " << numDCCur << endl;
		//
		set<int> empty;
		mdcTable.AddEntry(singletonCluster, empty, -1, empty, -1, numDCCur, 1);
	}

	// 
	// try all subsets of species
	for(int nr = 2; nr <= numSpecies; ++nr)
	{
//mdcTable.Dump();
//cout << "**************************nr = " << nr << endl;
		// Now, we enumerate all combinations of nr rows
		// We do this by first get the position vector of each such cases
		vector<int> posvec;
		GetFirstCombo( nr, numSpecies, posvec );
		while(true)
		{
//cout << "*****************************************************************\nNow treat cell for: " ;
//DumpIntVec( posvec );
			// Create a vector, and set the index vector
			//vector<int> subsetIndex;
			//vector<int> subsetRows;
			//mgtHelper.GetSpeciesAt( posvec, subsetRows  );
			//GetIntVec(numSpecies, posvec, subsetIndex );
			set<int> curCluster;
			PopulateSetByVec(curCluster, posvec );
			//subsetRows = posvec;
//cout << "subset of species = " ;
//DumpIntSet( curCluster );
//DumpIntVec( subsetRows );

			// this is the number of additional coalescent needed
			int numDCCur = mdcHelper.GetExtraLineageNum(curCluster);
			//int numDCCur = 1;
//cout << "numDCCur = " << numDCCur << endl;
			// now process this subset of species inside it
			// Note: we need to avoid duplicates
			for(int nrsub = 0; nrsub <= nr-2; ++nrsub)
			{
//cout << "**************************nrsub = " << nrsub << endl;
				// Now, we enumerate all combinations of nr rows
				// We do this by first get the position vector of each such cases
				vector<int> posvecsub;
				GetFirstCombo( nrsub, nr-1, posvecsub );
				while(true)
				{
					// get the subset of items (by appending the last one)
					set<int> posvecsubset;
					PopulateSetByVec(posvecsubset, posvecsub);
					vector<int> specurvec;
					GetSubsetVec( posvec, posvecsubset, specurvec );
					specurvec.push_back( posvec[ posvec.size()-1 ] );
					set<int> specur;
					PopulateSetByVec( specur, specurvec );
					set<int> specurRest = curCluster;
					SubtractSets(specurRest, specur);
//cout << "Processing subset: ";
//DumpIntSet( specur );
//cout << "and ";
//DumpIntSet( specurRest );

					// find their record
					for( int ss1 = 0; ss1<numLevels; ++ss1 )
					{
						int valss1 = mdcTable.GetValAt(specur, ss1);
                        if( valss1 < 0 )
                        {
                            continue;
                        }
						int numss1 = mdcTable.GetNumAt(specur, ss1);
						for( int ss2 = 0; ss2<numLevels; ++ss2 )
						{
							int valss2 = mdcTable.GetValAt(specurRest, ss2);
                            if( valss2 < 0 )
                            {
                                continue;
                            }
                            int numss2 = mdcTable.GetNumAt(specurRest, ss2);

							// 
							int valNew = valss1 + valss2 + numDCCur;
							int numNew = numss1 * numss2;
							if( numNew > 0 )
							{
								// add a record
								mdcTable.AddEntry( curCluster, specur, ss1, specurRest, ss2, valNew, numNew );
							}
						}
					}


					if( GetNextCombo( nrsub, nr-1, posvecsub ) == false )
					{
						break;
					}
				}
			}

//cout << "******************************************************\n"; 
//set<int> ssRows;
//PopulateSetByVec( ssRows, subsetRows );
//int numCfgs = storeLinCfgs.GetNumLinCfgsAt( ssRows );
//cout << "Total number of configs = " << numCfgs << endl;
//storeLinCfgs.DumpCfgsAt(ssRows);
//storeLinCfgs.DumpCfgsStatisticsAt(ssRows);
//cout << "******************************************************\n"; 
			// Get the next position
			if( GetNextCombo( nr, numSpecies, posvec ) == false )
			{
				break;
			}
		}

	}

	// not done yet
//	YW_ASSERT_INFO(false, "TBD");

//cout << "**** Final MDC table: \n";
//mdcTable.Dump();

	set<int> spRoot;
	PopulateSetWithInterval(spRoot, 0, numSpecies-1);
	int mdc = mdcTable.GetValAt(spRoot, 0);
	int mdcNum = mdcTable.GetNumAt(spRoot, 0);
//cout << "******************************************************\n";
//	cout << "Minimum deep coalescent is: " << mdc << ", and the number of MDC trees: " << mdcNum << endl;
//cout << "******************************************************\n";
	// dump out sub-optimal solutions
#if 0
	cout << "Suboptimal solutions:\n";
	for(int ii=0; ii<numLevels; ++ii)
	{
		cout << "objective value: " << mdcTable.GetValAt(spRoot, ii) 
			<< ", number of species trees with this objective: " << mdcTable.GetNumAt(spRoot, ii) << endl;
	}
#endif
	//MarginalTree mTreeRes;
	//storeLinCfgs.TraceBackBuildTree(sp, mleST);
//cout << "******************************************************\n"; 
	//cout << "Found optimal tree = ";
	//mleST.Dump();


	return mdc;

}

//int DeepCoalescence :: CountNearOptSTrees(int objVal)
//{
//	// TBD
//	YW_ASSERT_INFO(false, "TBD");
//}

int DeepCoalescence :: RetrieveNearOptSTrees( vector<set<string> > &listNewickTrees, int maxNeeded )
{
	listNewickTrees.clear();
	set<int> spRoot;
	PopulateSetWithInterval(spRoot, 0, numSpecies-1);
	int numTotTree = 0;
	for(int i=0; i<numLevels; ++i)
	{
		set<string> setNewNWs;
		mdcTable.RetrieveNearOptSTrees( spRoot, i, setNewNWs, maxNeeded );
		// add to output
		listNewickTrees.push_back(setNewNWs);
		numTotTree += setNewNWs.size();
		//for(set<string> :: iterator itt = setNewNWs.begin(); itt != setNewNWs.end(); ++itt)
		//{
		//	//YW_ASSERT_INFO( listNewickTrees.find( *itt ) == listNewickTrees.end(), "Can not have duplicate" );
		//	listNewickTrees.push_back(*itt);
		//}
		if( (int)numTotTree >= maxNeeded )
		{
			break;
		}
	}
	return numTotTree;
}



int DeepCoalescence :: FindMDCHeu()
{
	mdcTable.Init();

	// The difference here is that, we only consider clusters formed by 
	// clusters appeared in gene trees. That is, we will not create clusters not present
	// in gene trees
	// initiailize the base cases

	// collect the clusters
	vector< set< set<int> > > listGeneTreeClusters(numSpecies+1);
	ClusterCandidateFinder ccFinder( listGTreePtrs );
	ccFinder.SetMaxCladeNumLimit( maxNumExtraClades );
	const int CLUSTER_LEVEL = 1;
	int levevlToUse = CLUSTER_LEVEL;
	if( fIncMoreClades == false )
	{
		levevlToUse = 0;
	}
	ccFinder.FindCandidates(levevlToUse, listGeneTreeClusters);
	//for(int i=0; i<(int) listGTreePtrs.size(); ++i )
	//{
	//	// traverse the trees
	//	set< set<int> > clusters;
	//	mdcHelper.GetAllClusterFor(i, clusters);
	//
	//	// store in proper slots
	//	for( set< set<int> > :: iterator itt = clusters.begin(); itt != clusters.end(); ++itt )
	//	{
	//		YW_ASSERT_INFO( (int)itt->size() <= numSpecies, "Fail2" );
	//		listGeneTreeClusters[itt->size() ].insert(*itt);
	//	}
	//}


	// initialize the table
	for(int i=0; i<numSpecies; ++i)
	{
		// 
		MDC_Cluster singletonCluster;
		singletonCluster.insert( i );
		int numDCCur = mdcHelper.GetExtraLineageNum(singletonCluster);
//cout << "HEU: Singleton " << i << ": extra lineage = " << numDCCur << endl;
		//
		set<int> empty;
		mdcTable.AddEntry(singletonCluster, empty, -1, empty, -1, numDCCur, 1);
	}

	// 
	// try all subsets of species
	bool fFoundOneSolution = false;
	
	for(int nr = 2; nr <= numSpecies; ++nr)
	{
//mdcTable.Dump();
//cout << "**************************nr = " << nr << endl;
		// Now, we enumerate all combinations of nr rows
		// We do this by first get the position vector of each such cases
		// for each level, at least one item must be able to find matches, if not, then we fail
		for( set<set<int> > :: iterator itts = listGeneTreeClusters[nr].begin(); itts != listGeneTreeClusters[nr].end();  ++itts  )
		{
			// Create a vector, and set the index vector
			//vector<int> subsetIndex;
			//vector<int> subsetRows;
			//mgtHelper.GetSpeciesAt( posvec, subsetRows  );
			//GetIntVec(numSpecies, posvec, subsetIndex );
			set<int> curCluster = *itts;
			//subsetRows = posvec;
//cout << "subset of species = " ;
//DumpIntSet( curCluster );
//DumpIntVec( subsetRows );

			// this is the number of additional coalescent needed
			int numDCCur = mdcHelper.GetExtraLineageNum(curCluster);
//cout << "numDCCur = " << numDCCur << endl;
			// now process this subset of species inside it
			// Note: we need to avoid duplicates
			set< set<int> > doneSets;
			for(int nrsub = 1; nrsub <= nr-1; ++nrsub)
			{
//cout << "**************************nrsub = " << nrsub << endl;
				// Now, we enumerate all combinations of nr rows
				// We do this by first get the position vector of each such cases
				for( set<set<int> > :: iterator ittssub = listGeneTreeClusters[nrsub].begin(); ittssub != listGeneTreeClusters[nrsub].end();  ++ittssub  )
				{
					// get the subset of items (by appending the last one)
					set<int> specur = *ittssub;
					if(  IsSetContainer(curCluster, specur) == false  )
					{
						continue;
					}
					if( doneSets.find( specur ) != doneSets.end() )
					{
						continue;
					}

					set<int> specurRest = curCluster;
					SubtractSets(specurRest, specur);

					// mark both as done to avoid symmetry
					doneSets.insert( specur );
					doneSets.insert( specurRest );

//cout << "Processing subset: ";
//DumpIntSet( specur );
//cout << "and ";
//DumpIntSet( specurRest );

					// is the other part in the list
					if( listGeneTreeClusters[ specurRest.size() ].find( specurRest ) ==  listGeneTreeClusters[ specurRest.size() ].end() )
					{
						continue;
					}

					if( nr == numSpecies )
					{
						fFoundOneSolution = true;
					}

					// find their record
					for( int ss1 = 0; ss1<numLevels; ++ss1 )
					{
						int valss1 = mdcTable.GetValAt(specur, ss1);
                        if( valss1 < 0 )
                        {
                            // not valid, skip
                            continue;
                        }
						int numss1 = mdcTable.GetNumAt(specur, ss1);
						for( int ss2 = 0; ss2<numLevels; ++ss2 )
						{
							int valss2 = mdcTable.GetValAt(specurRest, ss2);
                            if( valss2 < 0 )
                            {
                                continue;
                            }
							int numss2 = mdcTable.GetNumAt(specurRest, ss2);

							// 
							int valNew = valss1 + valss2 + numDCCur;
							int numNew = numss1 * numss2;
							if( numNew > 0 )
							{
								// add a record
								mdcTable.AddEntry( curCluster, specur, ss1, specurRest, ss2, valNew, numNew );
							}
						}
					}

				}
			}

//cout << "******************************************************\n"; 
//set<int> ssRows;
//PopulateSetByVec( ssRows, subsetRows );
//int numCfgs = storeLinCfgs.GetNumLinCfgsAt( ssRows );
//cout << "Total number of configs = " << numCfgs << endl;
//storeLinCfgs.DumpCfgsAt(ssRows);
//storeLinCfgs.DumpCfgsStatisticsAt(ssRows);
//cout << "******************************************************\n"; 

		}


	}

	if( fFoundOneSolution == false )
	{
		return -1;
	}

	// not done yet
//	YW_ASSERT_INFO(false, "TBD");

//cout << "**** Final MDC table: \n";
//mdcTable.Dump();

	set<int> spRoot;
	PopulateSetWithInterval(spRoot, 0, numSpecies-1);
	int mdc = mdcTable.GetValAt(spRoot, 0);
	int mdcNum = mdcTable.GetNumAt(spRoot, 0);
//cout << "******************************************************\n"; 
//	cout << "Minimum deep coalescent is: " << mdc << ", and the number of MDC trees: " << mdcNum << endl;
//cout << "******************************************************\n"; 
	// dump out sub-optimal solutions
#if 0
	cout << "Suboptimal solutions:\n";
	for(int ii=0; ii<numLevels; ++ii)
	{
		cout << "objective value: " << mdcTable.GetValAt(spRoot, ii) 
			<< ", number of species trees with this objective: " << mdcTable.GetNumAt(spRoot, ii) << endl;
	}
#endif

	//MarginalTree mTreeRes;
	//storeLinCfgs.TraceBackBuildTree(sp, mleST);
//cout << "******************************************************\n"; 
	//cout << "Found optimal tree = ";
	//mleST.Dump();


	return mdc;
}


int DeepCoalescence :: GetMDCSolVal(int rank) const
{
	MDC_Cluster rootCluster;
	PopulateSetWithInterval(rootCluster, 0, this->numSpecies-1 );
	return mdcTable.GetValAt( rootCluster, rank );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////


ClusterCandidateFinder :: ClusterCandidateFinder(const vector<PhylogenyTreeBasic *> &listGTP) : listGTreePtrs(listGTP),
  maxNumExtraClades( DEF_MAX_NUM_NEW_CLUSTER_NUM )
{
}


void ClusterCandidateFinder :: FindCandidates(int level, vector< set< set<int> > > &clustersPerSize)
{
	clustersPerSize.clear();

	set<set<int> > clustersCurr;

	// first get the default clusters that appear exactly in the tree
	for(int i=0; i<(int) listGTreePtrs.size(); ++i )
	{
		// traverse the trees
		set< set<int> > clusters;
		GetAllClusterForTree(*listGTreePtrs[i], clusters);

		// store in proper slots
		for( set< set<int> > :: iterator itt = clusters.begin(); itt != clusters.end(); ++itt )
		{
			clustersCurr.insert(*itt);
		}
	}
//cout << "The number of tree clusters: " << clustersCurr.size() << endl;

	// now, we go to find more possible candidates
	for(int round = 0; round < level; ++round)
	{
		//
		set<set<int> > clustersMore;
		FindRelatedClusters( clustersCurr, clustersMore );

		// add it
		for( set< set<int> > :: iterator itt = clustersMore.begin(); itt != clustersMore.end(); ++itt )
		{
			clustersCurr.insert(*itt);
		}
	}

	// finally store the clusters
	for( set< set<int> > :: iterator itt = clustersCurr.begin(); itt != clustersCurr.end(); ++itt )
	{
		//YW_ASSERT_INFO( (int)itt->size() <= numSpecies, "Fail2" );
		clustersPerSize[itt->size() ].insert(*itt);
	}
//cout << "Total number of found clusters: " << clustersCurr.size() << endl;
}


void ClusterCandidateFinder :: FindRelatedClusters(const set<set<int> > &clustersCurr, set<set<int> > &clustersMore )
{
	// finding more clusters through set merge or set difference of two current sets
	clustersMore.clear();

	// keep track of how often a cluster occur
	map< set<int>, int> mapClusterCount;

	for( set<set<int> > :: iterator it1 = clustersCurr.begin(); it1 != clustersCurr.end(); ++it1 )
	{
		set<set<int> > :: iterator it2 = it1;
		++it2;
		for(; it2 != clustersCurr.end(); ++it2)
		{
			// does one contain another?
			if( IsSetContainer(*it1, *it2) == true  )
			{
//cout << "set1: ";
//DumpIntSet( *it1 );
//cout << "set2: ";
//DumpIntSet( *it2 );

				// 
				set<int> sdiff = *it1;
				SubtractSets( sdiff, *it2);
				YW_ASSERT_INFO( sdiff.size() > 0, "Fail" );
				//clustersMore.insert(sdiff);
//cout << "sdiff1: ";
//DumpIntSet( sdiff );
			}
			else if( IsSetContainer(*it2, *it1) == true  )
			{
//cout << "set1: ";
//DumpIntSet( *it1 );
//cout << "set2: ";
//DumpIntSet( *it2 );
				// 
				set<int> sdiff = *it2;
				SubtractSets( sdiff, *it1);
				//clustersMore.insert(sdiff);
				YW_ASSERT_INFO( sdiff.size() > 0, "Fail" );
//cout << "sdiff2: ";
//DumpIntSet( sdiff );
			}
			else
			{
				// are two disjoint?
				set<int> sint;
				JoinSets( *it1, *it2, sint);
				if( sint.size() == 0)
				{
//cout << "set1: ";
//DumpIntSet( *it1 );
//cout << "set2: ";
//DumpIntSet( *it2 );
					// contact
					set<int> sadd = *it1;
					UnionSets( sadd, *it2);
					//clustersMore.insert(sadd);
					// add to temp list
					if( clustersCurr.find( sadd ) == clustersCurr.end() )
					{
						if(mapClusterCount.find( sadd) == mapClusterCount.end() )
						{
							mapClusterCount.insert( map< set<int>, int> :: value_type(sadd, 1) );
						}
						else
						{
							mapClusterCount[sadd]++;
						}
					}


//cout << "sadd: ";
//DumpIntSet( sadd );

				}
			}
		}
	}


	// now find the largest items
	vector<int> listCounts;
	for( map< set<int>, int> :: iterator itt = mapClusterCount.begin(); itt != mapClusterCount.end(); ++itt )
	{
		listCounts.push_back( itt->second );
	}
	SortIntVec(listCounts);
	int cutoffVal = 0;
	if( mapClusterCount.size() > maxNumExtraClades )
	{
		cutoffVal = listCounts[ (int) mapClusterCount.size() - maxNumExtraClades ];
	}
//cout << "listCounts: ";
//DumpIntVec( listCounts );
//cout << "cutoff = " << cutoffVal << ", maxNumExtraClades  = " << maxNumExtraClades << ", mapClusterCount.size = " << mapClusterCount.size() << endl;

	// start putting it in
	// do it twice
	for( map< set<int>, int> :: iterator itt = mapClusterCount.begin(); itt != mapClusterCount.end(); ++itt )
	{
		if( itt->second > cutoffVal && clustersMore.size() <  maxNumExtraClades )
		{
			clustersMore.insert( itt->first );
		}
	}
	for( map< set<int>, int> :: iterator itt = mapClusterCount.begin(); itt != mapClusterCount.end(); ++itt )
	{
		if( itt->second == cutoffVal && clustersMore.size() <  maxNumExtraClades )
		{
			clustersMore.insert( itt->first );
		}
	}
}


