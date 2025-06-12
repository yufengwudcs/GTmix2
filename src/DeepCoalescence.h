#ifndef DEEP_COALESCENCE_H
#define DEEP_COALESCENCE_H

//#include "MarginalTree.h"
#include "PhylogenyTreeBasic.h"
#include <cmath>


// used for finding the optimal tree minimizing deep coalescence
// also can reconstruct all optimal (or sub-optimal) trees

typedef set<int> MDC_Cluster;
//typedef struct
//{
//	MDC_Cluster src1;
//	MDC_Cluster src2;
//	int level1;
//	int level2;
//} MDCEntrySrc;
typedef pair< pair<MDC_Cluster,int>, pair<MDC_Cluster,int> >  MDCEntrySrc;



class MDCTableEntry  
{
public:
	MDCTableEntry();

	// access
	void Dump() ;
	int GetNumItems() const { return mapNearpOptInfo.size(); }
	//void CheckSol(int val);
	void AddItem(int levelLimit, const MDC_Cluster &clusterSrc1, int ss1, const MDC_Cluster &clusterSrc2, int ss2, int val, int cnt) ;
	int GetVal(int i) const;
	int GetCount(int i) const;
	int GetSrcOptionsNum(int i) const;
	void GetSrcOptionsAt(int i, set<MDCEntrySrc> &subsetOptions);
	//{  
	//	subsetOptions = listSrcClusters[i];
	//}

	//void GetSrcOptionsAt(int i, int j, MDC_Cluster &subset1, MDC_Cluster &subset2) 
	//{  
	//	YW_ASSERT_INFO( j < GetSrcOptionsNum(i), "FATAL ERROR"  );  
	//	subset1 = listSrcClusters[i][j].first;
	//	subset2 = listSrcClusters[i][j].second;
	//}

private:
	void Init();

	map<int, pair<int, set<MDCEntrySrc> > > mapNearpOptInfo;
	bool fDirty;
	vector<int> vecOrderedNearOptVals;

	//vector<int> valsNearOpt;
	//vector<int> counstNearOpt;
	//vector< set<MDCEntrySrc> > listSrcClusters;
};

// this is the main DP table
class MDCTable
{
public:
	MDCTable(int numLevels );
	void Dump();
	void Init() {mdcTable.clear(); }

	void AddEntry( const MDC_Cluster &cluster, const MDC_Cluster &clusterSrc1, int ss1, const MDC_Cluster &clusterSrc2, int ss2, int val, int num );
	int GetValAt(const MDC_Cluster &cluster, int rank ) const;
	int GetNumAt(const MDC_Cluster &cluster, int rank );
	void RetrieveNearOptSTrees( const MDC_Cluster &cluster, int mdcLevel, set<string> &listNewickTrees, int maxNeeded );

private:
	map< MDC_Cluster, MDCTableEntry >  mdcTable;
	int numLevels;
};

// helper to figure out the tree specific
class DeepCoalescenceHelper
{
public:
	DeepCoalescenceHelper(const vector<PhylogenyTreeBasic *> &listGTreePtrs);

	void Init();
	int GetExtraLineageNum( const MDC_Cluster &cluster);
	int GetExtraLineageNumFor( int tr, const MDC_Cluster &cluster);
	void GetAllClusterFor(int tr, set< set<int> > &clusters);
	void SetMultiplictyofTrees(const vector<int> &vecMultis) { listGTreeMultiplicty = vecMultis; }

private:
	//void CollectSpecies(TreeNode *pn, multiset<int> &listSpecies);

	vector<PhylogenyTreeBasic *> listGTreePtrs;
	vector<int> listGTreeMultiplicty;
	map<TreeNode *, multiset<int> > mapNodesUnder;
};


// main interface class
class DeepCoalescence
{
public:
	// need a list of gene trees
	DeepCoalescence( const vector<PhylogenyTreeBasic *> &listGTreePtrs, int numLevels );
	void SetMultiplictyofTrees(const vector<int> &vecMultis) { listGTreeMultiplicty = vecMultis; mdcHelper.SetMultiplictyofTrees(vecMultis);}
	void SetIncMoreCladFlag(bool f) { fIncMoreClades = f; }
	void SetMaxCladeNumLimit(int limitNew) { maxNumExtraClades = limitNew; }

	int FindMDC();
	int FindMDCHeu();
	//int CountNearOptSTrees(int objVal);
	int RetrieveNearOptSTrees( vector<set<string> > &listNewickTrees, int maxNeeded = HAP_MAX_INT);
	int GetMDCSolVal(int rank) const;

private:
	// 
	void InsMDCEntry(const MDC_Cluster &item, int val);

	vector<PhylogenyTreeBasic *> listGTreePtrs;
	vector<int> listGTreeMultiplicty;
	int numLevels;		// how many sub-optimal you want to keep
	DeepCoalescenceHelper mdcHelper;
	MDCTable mdcTable;
	int numSpecies;
	bool fIncMoreClades;
	int maxNumExtraClades;
};



///////////////////////////////////////////////////////////////////////////////////////////////////

class ClusterCandidateFinder
{
public:
	ClusterCandidateFinder(const vector<PhylogenyTreeBasic *> &listGTreePtrs);

	void FindCandidates(int level, vector< set< set<int> > > &clustersPerSize);
	void SetMaxCladeNumLimit(int limitNew) { maxNumExtraClades = limitNew; }

private:
	void FindRelatedClusters(const set<set<int> > &clustersCurr, set<set<int> > &clustersMore );

	vector<PhylogenyTreeBasic *> listGTreePtrs;
	int maxNumExtraClades;
	//DeepCoalescenceHelper &mdcHelper;

};


#endif
