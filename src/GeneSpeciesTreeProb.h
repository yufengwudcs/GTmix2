#ifndef GENE_SPECIES_TREE_PROB_H
#define GENE_SPECIES_TREE_PROB_H

#include "MarginalTree.h"
#include "PhylogenyTreeBasic.h"
#include <cmath>

// Note: species tree is binary, while gene trees can be non-binary
// also, gene trees may contain multiple alleles per taxon
// while of course there are only distinct taxa species tree
int GetNumSpeciesFromGT( PhylogenyTreeBasic &treeGene);
int GetNumSpeciesFromGTs( const vector<PhylogenyTreeBasic *> &listGTreePtrs, set<int> *ptrSpecies=NULL );
void GetLeavesSpeciesFromGT( PhylogenyTreeBasic &treeGene, set<int> &species );
bool AreGeneTreesWithSameTaxa( const vector<PhylogenyTreeBasic *> &listGTP );
void GetSpeciesAtNode( TreeNode *pnode, set<int> &speciesUnder );
void SetMaxLineageCfgKept(int maxn);	// set the maximum number of cfg to keep
void SetMaxLineageAncesCfgKept(int maxn);  // set maximum ancestral cfgs to keep (different from the above, this limits the number of found ancestral cfgs during ancestral cfgs construction)
void SetVerbose(bool f);
bool IsVerbose();
bool GetSTELLSProbComputeWarningFlag();
void SetFixedCoalCalcMode(bool f);
bool IsFixedCoalMode();
void TestCalcProbBranch();

//////////////////////////////////////////////////////////////////////////////////
// useful data structure

// lineage cluster corresponds to a node in gene tree, and contain a number of lineages
class LineageCluster
{
public:
	LineageCluster() {}
	LineageCluster( int gnid, int gnidpar );
    LineageCluster(const LineageCluster &rhs) {gnid = rhs.gnid; numLins = rhs.numLins; gnidPar = rhs.gnidPar;}
    LineageCluster&  operator=(const LineageCluster &rhs) {gnid = rhs.gnid; numLins = rhs.numLins; gnidPar = rhs.gnidPar;  return *this;}
	bool operator <(const LineageCluster &lc2) const
	{
		if( gnid < lc2.gnid)
		{
			return true;
		}
		if( gnid == lc2.gnid && numLins < lc2.numLins)
		{
			return true;
		}
		if( gnid == lc2.gnid && numLins == lc2.numLins &&  gnidPar < lc2.gnidPar)
		{
			return true;
		}

		return false;
	}
    bool operator==(const LineageCluster &rhs) const
	{
		if( gnid == rhs.gnid && numLins == rhs.numLins && gnidPar == rhs.gnidPar)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	void Dump() const;

	// note: for the same idNode, we can NOT have two LC with identical idNodePar
	// since they would have merged already
	void Init(int idNode, int nl, int idNodePar) { gnid = idNode; numLins = nl; gnidPar = idNodePar; }
	void Append(const LineageCluster &lc);
	int GetGNId() const {return gnid;}
	void SetGNId(int id) { gnid = id; }
	int GetNumLineages() const {return numLins;}
	int GetGNParId() const  {return gnidPar;}
	void SetGNParId(int par)  { gnidPar = par; }
	void DecByOne() { YW_ASSERT_INFO( numLins > 1, "Can not decrease" ); numLins --; }

private:
	int gnid;
	int numLins;
	int gnidPar;
};

// configuration is a set of LCs, which is associated with a node in ST
class GeneSpeciesTreeHelper;

class LineageConfig
{
public:
	LineageConfig();
	LineageConfig(const LineageConfig &rhs)
	{
		// just call =
		snid = rhs.snid;
		//setLCs = lc2.setLCs;
		setLCs.clear();
		for(  set< LineageCluster > :: iterator itt = rhs.setLCs.begin(); itt != rhs.setLCs.end(); ++itt )
		{
			setLCs.insert(*itt);
		}
		valProb = rhs.valProb;
		numTotLins = rhs.numTotLins;		;
	}
	bool operator <(const LineageConfig &lc2) const;
	void operator=(const LineageConfig &lc2)
	{
		snid = lc2.snid;
		//setLCs = lc2.setLCs;
		setLCs.clear();
		for(  set< LineageCluster > :: iterator itt = lc2.setLCs.begin(); itt != lc2.setLCs.end(); ++itt )
		{
			setLCs.insert(*itt);
		}
		valProb = lc2.valProb;
		numTotLins = lc2.numTotLins;
	}
	bool operator==(const LineageConfig &lc2) const
	{
		return ((*this < lc2) == false && (lc2 < *this) == false);
	}
	bool AreLCsSame(const LineageConfig &lc2) const { return setLCs == lc2.setLCs; }

	// this initialize the config with a single lineage of GT
	void Clear() { snid = -1; setLCs.clear(); RsetLinNum(); }
	void SetNodeId(int id) { snid = id; }
	int GetNodeId() const { return snid; }
	void GetLCs( set< LineageCluster > &setLCsRes ) const { setLCsRes = setLCs; }
	void AddLC(const LineageCluster &lc) {  setLCs.insert(lc); RsetLinNum(); }
	void Consolidate( const GeneSpeciesTreeHelper &gstHelper );
	void Dump() const;
	void CopyReduceOne( const LineageCluster &lcToReduce, LineageConfig &lcNew, int gtnodeNew=-1, int gtnodeParNew=-1 ) const;
	double GetProb() const { return valProb; }
	void SetProb(double prob) { valProb = prob; }
	int GetTotNumLins() const;
	void GetSrcNodes(set<int> &nodeIds) const;
	void GetSrcNodesNumLins( map<int,int> &mapLinInfo) const;
	bool IsAncestralTo( const LineageConfig &cfgDest, const GeneSpeciesTreeHelper &gstHelper ) const;
	void FindCoalescableLins( set< set<LineageCluster> > &setCoalLins ) const;
	void MergeLCInto(const GeneSpeciesTreeHelper &gstHelper, const LineageCluster &lin1, const LineageCluster &lin2, LineageConfig &lcNew ) const;
	void RemoveLineage(const GeneSpeciesTreeHelper &gstHelper,   const LineageCluster &lin1, LineageConfig &lcNew ) const;
	void GetParLins(multiset<int> &parids) const;

private:
	void RsetLinNum() {numTotLins = -1;}
	void SetTotNumLins(int nl) { numTotLins = nl; }

	int snid;		// which node this config is for (meaning these lineages survive up to this node) in SPECIES tree
	set< LineageCluster > setLCs;
	double valProb;
	int numTotLins;
};


// Helper class for determinting constraints imposed by
// species tree about which transitions are allowed
class CoalEvent
{
public:
	// by default, there are two lineages that can coalesce (like that in binary trees)
	CoalEvent(int nid) : nodeId(nid), numLinsToChooseFrom(2) {}
	int GetID() const { return nodeId; }
	int GetCoalLinChoiceNum() const { return numLinsToChooseFrom; }
	void SetCoalLinChoiceNum(int t) { numLinsToChooseFrom = t; }
	void Dump() const { cout << "event: " << nodeId << "  "; }

private:
	int nodeId;
	int numLinsToChooseFrom;		// at mulfurcating nodes, we may have more than wwo lineages that can coalesce
};

struct LTLinCfg
{
  bool operator()(const LineageConfig &lin1, const LineageConfig &lin2) const
  {
    return lin1 < lin2;
  }
};

// used to represent ancestral configs: what events have occured since its source config
struct AncLinConfig
{
public:
	AncLinConfig() { pSrcConfig = NULL; }
	AncLinConfig( const LineageConfig *pscfg,  const LineageConfig &linCfg, const vector<CoalEvent>& evts ) : 
		listEvts( evts), cfgCurr(linCfg), pSrcConfig(pscfg)   { }
	double GetProb() const {return cfgCurr.GetProb();}
	void SetProb(double val) {cfgCurr.SetProb( val );}
	void AddCoalEvent(const CoalEvent &evt) { listEvts.push_back(evt); }
	LineageConfig &GetLinCfg() {return cfgCurr; }
	const LineageConfig *GetSrcLinCfg() { return pSrcConfig; }
	void UpdateProb(double len, GeneSpeciesTreeHelper *gstHelper);
	void AppendAC(const AncLinConfig &cfgAdd);
	void GetCoalEvts(vector<CoalEvent> &evtsCoal) const { evtsCoal = listEvts; }
	void SetCoalEvts(const vector<CoalEvent> &evtsCoal) { listEvts = evtsCoal; }
	void Dump() const;

private:
	vector<CoalEvent> listEvts;
	LineageConfig cfgCurr;
	//double valProb;
	const LineageConfig *pSrcConfig;
};

// useful funcion
int GetSpeciesFromLabel(const string &label);		// what this is for a label?
double UtilsCalcCoalProbBranch( int u, int v, double len );	// prob of coalsecing from u lineage to v lineages with len time
double UtilsCalcCoalProbBranch2( int u, int v, double len );    // prob of coalsecing from u lineage to v lineages with len time
double UtilsCalcCoalProbBranchPoissonApprox( int u, int v, double len );	// prob of coalsecing from u lineage to v lineages with len time
double UtilsCalcCoalProbBranchPoissonApprox2( int u, int v, double len );    // prob of coalsecing from u lineage to v lineages with len time
double UtilsCalcCoalProbBranchNormalApprox( int u, int v, double len );	// prob of coalsecing from u lineage to v lineages with len time

// utility class for mulfurcating tree prob
// this class handle the task of accessing coalescent coefficient w/ non-binary gene tree
class MultiLineageHelper
{
public:
	static MultiLineageHelper & Instance() { return instance; }
	double GetMultiLinsCoeff(int numLinsSrc, int numLinsDest);

private:
	MultiLineageHelper( );

	static MultiLineageHelper instance;
	map<pair<int,int>, double> mapPreCMultiLinsCoeff;
};




// utility clas for gene tree prob
class GeneSpeciesTreeHelper
{
public:
	GeneSpeciesTreeHelper( MarginalTree &treeSpecies, PhylogenyTreeBasic &treeGene );

	void Init();
	void Reset() {cacheBranchProb.clear(); cacheLCAncs.clear();}
	void ClearBranchLenCache() {  cacheBranchProb.clear(); }

	// support query
	void GetAllSubtreeST( set<int> &snodes );		// get all the nodes in species tree
	int GetTaxaAtLeafST(int snid);					// get the taxa id from a particular leaf node of species tree
	void GetLeavesSpeciesGT( set<int> &geneAlleles );		// get all the nodes in species tree
	void GetGeneAllelesForSpecies( int taxon, set<int> &geneAlleles);	// get all alleles in GT for a gene
    void GetGeneAllelesForSTNode( int stNode, set<int> &geneAlleles, const map<int,int> *pMapOldLeavesToNew=NULL );  // get all alleles from species under a specific ST branch
	void FormLinConfig( int snid, const set<int> &linsGT, LineageConfig &linConfig ) const;
	void FormLinConfig( int snid, const LineageConfig &linsGT1, const LineageConfig &linsGT2, LineageConfig &linConfig ) const;
	//void FindAncestralCoalConfig( const LineageConfig &linConfigInit, set< LineageConfig > &setAncesLinConfigs, set<LineageConfig> &cacheDoneLins );
	int GetMRCATwoLineagesGT( int lin1, int lin2 ) const;
	int GetParentLineage(int lin) const;
	double CalcCoalCoeffForTwoCfg( const LineageConfig &linCfg1, const LineageConfig &linCfg2 ) const;
	int GetRootIdGT();
	double CalcCoalProbBranch( int u, int v, double len );	// prob of coalsecing from u lineage to v lineages with len time
	//void CollectAncestralCfgs( int nodeId, const set<LineageConfig,LTLinCfg> &setOrigCfgs, set<LineageConfig,LTLinCfg> &setAncestralCfgs  );
	//void CollectAncestralCfgs( double branchLen, const set<LineageConfig,LTLinCfg> &setOrigCfgs, set<LineageConfig,LTLinCfg> &setAncestralCfgs  );
	bool IsNodeAncestralToGT(int anc, int des) const;
	double CalcdbVal(int ub, int cb) const;
    double CalcProddbValAndCoalFac(int ub, int cb) const;
    double CalcLogProddbValAndCoalFac(int ub, int cb) const;
	double CalcwbVal( const vector<CoalEvent> &listEvts, const LineageConfig &lcSrc, int ub = -1, int cb = -1 ) const;
	double CalcwbValMulti( const LineageConfig &linCfgSrc, const LineageConfig &lcDest );
    double CalcwbValMultiHeu( const vector<CoalEvent> &listEvts, const LineageConfig &lcSrc, int ub = -1, int cb = -1 ) const;

	// species tree/gene tree related
    int GetSpeciesForGTLeaf(int lin) const;     // for a gene tree leaf, find the species (caution: GetSpeciesForGTLin does not work)
	int GetSpeciesForGTLin(int lin) const;		// for a GT lineage, what species it comes from; YW: this does not work yet
	void SetupSpeciesIdForGTLin( const vector< pair<int,int> > &listRecords);
	//int GetDescendentsNumGT(int nodeGT) const;
	void GetGTNodeDescendents( int nodeId, set<int> &listDescs) const;
	int GetGTNodeNumChildren(int gnid) const;

	// new functions: hopefully faster
	void FastCollectAncConfig( int nodeId, const set<LineageConfig,LTLinCfg> &setOrigCfgs, vector< AncLinConfig *> &listAncCfgsNew );
	void FastCollectAncConfig( double lenOrig, const set<LineageConfig,LTLinCfg> &setOrigCfgs, vector< AncLinConfig *> &listAncCfgsNew );
	//string FindGNShapeLabel( const set<int> &ids );
	string FindGNShapeLabel( const LineageConfig &cfgCurr );
	bool IsGeneTreeMulfurcate() const {  return fMulfurcateTree; }
    void GetNodesForIds(const set<int> &nids, set<TreeNode *> &setNodes) const;

private:
	// for quick query
	int GetNodeGTPar(int gnid) const;
	void FindDoneLineages( const set<int> &nodeIdsSrc, const map<int,int> &srcLinNums, 
		const set<int> &nodeIdsDest, const map<int,int> &destLinNums, set< pair<int,int> > &nodesDone ) const;
	void FindCoalEventsBtwConfigs( const LineageConfig &linCfg1, const LineageConfig &linCfg2, vector<CoalEvent> &listEvts  ) const;
	//int FindNumLinsFromCfgToNode(const set<int> &nodeIdsSrc, int gtNodeId) const;
	double CountParFromChildCoal( const vector< pair<int,double> > &listChildCounts, int numRemLins) const;
	bool AreSrcDestCfgsCompatible( const LineageConfig &linCfgSrc, const LineageConfig &linCfgDest );
	void AddToCache( const LineageConfig &setOrigCfgs, const set<LineageConfig> &setAncestralCfgs);
	bool AreLinsHalfSibling( const LineageCluster &lin1, const LineageCluster &lin2 ) const;
	void FastCollACSingleCfg( double lenBranch, const LineageConfig &origCfgs, 
		vector< AncLinConfig *> &listAncCfgsNew, map<string, AncLinConfig*> &mapExistACs );
	void AddCoalEvt( vector<CoalEvent> &evtsCoal, const LineageCluster &lin1, const LineageCluster &lin2 ) const;
	void InitDefaultGTLinMap();
	double CalcBranchFactorAt(int numDupEvts, const set<int> &sDescs, int numActiveLins, map<int,int> &mapNodesUnder );
	//double CalcBranchFactorNoDesc(int numDupEvts);
	//void InitMultiNoDescCoeff();
	double GetMultiNoDescCoeff(int numLinsOrig, int numDupEvts);
	void FindNumLinsForEvtsCfg(const vector<int> &listEvts, map<int,int> &mapEvtDupNums, 
		const LineageConfig &linCfgSrc, vector<int> &listLinsNumForEvt) const;
	int GetNumLinsAtNode(int gnid) const;
	double CalcMultiBranchComboNum( const vector< pair<int,double> > &listCounts ) const;
	void GetLinCfgParLinsSet(const LineageConfig &linCfgSrc, multiset<int> &sint) const;

	// define how the node counting work
	typedef struct
	{
		vector<int> lineagesDesc;		// which lineages are under it
		int numEvtsBelow;				// number of events under this lineage (including those right at this node)
		double numWays;					// number of legal ways to arrange these events
	} NODE_COUNT_INFO;


	MarginalTree &treeSpecies;
	PhylogenyTreeBasic &treeGene;
	bool fMulfurcateTree;

	// mapping from GT lin id to species id
	map<int,int> mapGTLinToSpecies;

	// data structure for enforcing the coalescent constraints by ST/GT
	// The most important thing is to decide which nodes can merge and when
	// For this, we create a data structure called lineage cluster (which represent several lineages
	// entering the same parent node in gene tree).

	// for speedup purpose
	vector< vector< bool> > vecNodeAncestralGT;	// [i][j] = true if node i is ancestral to j in GT
	//vector< vector< int> > vecMRCANodePairGT;	// [i][j] = MRCA node of two nodes i and j in GT
	map<int, int> mapNodeParGT;					// in GT, which parent node this node has?
	map<int, TreeNode*> mapIdNode;				// in GT, which id is for which TreeNode
	map< int, int > mapTaxonLeafNodeIdST;		// for a given taxon, which leaf node in ST
	map<int, int> mapNodeParST;					// in ST, which parent node this node has?
	set< int > setLeafIdGT;						// in GT, what are the ids for leaves?
	vector< vector<int> > vecMRCATwoLinsGT;		// for two lineages a and b, their MRCA
	//map<int,int> mapNodesDescendentNumGT;		// in GT, for a node, how many descendent it has
	map<int, set<int> > mapNodesDescendentsGT;	// in GT, for a node, what descendents it has
	vector<double> vecMultiNoDescCoeff;			// this is a table showing how to compute no-descendent multiple lineage coalesceing coefficient (similar to Wb factor in Degnan's paper)
	//int szMultiNoDescAllowed;

	// for speedup purpose
	map< pair<pair<int,int>,double>, double > cacheBranchProb; 
	map<LineageConfig, set< LineageConfig > > cacheLCAncs;
	map< pair<multiset<int>, multiset<int> >, double> cacheWbValMulti;
};

// store the LineageConfigs
class LineageConfigSore;
//typedef map<int, set< LineageConfig,LTLinCfg > >  MAP_NUMLINS_LINS;
//typedef map<int, MAP_NUMLINS_LINS >  MAP_ID_LINS;


class LinCfgSoreIterator
{
public:
	LinCfgSoreIterator( int nodeId, LineageConfigSore &lcStore );
	void Init();
	bool IsDone();
	void Next();
	const LineageConfig &GetLinCfg() { return *it; }

private:
	int nodeId;
	LineageConfigSore &lcStore;
	set< LineageConfig > :: iterator it;
};

#if 0
class LinCfgSoreIterator
{
public:
	LinCfgSoreIterator( int nodeId, LineageConfigSore &lcStore );
	void Init();
	bool IsDone();
	void Next();
	const LineageConfig &GetLinCfg() { return *it; }

private:
	int nodeId;
	LineageConfigSore &lcStore;
	set< LineageConfig, LTLinCfg> :: iterator it;
	MAP_NUMLINS_LINS :: iterator itMap;
};
#endif

class LineageConfigSore
{
	friend class LinCfgSoreIterator;
public:
	LineageConfigSore(MarginalTree &treeSpecies, PhylogenyTreeBasic &treeGene, GeneSpeciesTreeHelper &gstHelper);

	// add a record of LinCfg: for the purpose of computing the probablity, need to know which
	// LinCfg it derives from. For leaves, we also support a plain version: with no probility
	void AddLinCfg( int nodeId, LineageConfig &lcNew );
	//void AddLinCfg( int nodeId, LineageConfig &lcNew, const LineageConfig &lcLeftAnc, const LineageConfig &lcLeftSrc,
	//	const LineageConfig &lcRightAnc, const LineageConfig &lcRightSrc );
	int GetNumLinCfgsAt(int nodeId);
	double CalcTotProbAt(int nodeId);
	void DumpConfigsAt(int nodeId);
	//void GetLCSetAt(int node, set< LineageConfig,LTLinCfg > &res)
	//{
	//	YW_ASSERT_INFO( statesDepot.find(node) != statesDepot.end(), "Fail to find" );
	//	res = statesDepot[node];
	//}
	set< LineageConfig,LTLinCfg > & GetLCSetAt(int node)
	{
//cout << "GetLCSetAt: node: " << node << endl;
        if( statesDepot.find(node) == statesDepot.end() )
        {
            cout << "FATAL ERROR: this node is not found in statesDepot: " << node  << ". Possibly caused by numerical underflow."<< endl;
        }
		YW_ASSERT_INFO( statesDepot.find(node) != statesDepot.end(), "Fail to find1.2" );
		return statesDepot[node];
	}
	void SetLCSetAt(int node, set< LineageConfig,LTLinCfg > &setCfgs)
	{
		if( statesDepot.find(node) == statesDepot.end() )
		{
			statesDepot.insert( map<int, set< LineageConfig,LTLinCfg > > :: value_type( node, setCfgs )  );
		}
		else
		{
			statesDepot[node] = setCfgs;
		}
	}
	void SetLinCfgProbAt(int nodeId, LineageConfig &lcNew, double probNew);
	void Reset() { statesDepot.clear();  }
	void ResetAt( int node ) 
	{ 
		statesDepot.erase(node);  
		YW_ASSERT_INFO( statesDepot.find(node) == statesDepot.end(), "Error in ResetAt" );
	}

private:
	//LineageConfig *FindLinCfgRep(int nid, const LineageConfig &linCfg);
	void Prune(int nodeId);
	//string FindGNShapeLabel( const set<int> &ids );


	MarginalTree &treeSpecies;
	PhylogenyTreeBasic &treeGene;
	GeneSpeciesTreeHelper &gstHelper;

	map<int, set< LineageConfig,LTLinCfg > > statesDepot;

//	MAP_ID_LINS statesDepot;
};

//////////////////////////////////////////////////////////////////////////////////
// Base class for prob compuatation
    
class GenericGeneSpeciesTreeProb
{
public:
  virtual ~GenericGeneSpeciesTreeProb() {}
  virtual double CalcProb() = 0;
  virtual double CalcLogProb() = 0;
  // attempt to change a branch to a new length, and update the prob. Since we may need to recover the original
  // unchanged probs, we just save those changed one (i.e. keep a copy of what was there before)
  virtual double TestNewBranchLenAt(int threadId, int branch, double lenNew, map<int,set< LineageConfig,LTLinCfg > > &origLinCfgs, bool fSetBrlen ) = 0;
  virtual void SetLinCfgs(map<int,set< LineageConfig,LTLinCfg > > &linCfgsToSet) = 0;	// useful when we want to un-do the effect of the previous tweaking of branch length
  virtual double GetSTBrachLen(int br) const = 0;
  virtual void SetSTBranchLen(int br, double brLen) = 0;
  virtual void GetSpeciesTree(MarginalTree &spTree) const = 0;
  virtual int GetTaxaAtLeafST(int snid) = 0;
  virtual void GetGeneAllelesForSpecies( int taxon, set<int> &geneAlleles) = 0;
    virtual void JumpStartAtSTNodes( GenericGeneSpeciesTreeProb *pProbPrev, const set<int> &setSTNodesToInit, map<int,int> &mapSTNewToPrevST )
    {
        // by default, do nothing
    }
    virtual void SetThreadId(int tid)
    {
        // by default, do nothing
    }
    virtual std::string GetSpeciesTreeNW() const = 0;
};
      
// simply return prob 1.0 (or log 0.0); used as a place holder
class DummyGeneSpeciesTreeProb : public GenericGeneSpeciesTreeProb
{
public:
  DummyGeneSpeciesTreeProb() {}
  virtual ~DummyGeneSpeciesTreeProb() {}
  virtual double CalcProb() {return 1.0;}
  virtual double CalcLogProb() {return 0.0; }
  // attempt to change a branch to a new length, and update the prob. Since we may need to recover the original
  // unchanged probs, we just save those changed one (i.e. keep a copy of what was there before)
  virtual double TestNewBranchLenAt(int threadId, int branch, double lenNew, map<int,set< LineageConfig,LTLinCfg > > &origLinCfgs, bool fSetBrlen ) { return CalcLogProb(); }
  virtual void SetLinCfgs(map<int,set< LineageConfig,LTLinCfg > > &linCfgsToSet) {}
  virtual double GetSTBrachLen(int br) const { return 0.0; }
  virtual void SetSTBranchLen(int br, double brLen) {}
  virtual void GetSpeciesTree(MarginalTree &spTree) const {}
  virtual int GetTaxaAtLeafST(int snid) {return 0; }
  virtual void GetGeneAllelesForSpecies( int taxon, set<int> &geneAlleles) {}
    virtual std::string GetSpeciesTreeNW() const { return std::string(""); }
};
      
      
      
// main class for computing a given gene tree on species tree probability

class GeneSpeciesTreeProb : public GenericGeneSpeciesTreeProb
{
public:
	// need two trees
	GeneSpeciesTreeProb(MarginalTree &treeSpecies, PhylogenyTreeBasic &treeGene);
    virtual ~GeneSpeciesTreeProb();
	//void ReadGTLinSpecies( const char *filename  );
	virtual double CalcProb();
    virtual double CalcLogProb() { return log( CalcProb() ); }
	// attempt to change a branch to a new length, and update the prob. Since we may need to recover the original
	// unchanged probs, we just save those changed one (i.e. keep a copy of what was there before)
	virtual double TestNewBranchLenAt(int threadId, int branch, double lenNew, map<int,set< LineageConfig,LTLinCfg > > &origLinCfgs, bool fSetBrlen );
	virtual void SetLinCfgs(map<int,set< LineageConfig,LTLinCfg > > &linCfgsToSet);	// useful when we want to un-do the effect of the previous tweaking of branch length
    virtual double GetSTBrachLen(int br) const;
	virtual void SetSTBranchLen(int br, double brLen);
	virtual void GetSpeciesTree(MarginalTree &spTree) const { spTree = treeSpecies; }
	virtual int GetTaxaAtLeafST(int snid) { return pgstHelper->GetTaxaAtLeafST(snid); }
	virtual void GetGeneAllelesForSpecies( int taxon, set<int> &geneAlleles) { pgstHelper->GetGeneAllelesForSpecies( taxon, geneAlleles); }
    virtual std::string GetSpeciesTreeNW() const
    {
        return treeSpeciesUse.GetNewickSorted(true);
    }

private:
	void Reset() { pstoreLinCfgs->Reset(); pgstHelper->Reset(); }

	MarginalTree &treeSpecies;
    MarginalTree treeSpeciesUse;        // sometimes some taxa are missing from gene trees; if this occur, we only use the subtree with the subset of taxa in the gene tree
	PhylogenyTreeBasic &treeGene;
	GeneSpeciesTreeHelper *pgstHelper;
	LineageConfigSore *pstoreLinCfgs;
    map<int,int> mapNewNodeToOldNode;
};



#endif
