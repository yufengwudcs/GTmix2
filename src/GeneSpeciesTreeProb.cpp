#include "GeneSpeciesTreeProb.h"
#include <cmath>
#include <iomanip>
#include <stack>
#include "UtilsNumerical.h"
//#include "partition-utils.h"
//#include <multiset>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// common function

static int maxLinCfgNum = HAP_MAX_INT;
static int maxLinAncCfgNum = HAP_MAX_INT;
static bool fHeuModeSet = false;
static bool fVerboseMode = false;
bool fSTELLSProbComputeWarningFlag = false;
bool fFixedCoalMode = false;

bool GetSTELLSProbComputeWarningFlag()
{
    //
    return fSTELLSProbComputeWarningFlag;
}

void SetFixedCoalCalcMode(bool f)
{
    fFixedCoalMode = f;
}

bool IsFixedCoalMode()
{
    return fFixedCoalMode;
}

// get a distinct label
static string GetDistinctString(const string &strSeed)
{
    return strSeed;
#if 0
    static int indexString = 1;
    char buf[100000];
    sprintf(buf, "%d", indexString++);
    string res = strSeed + buf;
    return res;
#endif
}

// test brranch prob
void TestCalcProbBranch()
{
    //
    int numLins = 10;
    double listLens[4] = {0.001, 0.01, 0.1, 1.0};
    for(int i=0; i<4; ++i)
    {
        double len = listLens[i];
cout << "*************************** length: " << len << endl;
        for(int v = 1; v<= numLins; ++v)
        {
            double probAcc = UtilsCalcCoalProbBranch( numLins, v, len );
            double probApproxPoisson = UtilsCalcCoalProbBranchPoissonApprox( numLins, v, len );
            double probApproxNormal = UtilsCalcCoalProbBranchNormalApprox( numLins, v, len );
            cout << "v: " << v << ", probAcc = " << probAcc << ", probApproxPoisson = " << probApproxPoisson << ", probApproxNormal = " << probApproxNormal << endl;
        }
    }
}

#if 0

void EraseCommonItemsFrom( vector<int> &listItems1, vector<int> &listItems2)
{
	// remove shared common items
	// first sort the list
	SortIntVec(listItems1);
	SortIntVec(listItems2);
//cout << "Before EraseCommonItemsFrom: \n";
//DumpIntVec(listItems1);
//DumpIntVec(listItems2);
	vector<int> listItemNew1, listItemNew2;
	// iterate through the two list concurrently, and avoid one common item when needed
	int pos1 = 0, pos2 = 0;
	while( pos1 <(int)listItems1.size() && pos2 <(int)listItems2.size() )
	{
		// if one item is bigger than move it
		if( listItems1[pos1] < listItems2[pos2] )
		{
			// put the item to new list
			listItemNew1.push_back( listItems1[pos1] );
			pos1++;
		}
		else if( listItems1[pos1] > listItems2[pos2]  )
		{
			listItemNew2.push_back( listItems2[pos2] );
			pos2++;
		}
		else
		{
			// move together but skip the common items
			pos1++;
			pos2++;
		}
	}
	// now add whatever left over to the two list
	for(int i=pos1; i<(int)listItems1.size(); ++i)
	{
		listItemNew1.push_back( listItems1[i] );
	}
	for(int i=pos2; i<(int)listItems2.size(); ++i)
	{
		listItemNew2.push_back( listItems2[i] );
	}
	listItems1 = listItemNew1;
	listItems2 = listItemNew2;
//cout << "AFTER EraseCommonItemsFrom: \n";
//DumpIntVec(listItems1);
//DumpIntVec(listItems2);
}

#endif

void SetMaxLineageCfgKept(int maxn)
{
	maxLinCfgNum = maxn;
    // if this is ever called, then heuristic
    fHeuModeSet = true;
}

void SetMaxLineageAncesCfgKept(int maxn)
{
	// set maximum ancestral cfgs to keep (different from the above, this limits the number of found ancestral cfgs during ancestral cfgs construction)
	maxLinAncCfgNum = maxn;
    fHeuModeSet = true;
}

void SetVerbose(bool f)
{
	fVerboseMode = f;
}

bool IsVerbose()
{
    return fVerboseMode;
}

double CalcCombNum(int n, int k)
{
	// number of choosing k items among n (distinct) items
	double res = 1.0;
	// 
	for(int i=1; i<=k; ++i)
	{
		res *= ((double) n-k+i )/((double) i);
	}

	return res;
}

double UtilsCalcCoalProbBranch( int u, int v, double len )
{
    // YW: test approx
    const int THRES_APPROX_U_VAL = 30;
    const double THRES_APPROX_LEN = 0.0001;
    if( u > THRES_APPROX_U_VAL || len < THRES_APPROX_LEN )
    {
        // use approximation if approximation is needed
        return UtilsCalcCoalProbBranchPoissonApprox( u, v, len );
    }
    
//cout << "CalcCoalProbBranch: u=" << u << ", v=" << v << ", len=" << len;
	// prob of coalsecing from u lineage to v lineages with len time
	// CAUTION: we assume u >= v > =1
	YW_ASSERT_INFO( u >= v && v >= 1, "Here, it must be u>=v>=1" );
	long double res = 0.0;
	//double kfac=1.0;
	//for(int i=1; i<=v; ++i)
	//{
	//	kfac *= i;
	//}
	//long double resfac = exp(  (long double)(-0.5*(v-1)*v*len) );
	//vector<long double> listFacsStep;
	//vector<long double> listExpoStep;
	//vector<long double> listValsToAdd;

	//int expoOrderAccum = 0;
	for(int k=u; k>=v; --k)
	{
		vector<int> listNumItems;
		vector<int> listDenomItems;

		long double vexp = exp( (long double)(-0.5*((k-1)*k)*len) );
		//listExpoStep.push_back(vexp);
		//expoOrderAccum += (k-1)*k;
		//listNumItems.push_back( vexp );

//cout << "k=" << k << ", exp=" << vexp ;
		int fsign = 1;
		if( (k-v) % 2 == 1 )
		{
			fsign = -1;
		}
		//listNumItems.push_back( fsign);

		listNumItems.push_back( 2*k-1 );
		listDenomItems.push_back( v+k-1 );

		for(int y=0; y<=k-1; ++y)
		{
			listNumItems.push_back( v+y );
			listNumItems.push_back( u-y );
			listDenomItems.push_back(u+y);
		}
		for(int jj=1; jj<=k-v; ++jj)
		{
			listDenomItems.push_back( jj );
		}
		for(int jj=1; jj<=v; ++jj)
		{
			listDenomItems.push_back( jj );
		}
		EraseCommonItemsFrom(listNumItems, listDenomItems);
		//SortIntVec( listNumItems );
		//SortIntVec(listDenomItems);
		//ReverseIntVec(listNumItems);
		//ReverseIntVec(listDenomItems);
		long double valc = fsign*vexp;
		long double valc2 = 1.0;
		int ind = 0;
		for( ; ind <(int)listNumItems.size() && ind <(int)listDenomItems.size(); ++ind )
		{
			valc2 *= ((long double)listNumItems[ind]/listDenomItems[ind]);
		}
		for(int jj=ind; jj<(int)listNumItems.size(); ++jj)
		{
			valc2 *= listNumItems[jj];			
		}
		for(int jj=ind; jj<(int)listDenomItems.size(); ++jj)
		{
			valc2 = valc2 / listDenomItems[jj];			
		}
		//listFacsStep.push_back( fsign*valc2 );
		res += valc*valc2;
		//listValsToAdd.push_back(valc*valc2);
//if( u == 32 && v == 10)
//if( res < 0.0 && k==u)
//{
//cout << "k = " << k << ", u = " << u << ", v = " << v << ", len = " << len << ", valc = " << valc << ", res = " << res << endl;
//cout << "ListNumerators: ";
//DumpIntVec(listNumItems);
//cout << "ListDenomItems: ";
//DumpIntVec(listDenomItems);
//}
//cout << ", valc = " << valc << endl;
	}
#if 0
	int pos=0;
	for(; pos<(int)listValsToAdd.size()-1; pos = pos+2)
	{
		res += (listValsToAdd[pos] + listValsToAdd[pos+1]);
	}
	if( pos == (int)listValsToAdd.size()-1 )
	{
		res += listValsToAdd[pos];
	}
#endif
	// now combine to get the results using a convolutionary way. This helps to obtian more numerical stable computation
	//YW_ASSERT_INFO( listFacsStep.size() == listExpoStep.size(), "Size mismatch" );
	//res = 0.0;
	//for( int i=(int)listFacsStep.size()-1; i>=0; --i )
	//{
	//	res = (res + listFacsStep[i])*listExpoStep[i];
	//}


	// in case numerical underflow, let it be 0.0
	if( res < 0.0)
	{
 		if( fVerboseMode == true)
		{
			cout << "*Caution* A minor numerical issue detected: change value " << res << " to 0.0, u=" << u << ", v=" << v << ", len=" << len << "\n";
		}
		res = 0.0;
	}

    if( res > 1.0)
    {
        res = 1.0;
        //cout << "u: " << u << ",v: " << v << ", len: " << len << ", prob=" << res << endl;
    }
    

	//YW_ASSERT_INFO( res <= 1.0, "UtilsCalcCoalProbBranch: prob can not be larger than 1.0" );
#if 0
set<pair<int,int> > setDonePairs;
pair<int,int> pp(u,v);
if( setDonePairs.find( pp ) == setDonePairs.end() )
{
cout << "CalcCoalProbBranch: u=" << u << ", v=" << v << ", len=" << len;
cout << ": prob = " << res << endl;
setDonePairs.insert(pp);
}
#endif
	return (double)res;
}


double UtilsCalcCoalProbBranch2( int u, int v, double len )
{
    // YW: test approx
    const int THRES_APPROX_U_VAL = 50;
    //const int THRES_APPROX_U_VAL = 30;
    const double THRES_APPROX_LEN = 0.0001;
    if( u > THRES_APPROX_U_VAL || len < THRES_APPROX_LEN )
    {
        // use approximation if approximation is needed
        return UtilsCalcCoalProbBranchPoissonApprox2( u, v, len );
    }
    
//cout << "CalcCoalProbBranch: u=" << u << ", v=" << v << ", len=" << len;
    // prob of coalsecing from u lineage to v lineages with len time
    // CAUTION: we assume u >= v > =1
    YW_ASSERT_INFO( u >= v && v >= 1, "Here, it must be u>=v>=1" );
    long double res = 0.0;
    //double kfac=1.0;
    //for(int i=1; i<=v; ++i)
    //{
    //    kfac *= i;
    //}
    //long double resfac = exp(  (long double)(-0.5*(v-1)*v*len) );
    //vector<long double> listFacsStep;
    //vector<long double> listExpoStep;
    //vector<long double> listValsToAdd;

    //int expoOrderAccum = 0;
    for(int k=u; k>=v; --k)
    {
        vector<int> listNumItems;
        vector<int> listDenomItems;

        long double vexp = exp( (long double)(-0.5*((k-1)*k)*len) );
        //listExpoStep.push_back(vexp);
        //expoOrderAccum += (k-1)*k;
        //listNumItems.push_back( vexp );

//cout << "k=" << k << ", exp=" << vexp ;
        int fsign = 1;
        if( (k-v) % 2 == 1 )
        {
            fsign = -1;
        }
        //listNumItems.push_back( fsign);

        listNumItems.push_back( 2*k-1 );
        listDenomItems.push_back( v+k-1 );

        for(int y=0; y<=k-1; ++y)
        {
            listNumItems.push_back( v+y );
            listNumItems.push_back( u-y );
            listDenomItems.push_back(u+y);
        }
        for(int jj=1; jj<=k-v; ++jj)
        {
            listDenomItems.push_back( jj );
        }
        for(int jj=1; jj<=v; ++jj)
        {
            listDenomItems.push_back( jj );
        }
        EraseCommonItemsFrom(listNumItems, listDenomItems);
        //SortIntVec( listNumItems );
        //SortIntVec(listDenomItems);
        //ReverseIntVec(listNumItems);
        //ReverseIntVec(listDenomItems);
        long double valc = fsign*vexp;
        long double valc2 = 1.0;
        int ind = 0;
        for( ; ind <(int)listNumItems.size() && ind <(int)listDenomItems.size(); ++ind )
        {
            valc2 *= ((long double)listNumItems[ind]/listDenomItems[ind]);
        }
        for(int jj=ind; jj<(int)listNumItems.size(); ++jj)
        {
            valc2 *= listNumItems[jj];
        }
        for(int jj=ind; jj<(int)listDenomItems.size(); ++jj)
        {
            valc2 = valc2 / listDenomItems[jj];
        }
        //listFacsStep.push_back( fsign*valc2 );
        res += valc*valc2;
        //listValsToAdd.push_back(valc*valc2);
//if( u == 32 && v == 10)
//if( res < 0.0 && k==u)
//{
//cout << "k = " << k << ", u = " << u << ", v = " << v << ", len = " << len << ", valc = " << valc << ", res = " << res << endl;
//cout << "ListNumerators: ";
//DumpIntVec(listNumItems);
//cout << "ListDenomItems: ";
//DumpIntVec(listDenomItems);
//}
//cout << ", valc = " << valc << endl;
    }
#if 0
    int pos=0;
    for(; pos<(int)listValsToAdd.size()-1; pos = pos+2)
    {
        res += (listValsToAdd[pos] + listValsToAdd[pos+1]);
    }
    if( pos == (int)listValsToAdd.size()-1 )
    {
        res += listValsToAdd[pos];
    }
#endif
    // now combine to get the results using a convolutionary way. This helps to obtian more numerical stable computation
    //YW_ASSERT_INFO( listFacsStep.size() == listExpoStep.size(), "Size mismatch" );
    //res = 0.0;
    //for( int i=(int)listFacsStep.size()-1; i>=0; --i )
    //{
    //    res = (res + listFacsStep[i])*listExpoStep[i];
    //}


    // in case numerical underflow, let it be 0.0
    if( res < 0.0)
    {
         if( fVerboseMode == true)
        {
            cout << "*Caution* A minor numerical issue detected: change value " << res << " to 0.0, u=" << u << ", v=" << v << ", len=" << len << "\n";
        }
        res = 0.0;
    }

    if( res > 1.0)
    {
        res = 1.0;
        //cout << "u: " << u << ",v: " << v << ", len: " << len << ", prob=" << res << endl;
    }
    

    //YW_ASSERT_INFO( res <= 1.0, "UtilsCalcCoalProbBranch: prob can not be larger than 1.0" );
#if 0
set<pair<int,int> > setDonePairs;
pair<int,int> pp(u,v);
if( setDonePairs.find( pp ) == setDonePairs.end() )
{
cout << "CalcCoalProbBranch: u=" << u << ", v=" << v << ", len=" << len;
cout << ": prob = " << res << endl;
setDonePairs.insert(pp);
}
#endif
    
//cout << "Exact coalescent prob: " << res << " for u:" << u << ", v:" << v << ", len:" << len << endl;
//cout << "For comparison: poisson approximation: " << UtilsCalcCoalProbBranchPoissonApprox2( u, v, len ) << endl;
//cout << "For comparison: normal approximation: " << UtilsCalcCoalProbBranchNormalApprox(u, v, len) << endl;

    return (double)res;
}


// *************************************************************************************
// ways to approximate
double UtilsCalcCoalProbBranchPoissonApprox( int u, int v, double len )
{
    // perform the simple Poisson approx as in Griffiths, 1984
    double mv = 0.5*u*(u-1)*len;
    int diff = u-v;
    YW_ASSERT_INFO(diff >=0, "Cannot be negative");
    double res = pow( mv, diff ) * exp( -1.0*mv );
    for(int i=1; i<=diff; ++i)
    {
        res = res/i;
    }
    return res;
}
double UtilsCalcCoalProbBranchPoissonApprox2( int u, int v, double len )
{
    // address underflow issue
    // perform the simple Poisson approx as in Griffiths, 1984
    double mv = 0.5*u*(u-1)*len;
    int diff = u-v;
    YW_ASSERT_INFO(diff >=0, "Cannot be negative");
    
    // determine if it is too large
    double res = 1.0;
    double ftest = diff*log(mv);
    const double MAXLOG = log(1.0e300);
    if( ftest >= MAXLOG )
    {
        // too large! set to a large number
        res = 1.0e300;
    }
    else
    {
        res = pow( mv, diff );
    }
//cout << "UtilsCalcCoalProbBranchPoissonApprox2: res1:" << res << endl;
    res *= exp( -1.0*mv );
//cout << "res2:" << res << endl;
    for(int i=1; i<=diff; ++i)
    {
        res = res/i;
    }
//cout << "UtilsCalcCoalProbBranchPoissonApprox2: res:" << res << " for u:" << u << ", v:" << v << ", len:" << len << endl;
//cout << "For comparison: normal approximation: " << UtilsCalcCoalProbBranchNormalApprox(u, v, len) << endl;
    return res;
}

double UtilsCalcCoalProbBranchNormalApprox( int u, int v, double len )
{
    // from Griffiths, 1984
    double a = 0.5*u*len;
    double b = -0.5*len;
    double c = a*b/((a+b)*exp(b)-a);
    double valMean = 2.0*c/len;
    double valVar = valMean*(c+b)*(c+b)*(1.0+  c/(c+b) - c/a -c/(a+b) -2.0*c  )/(b*b);
    double valSD = sqrt(valVar);
    
//cout << "UtilsCalcCoalProbBranchNormalApprox: a=" << a << ", b=" << b << ", c=" << c << ", valMean=" << valMean << ", valVar=" << valVar << endl;
    
    // now convert to standard nomral
    double v1 = ( v - valMean )/valSD;
    
    // use CDF of normal approximation
    double v2 = (v-1.0-valMean)/valSD;
    
    return CalcApproxCDFStdNormal(v1) - CalcApproxCDFStdNormal(v2);
}

// *************************************************************************************

int GetNumSpeciesFromGT( PhylogenyTreeBasic &treeGene)
{
	set<int> species;
	GetLeavesSpeciesFromGT(treeGene, species);
	return species.size();
}

int GetNumSpeciesFromGTs( const vector<PhylogenyTreeBasic *> &listGTreePtrs, set<int> *ptrSpecies )
{
    //
    set<int> speciesAll;
    for(int i=0; i<(int)listGTreePtrs.size(); ++i)
    {
        set<int> species;
        GetLeavesSpeciesFromGT(*listGTreePtrs[i], species);
        UnionSets( speciesAll, species );
    }
    if( ptrSpecies != NULL )
    {
        *ptrSpecies = speciesAll;
    }
    return speciesAll.size();
}


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

bool AreGeneTreesWithSameTaxa( const vector<PhylogenyTreeBasic *> &listGTP )
{
    // check if all gene trees has the same set of taxa
    bool res = true;
    set<int> species;
	for(int i=0;i<(int)listGTP.size(); ++i)
	{
		set<int> speciesStep;
		GetLeavesSpeciesFromGT( *listGTP[i], speciesStep );
		if( species.size() == 0 )
		{
			species = speciesStep;
		}
		if( species != speciesStep )
        {
            res = false;
            break;
        }
	}
    return res;
}

// useful functions
void GetSpeciesAtNode( TreeNode *pnode, set<int> &speciesUnder )
{
	YW_ASSERT_INFO( pnode != NULL, "Fail" );
	vector<string> listLeafLabels;
	pnode->GetAllLeafLabeles( listLeafLabels );
	speciesUnder.clear();
	for( int jj=0; jj<(int)listLeafLabels.size(); ++jj )
	{
		// what species is this label?
		int sp = GetSpeciesFromLabel( listLeafLabels[jj] );
		speciesUnder.insert( sp );
	}
}

// are we in heuristic mode?
//static bool GSTPIsHeuModeOn()
//{
//    return fHeuModeSet;
//}

///////////////////////////////////////////////////////////////////////////////////////////////////


LineageCluster :: LineageCluster( int gnid, int gnparid )
{
	Init( gnid, 1, gnparid );
}

void LineageCluster :: Append(const LineageCluster &lc)
{
	// make sure they are onto the same par
	YW_ASSERT_INFO( this->gnidPar == lc.gnidPar, "Can not append" );
	// add up num lineages
	this->numLins += lc.numLins;
}
void LineageCluster :: Dump() const
{
//	cout << "LC: gnid = " << gnid << ", numLins = " << numLins << ", par = " << gnidPar << endl;
	cout << gnid << " (" << gnidPar  << ")";
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// this initialize the config with a single lineage of GT
LineageConfig :: LineageConfig(  )
{
	this->snid = -1;
	this->valProb = 0.0;
	RsetLinNum();
}

//LineageConfig :: LineageConfig( int id, int gnid, int gnparid )
//{
//	this->snid = id;
//	LineageCluster lc;
//	lc.Init( gnid, 1, gnparid );
//	setLCs.insert(lc);
//}

bool LineageConfig :: operator <(const LineageConfig &lc2) const
{
	if( snid < lc2.snid )   //|| numTotLins < lc2.numTotLins )
	{
		return true;
	}
	if( snid == lc2.snid && setLCs.size() < lc2.setLCs.size() )
	{
		return true;
	}
	if( snid == lc2.snid && setLCs.size() == lc2.setLCs.size() /*&& setLCs < lc2.setLCs*/ )
	{

		for(  set< LineageCluster > :: iterator itt1  = setLCs.begin(), itt2 = lc2.setLCs.begin(); 
			itt1 !=setLCs.end() && itt2 != lc2.setLCs.end(); ++itt1, ++itt2 )
		{
			if( *itt1 < *itt2 )
			{
				return true;
			}
			else if( *itt2 < *itt1  )
			{
				return false;
			}
		}

		return false;
/*
		if( setLCs < lc2.setLCs )
		{
			return true;
		}
*/
	}
	// now compare each of the components
	return   false;
}

void LineageConfig :: Consolidate( const GeneSpeciesTreeHelper &gstHelper )
{
//cout << "LineageConfig :: Consolidate: num lineage clusters = " << setLCs.size() << endl;
	// consolidate the lineages into clusters if we can
	// Idea: just compare the parent node to see if they match. If they do, combine it
	map<int, set<LineageCluster> > mapLCs;
	for( set<LineageCluster> :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
		int pid = it->GetGNParId();
		if( mapLCs.find( pid ) == mapLCs.end() )
		{
			set<LineageCluster> setLCsEmpty;
			mapLCs.insert( map<int,set< LineageCluster> > :: value_type( pid, setLCsEmpty ) );
		}
		mapLCs[pid].insert( *it );
	}
	// do we have duplicates?
	set< LineageCluster > setLCsNew;
	for( map<int, set<LineageCluster> > :: iterator it = mapLCs.begin(); it != mapLCs.end(); ++it )
	{
		if( it->second.size() == 1 )
		{
			// no consolidate
			setLCsNew.insert( *(it->second).begin() );
		}
		else
		{
			// consolidate
			LineageCluster lcNew = *(it->second.begin() );
			set< LineageCluster > :: iterator it2 = it->second.begin();
			it2++;
			//set<int> nodeIds;
			int mrcaOverall = -1;
			int nodeLast = lcNew.GetGNId();
			for( ; it2 != it->second.end(); ++it2 )
			{
				lcNew.Append( *it2 );
				mrcaOverall = gstHelper.GetMRCATwoLineagesGT(it2->GetGNId(), nodeLast );
				nodeLast =  it2->GetGNId();
			}
			// update the node id for this new lineage
			YW_ASSERT_INFO( mrcaOverall >= 0, "MRCA: not found" );
			lcNew.SetGNId( mrcaOverall );
			// if node id changed, then also change its parent
			lcNew.SetGNParId( gstHelper.GetParentLineage(mrcaOverall) );

			// need to assign the MRCA of the two lineages
			setLCsNew.insert( lcNew );
		}
	}
	// setup new
	this->setLCs = setLCsNew;
	//RsetLinNum();
//cout << "LineageConfig :: Consolidate: num lineage clusters = " << setLCs.size() << endl;
}

void LineageConfig :: FindCoalescableLins( set< set<LineageCluster> > &setCoalLins ) const
{
//cout << "FindCoalescableLins: begin\n";
	setCoalLins.clear();
	//
	// Idea: just compare the parent node to see if they match. If they do, combine it
	map<int, set<LineageCluster> > mapLCs;
	for( set<LineageCluster> :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
//it->Dump();
		int pid = it->GetGNParId();
		if( mapLCs.find( pid ) == mapLCs.end() )
		{
			set<LineageCluster> setLCsEmpty;
			mapLCs.insert( map<int,set< LineageCluster> > :: value_type( pid, setLCsEmpty ) );
		}
		mapLCs[pid].insert( *it );
	}
	// do we have duplicates?
	for( map<int, set<LineageCluster> > :: iterator it = mapLCs.begin(); it != mapLCs.end(); ++it )
	{
		setCoalLins.insert(it->second);
	}
//cout << "FindCoalescableLins: end\n";
}

void LineageConfig :: Dump() const
{
	cout << "Species tree node: " << snid << ", num of lineages = " << setLCs.size() << ", prob = " << this->valProb << endl;
	int numDump=0;
	bool fRet = false;
	for( set< LineageCluster > :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
		it->Dump();
		cout << "      ";
		if( ++numDump == 8)
		{
			numDump = 0;
			cout << endl;
			fRet = true;
		}
		else
		{
			fRet = false;
		}
	}
	if( fRet == false)
	{
		cout << endl;
	}
}


void LineageConfig :: CopyReduceOne( const LineageCluster &lcToReduce, LineageConfig &lcNew, int gtnodeNew, int gtnodeParNew ) const
{
//cout << "CopyReduceOne: gtnodeNew = " <<  gtnodeNew << ", gtnodeParNew = " << gtnodeParNew << ". lcToReduce = ";
//lcToReduce.Dump();
	//lcNew.setLCs =  this->setLCs;
	lcNew.snid = this->snid;
	lcNew.valProb = this->valProb;
	// now search for the item to reduce
	for( set< LineageCluster > :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
//cout << "current linclust: ";
//it->Dump();
		if( *it == lcToReduce )
		{
//cout << "Now reduce it\n";
			//
			LineageCluster linClusRed = *it;
			linClusRed.DecByOne();
			//lcNew.setLCs.erase( *it );
			//break;
			// do we need to setup new node info
			if( gtnodeNew >= 0)
			{
				linClusRed.Init( gtnodeNew, linClusRed.GetNumLineages(), gtnodeParNew );
			}
			lcNew.setLCs.insert( linClusRed );
		}
		else
		{
			lcNew.setLCs.insert(*it);
		}
	}
//cout << "Current lcNew: ";
//lcNew.Dump();
	// add back by
	//LineageCluster lcNew2 = lcToReduce;
	//lcNew2.DecByOne();

	// do we need to setup new node info
	//if( gtnodeNew >= 0)
	//{
	//	lcNew2.Init( gtnodeNew, lcNew2.GetNumLineages(), gtnodeParNew );
	//}

	//lcNew.setLCs.insert( lcNew2 );
}

int LineageConfig :: GetTotNumLins() const
{
	// we only count the number of lin cluster (but we do not add up the size of cluster)
	// this is because the cluster size meaning has been changed: it means there 
	// are some lineages should be merged
	return setLCs.size();
}

void LineageConfig :: GetSrcNodes(set<int> &nodesIds) const
{
	// 
	nodesIds.clear();
	for( set< LineageCluster > :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
		int id = it->GetGNId();
		YW_ASSERT_INFO( nodesIds.find( id ) == nodesIds.end(), "Id should not be already in set" );
		nodesIds.insert( id );
	}
}

void LineageConfig :: GetSrcNodesNumLins( map<int,int> &mapLinInfo) const
{
	// this one return how many COALSECENT events allowed at each node
	// thus, at a leaf, it is 0
	// at an internal node with k downwards lineage, it is k-1

	mapLinInfo.clear();
	for( set< LineageCluster > :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
		int id = it->GetGNId();
		int nl = it->GetNumLineages();
		// YW: changed at 4/11 because num lineage of cluster has been changed to mean collapsed lineages at that lins
		//int nl = 1;	
		//YW_ASSERT_INFO(nl != 0, "Do not allow intermediate node with a single downward edge");
		if( nl < 0 )
		{
			nl = 0;
		}
		//YW_ASSERT_INFO( nodesId.find( id ) == nodesId.end(), "Id should not be already in set" );
		mapLinInfo.insert( map<int,int> :: value_type( id, nl )  );
	}
}

bool LineageConfig :: IsAncestralTo( const LineageConfig &cfgOrig, const GeneSpeciesTreeHelper &gstHelper ) const
{
	// can cfgOrig coalesce into the current one
	// the rule is, if so, then each lineage in this must be ancestral to some src (orig) lineages
	for( set< LineageCluster > :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
		bool fAnces = false;
		int linAnc = it->GetGNId();
		for( set< LineageCluster > :: iterator it2 = cfgOrig.setLCs.begin(); it2 != cfgOrig.setLCs.end(); ++it2 )
		{
			// 
			int linOrig = it2->GetGNId();
			if( linAnc == linOrig || gstHelper.IsNodeAncestralToGT(linAnc, linOrig) == true)
			{
				fAnces = true;
				break;
			}
		}
		if( fAnces == false)
		{
			return false;
		}
	}
	return true;
}

void LineageConfig :: MergeLCInto(const GeneSpeciesTreeHelper &gstHelper, 
								  const LineageCluster &lin1, const LineageCluster &lin2, LineageConfig &lcNew ) const
{
//YW_ASSERT_INFO(false, "Should not be here");
	// NOTE: this does not quite work with non-binary gene tree!!!!!!!!!!!!!!!!!!
	// Potential future BUG. YW: TBD
	// 
	// consolidate
	//lcNew.snid = this->snid;
	//lcNew.valProb = this->valProb;
	//lcNew.numTotLins = this->numTotLins;
	//lcNew.setLCs = this->setLCs;

	lcNew = *this;
	
	YW_ASSERT_INFO( lcNew.setLCs.find( lin1 ) != lcNew.setLCs.end(), "Wrong1" );
	YW_ASSERT_INFO( lcNew.setLCs.find( lin2 ) != lcNew.setLCs.end(), "Wrong2" );
	lcNew.setLCs.erase(lin1);
	lcNew.setLCs.erase(lin2);
	// this is how many accumulated coalesced lineages
	// when this num reaches the num of total descent num, we will then move to its parent lin (meaning we are done with
	// all the needed lineages at this node)
	int numLins = lin1.GetNumLineages() + lin2.GetNumLineages();
	int nodeParLin1 = lin1.GetGNParId();
	int numDescPar = gstHelper.GetGTNodeNumChildren(nodeParLin1);
	YW_ASSERT_INFO( numDescPar >= 0, "numDescPar: not found" );
	LineageCluster linNew;
	if( numLins < numDescPar )
	{
		// continue with the first lin but with combined lin number
		linNew.Init(lin1.GetGNId(), numLins, lin1.GetGNParId() );
	}
	else
	{
		// in this case we always reduce to 1 lineage (since all lineages are coalesced
		//int mrca = gstHelper.GetMRCATwoLineagesGT(  lin1.GetGNId(), lin2.GetGNId() );
		linNew.Init(nodeParLin1, 1, gstHelper.GetParentLineage(nodeParLin1) );
	}
	lcNew.AddLC( linNew );


	// append them. NOTE: YW: 040110: do not work with non-binary tree YET
	// TBD
	//LineageCluster linNew;
	//linNew.Init(mrca, numLins, gstHelper.GetParentLineage(mrca) );
	//lcNew.AddLC( linNew );
}

void LineageConfig :: RemoveLineage(const GeneSpeciesTreeHelper &gstHelper, 
								  const LineageCluster &lin1, LineageConfig &lcNew ) const
{
	// remove one lineage from the current config
	//lcNew.snid = this->snid;
	//lcNew.valProb = this->valProb;
	//lcNew.numTotLins = this->numTotLins;
	//lcNew.setLCs = this->setLCs;

	lcNew = *this;
	
	YW_ASSERT_INFO( lcNew.setLCs.find( lin1 ) != lcNew.setLCs.end(), "Wrong1" );
	lcNew.setLCs.erase(lin1);
	const_cast<LineageConfig *>(this)->numTotLins = -1;
	YW_ASSERT_INFO( lcNew.setLCs.size() > 0, "Can not reduce to empty" );
	// append them. NOTE: YW: 040110: do not work with non-binary tree YET
	// This should work with non-binary tree. Need test-out
}

void LineageConfig :: GetParLins(multiset<int> &parids) const
{
	parids.clear();
	for( set< LineageCluster > :: iterator it = setLCs.begin(); it != setLCs.end(); ++it )
	{
		parids.insert( it->GetGNParId() );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////

void AncLinConfig :: UpdateProb(double len, GeneSpeciesTreeHelper *gstHelper)
{
	YW_ASSERT_INFO(pSrcConfig != NULL, "src cfg can not be NULL");
	
	// update its prob from ancestral prob by multiplying the factors due to num of coal and evts
	int numEvts = listEvts.size();
	int ub = pSrcConfig->GetTotNumLins();
	int cb = numEvts;
if( ub-cb < 1 )
{
cout << "ub = " << ub << ", cb = " << cb << endl;
Dump();
}
	YW_ASSERT_INFO( ub-cb>=1, "Can not vanish" );
	//double coeff=1.0/gstHelper->CalcdbVal(ub, cb);
	double coeff = 1.0;
    //if( GSTPIsHeuModeOn() == true )
    //{
    //    coeff *= gstHelper->CalcwbValMultiHeu( listEvts, *pSrcConfig, ub, cb );
    //}
    //else
    //{
    if( gstHelper->IsGeneTreeMulfurcate() == false)
    {
        coeff *= gstHelper->CalcwbVal(listEvts, *pSrcConfig, ub, cb);
    }
    else
    {
        coeff=1.0/gstHelper->CalcdbVal(ub, cb);
        coeff *= gstHelper->CalcwbValMulti( *pSrcConfig, cfgCurr );
    }
	//}
//	double probCoeff = gstHelper.CalcCoalCoeffForTwoCfg(lcRightSrc, lcRightAnc);
	double probBranch = gstHelper->CalcCoalProbBranch(ub, ub-cb, len);
	// update prob
	double srcProb = pSrcConfig->GetProb();
	SetProb( srcProb*probBranch*coeff );
	if( srcProb < 0.0 )
	{
		cout << "AncLinConfig :: UpdateProb: src prob is negative\n";
	}
	if( probBranch < 0.0 )
	{
		cout << "AncLinConfig :: UpdateProb: probBranch is negative\n";
		cout << "ub: " << ub << ", cb = " << cb << ", len = " << len << endl;
	}
	if( coeff < 0.0 )
	{
		cout << "AncLinConfig :: UpdateProb: coeff is negative\n";
	}
    if( fVerboseMode == true )
    {
        cout << "UpdateProb: numEvts=" << numEvts << ", ub=" << ub << ",cb=" << cb << ",coeff=" << coeff << ", probBranch=" << probBranch << ",srcProb=" << srcProb << endl;
        Dump();
    }
}

void AncLinConfig :: AppendAC(const AncLinConfig &cfgAdd)
{
	// add the prob to the current prob
	SetProb( this->GetProb() + cfgAdd.GetProb() );
}

void AncLinConfig :: Dump() const
{
	cout << "List of events: ";
	for(int i=0; i<(int)listEvts.size(); ++i)
	{
		listEvts[i].Dump();
	}
	cout << endl;
	cout << "Lin cfgs: ";
	cfgCurr.Dump();
	YW_ASSERT_INFO( pSrcConfig != NULL, "Can not be NULL" );
	cout << "Source cfg: ";
	pSrcConfig->Dump();
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

const int MAX_NUM_LINS_PER_SPECIES = 15;

MultiLineageHelper  MultiLineageHelper :: instance;


MultiLineageHelper :: MultiLineageHelper( )
{
	// pre-compute certain number of table entries for faster lookup
	for(int nsrc =2; nsrc <= MAX_NUM_LINS_PER_SPECIES; ++nsrc)
	{
		double res = 1.0;
		for(int ndest = nsrc; ndest >=1; --ndest)
		{
			int nevt = nsrc-ndest;
			if( nevt >=1)
			{
				res *= CalcCombNum(ndest+1, 2)/nevt;
			}
			pair<int,int> pp(nsrc,ndest);
			mapPreCMultiLinsCoeff.insert( map<pair<int,int>, double> :: value_type(pp, res));
//cout << "nsrc = " << nsrc << ", ndest = " << ndest << ", mapPreCMultiLinsCoeff[pp] = " << res << endl;;
		}
	}
}

double MultiLineageHelper :: GetMultiLinsCoeff(int numLinsSrc, int numLinsDest)
{
	YW_ASSERT_INFO(numLinsSrc >= numLinsDest, "Wrong");
	if( numLinsSrc == 1 )
	{
		return 1.0;
	}
	if( numLinsSrc > MAX_NUM_LINS_PER_SPECIES )
	{
		// for now, abort
		YW_ASSERT_INFO(false, "Too large mulfurcation: Not supported yet");
		return 0.0;
	}
	pair<int,int> pp(numLinsSrc, numLinsDest);
	YW_ASSERT_INFO(mapPreCMultiLinsCoeff.find(pp) != mapPreCMultiLinsCoeff.end(), "Must find it");
	return mapPreCMultiLinsCoeff[pp];
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Helper class for determinting constraints imposed by
// species tree about which transitions are allowed

// useful function
int GetSpeciesFromLabel(const string &label) 
{
	int id = -1;
	sscanf( label.c_str(), "%d", &id );
	return id;
}



GeneSpeciesTreeHelper :: GeneSpeciesTreeHelper(MarginalTree &ts, PhylogenyTreeBasic &tg) :
treeSpecies( ts), treeGene( tg ), fMulfurcateTree(false)
//, szMultiNoDescAllowed(MAX_NUM_LINS_PER_SPECIES)
{
	fMulfurcateTree = treeGene.IsMulfurcate();
	Init();
}

void GeneSpeciesTreeHelper :: Init()
{
	// ???
	//YW_ASSERT_INFO(false, "");

	// collect gene tree parent
	vector<int> nidsList, nparsList;
	treeGene.GetNodeParInfoNew(nidsList, nparsList);
	YW_ASSERT_INFO(nidsList.size() == nparsList.size(), "Fail");
	// setup info: we will simply use node position as ID
	for(int i=0; i<(int)nparsList.size(); ++i)
	{
//cout << "Node " << nidsList[i] << ", has parent: " << nparsList[i] << endl;
		mapNodeParGT.insert( map<int,int> :: value_type(nidsList[i], nparsList[i] ) );

		// also maintain num of descendents
		if( mapNodesDescendentsGT.find( nparsList[i] ) == mapNodesDescendentsGT.end() )
		{
			set<int> sempty;
			mapNodesDescendentsGT.insert(map<int,set<int> > :: value_type(nparsList[i], sempty) );
		}
		mapNodesDescendentsGT[ nparsList[i] ].insert( nidsList[i]  );
//cout << "In GT, node " << nparsList[i] << " has descendents so far: ";
//DumpIntSet( mapNodesDescendentsGT[ nparsList[i]  ] );
	}

// ensure the num of descent is correct

	//treeGene.InitPostorderWalk();
	//while(true)
	//{
	//	// get the node
	//	TreeNode *pnode = treeGene.InitPostorderWalk();
	//}

	// pre-compute MRCA info
	int totLinNumGT = treeGene.GetNumVertices();
	vecMRCATwoLinsGT.resize( totLinNumGT );		// for two lineages a and b, their MRCA
	vecNodeAncestralGT.resize( totLinNumGT );	// query like is a is b's ancester?
	for(int i=0; i<totLinNumGT; ++i)
	{
		vecMRCATwoLinsGT[i].resize(totLinNumGT);
		vecNodeAncestralGT[i].resize(totLinNumGT);
		for(int jj=0; jj<totLinNumGT; ++jj)
		{
			vecMRCATwoLinsGT[i][jj] = -1;
			vecNodeAncestralGT[i][jj] = false;
		}
	}
	// now traverse the tree and find the MRCA
	// simple approach: traverse the tree; for each node,
	// find its left and right and then set their MRCA to the current one
	PhylogenyTreeIterator itorTree(treeGene);
	itorTree.Init();
	//treeGene.InitPostorderWalk();
	while(itorTree.IsDone() == false )
	{
//		TreeNode *pn = treeGene.NextPostorderWalk( ) ;
		TreeNode *pn = itorTree.GetCurrNode( ) ;
		itorTree.Next();
		if( pn == NULL )
		{
			break;      // done with all nodes
		}
		int idSelf = pn->GetID();
		YW_ASSERT_INFO( idSelf <(int)vecMRCATwoLinsGT.size(), "ID out of range error" );

		// keep track of id-->node info
//cout << "In Gene tree traversal, find a node idSelf: " << idSelf << endl;
		mapIdNode.insert( map<int, TreeNode*> :: value_type( idSelf, pn )  );

		if( pn->IsLeaf() == true)
		{
			continue;
		}


		// collect children for all branches underneath
		vector< set<int> > listDescNodeIds;
		for(int i=0; i<pn->GetChildrenNum(); ++i)
		{
			set< TreeNode *> setDescNodes;
			pn->GetChild(i)->GetAllDescendents(setDescNodes);
			set<int> idList;
			for(set<TreeNode *> :: iterator itt=setDescNodes.begin(); itt != setDescNodes.end(); ++itt)
			{
				int idToAdd = (*itt)->GetID();
				YW_ASSERT_INFO( idToAdd <(int)vecMRCATwoLinsGT.size(), "ID out of range error" );
				idList.insert( idToAdd );
			}
			listDescNodeIds.push_back( idList );
		}
		// now setup everything
		for(int i=0; i<pn->GetChildrenNum(); ++i)
		{
			// first, i to self is self
			for(set<int> :: iterator itki = listDescNodeIds[i].begin(); itki != listDescNodeIds[i].end(); ++itki)
			{
				vecMRCATwoLinsGT[ *itki ][ idSelf ] = idSelf;
				vecMRCATwoLinsGT[idSelf][ *itki ] = idSelf;
				YW_ASSERT_INFO( *itki != idSelf, "Can not the same" );
				vecNodeAncestralGT[idSelf][*itki] = true;
//cout << "***Setting " << idSelf << " to ancestral to " << *itki << endl;
				for(int j=i+1; j<pn->GetChildrenNum(); ++j)
				{
					for(set<int> :: iterator itkj = listDescNodeIds[j].begin(); itkj != listDescNodeIds[j].end(); ++itkj)
					{
						vecMRCATwoLinsGT[ *itki ][ *itkj ] = idSelf;
						vecMRCATwoLinsGT[*itkj][ *itki ] = idSelf;
					}
				}
			}

		}


	}

	// now initialize the mapping of species id from gene tree id
	InitDefaultGTLinMap();

	// initialize the coefficient map
	//InitMultiNoDescCoeff();

}

#if 0
void GeneSpeciesTreeHelper :: InitMultiNoDescCoeff()
{
	YW_ASSERT_INFO( szMultiNoDescAllowed >0, "szMultiNoDescAllowed must be at least one" );
	vecMultiNoDescCoeff.clear();
	vecMultiNoDescCoeff.resize(szMultiNoDescAllowed+1);
	// this is a DP approach. CAUTION: for ease of access, I use T[i] = i events are at the same
	// node (that is i+1 lineages coalescing)
	vecMultiNoDescCoeff[0] = 1.0;
	vecMultiNoDescCoeff[1] = 1.0;
	for(int numEvts =2; numEvts <= szMultiNoDescAllowed; ++numEvts)
	{
		int numLins = numEvts + 1;	// this is how many lineages we have to consider
		int limit = numLins/2;		// caution: over-counting for the middle item if numLins is even 
// cout << "numEvts: " << numEvts  << ", limit = " << limit << endl;

		// we need to consider the number of ways of spliting the nodes into (k,n-k) partitions
		double resStep = 0.0;	// Note the topmost is always 1.0 divided by numEvts, which will be considered later
		for(int i=1; i<= limit; ++i)
		{
			// compute the stepwise
			double numPartWays = CalcCombNum(numLins, i);
//cout << "numPartWays = " << numPartWays << " for i=" << i << endl;
			int evtLeft = i-1;
			int evtRight = numLins - i - 1;
			YW_ASSERT_INFO( evtRight >= 0, "Can not be over" );
			double resStepi = numPartWays*vecMultiNoDescCoeff[evtLeft]*vecMultiNoDescCoeff[evtRight];
			if( (numLins % 2) == 0 && i==limit )
			{
				resStepi *= 0.5;
			}
			resStep += resStepi;
//cout << "resStepi = " << resStepi << endl;
		}
		resStep *= 1.0/numEvts;
		vecMultiNoDescCoeff[numEvts] = resStep;
cout << "numEvts = " << numEvts << ", vecMultiNoDescCoeff[i] = " << resStep << endl;;
	}
}
#endif

double GeneSpeciesTreeHelper :: GetMultiNoDescCoeff(int numLinsOrig, int numDupEvts)
{
	YW_ASSERT_INFO( numDupEvts <= numLinsOrig-1, "Not enough to coalesce" );
	// 
	//YW_ASSERT_INFO( numEvts <= szMultiNoDescAllowed, "Num of intra-species lineages too large. Need to enlarge the settings" );
	int numDestLins = numLinsOrig - numDupEvts;
	//return vecMultiNoDescCoeff[numEvts];
	return MultiLineageHelper :: Instance().GetMultiLinsCoeff(numLinsOrig, numDestLins);
}


// support query
void GeneSpeciesTreeHelper :: GetAllSubtreeST( set<int> &snodes )
{		
	// get all the nodes in species tree
	//YW_ASSERT_INFO(false, "");
    treeSpecies.GetLeavesUnder( treeSpecies.GetRoot(), snodes );
}

int GeneSpeciesTreeHelper :: GetTaxaAtLeafST(int snid)
{
    if( treeSpecies.IsLeaf(snid) == false )
    {
        return -1;
    }
//cout << "GetTaxaAtLeafST: " << snid << endl;
	// get the taxa id from a particular leaf node of species tree
	// treat the species tree label
	//YW_ASSERT_INFO( treeSpecies.IsLeaf(snid) == true, "Must be a leaf1" );
	return treeSpecies.GetLabel(snid);
}


void GeneSpeciesTreeHelper :: GetLeavesSpeciesGT( set<int> &species )
{
	// get all the species/labels in gene tree
	// 
	GetLeavesSpeciesFromGT(treeGene, species);
}


void GeneSpeciesTreeHelper :: GetGeneAllelesForSpecies( int taxon, set<int> &geneAlleles)
{
	geneAlleles.clear();
//cout << "GetGeneAllelesForSpecies: taxon = " << taxon << endl;
	// get all alleles in GT for a gene
	char idbuf[100];
	snprintf(idbuf, 100, "%d", taxon);
	string idstr = idbuf;
	// search for each leaf, to find match
	treeGene.GetLeavesIdsWithLabel(  idstr, geneAlleles );
}

int GeneSpeciesTreeHelper :: GetSpeciesForGTLeaf(int lin) const
{
    //
    YW_ASSERT_INFO(mapIdNode.find(lin) != mapIdNode.end(), "Cannot find the node");
    GeneSpeciesTreeHelper &pthis = const_cast<GeneSpeciesTreeHelper &>(*this);
    TreeNode *pGT = pthis.mapIdNode[lin];
    YW_ASSERT_INFO(pGT->IsLeaf() == true, "Must be a leaf");
    string lbl = pGT->GetLabel();
    int res = -1;
    sscanf( lbl.c_str(), "%d", &res );
    YW_ASSERT_INFO( res >=0, "Fail" );
    return res;
}

void GeneSpeciesTreeHelper :: GetGeneAllelesForSTNode( int stNode, set<int> &geneAlleles, const map<int,int> *pMapOldLeavesToNew )
{
    geneAlleles.clear();
    // get all alleles from species under a specific ST branch
    // first get all leaves in ST under stNode
    set<int> stLeaves;
    treeSpecies.GetLeavesUnder( stNode, stLeaves);
    if( pMapOldLeavesToNew != NULL)
    {
        set<int> stLeavesConv;
        for(set<int> :: iterator it = stLeaves.begin(); it != stLeaves.end(); ++it)
        {
            YW_ASSERT_INFO( pMapOldLeavesToNew->find(*it) != pMapOldLeavesToNew->end(), "Fail to find" );
            stLeavesConv.insert( (*( pMapOldLeavesToNew->find(*it) )).second );
        }
        stLeaves = stLeavesConv;
    }
//cout << "GetGeneAllelesForSTNode: Under stnode: " << stNode << ", stLeaves: ";
//DumpIntSet(stLeaves);
    for(set<int> :: iterator it = stLeaves.begin(); it != stLeaves.end(); ++it)
    {
        set<int> gAlleles;
        GetGeneAllelesForSpecies(*it, gAlleles);
        UnionSets( geneAlleles, gAlleles );
    }
//cout << "Set of gene alleles: ";
//DumpIntSet( geneAlleles );
}


// Create a linConfig
void GeneSpeciesTreeHelper :: FormLinConfig( int snid, const set<int> &linsGT, LineageConfig &linConfig ) const
{
	linConfig.Clear();

	// generate trivial lineage cluster of single lineage
	// find whether there are lineages that can be merged
	linConfig.SetNodeId(snid);
	for( set<int> :: iterator it = linsGT.begin(); it != linsGT.end(); ++it  )
	{
		// 
		int nid = *it;
		int nidpar = (const_cast<GeneSpeciesTreeHelper *>(this))->GetNodeGTPar(nid);
		LineageCluster lc(nid, nidpar);
		linConfig.AddLC(lc);
	}
//	linConfig.Consolidate(*this);
}

void GeneSpeciesTreeHelper :: FormLinConfig( int snid, const LineageConfig &linsGT1, const LineageConfig &linsGT2, LineageConfig &linConfig ) const
{
	// merging two configs into one
	linConfig.Clear();

	// generate trivial lineage cluster of single lineage
	// find whether there are lineages that can be merged
	linConfig.SetNodeId(snid);
	set< LineageCluster > setLCs1;
	linsGT1.GetLCs(setLCs1);
	set< LineageCluster > setLCs2;
	linsGT2.GetLCs(setLCs2);
	// add to it
	for( set< LineageCluster > :: iterator it=  setLCs1.begin(); it != setLCs1.end(); ++it )
	{
		linConfig.AddLC(*it);
	}
	for( set< LineageCluster > :: iterator it=  setLCs2.begin(); it != setLCs2.end(); ++it )
	{
		linConfig.AddLC(*it);
	}

//	linConfig.Consolidate(*this);
}

int GeneSpeciesTreeHelper :: GetNodeGTPar(int gnid) const
{
	YW_ASSERT_INFO( mapNodeParGT.find( gnid ) != mapNodeParGT.end(), "Can not find gnid's parent" );
	return const_cast<GeneSpeciesTreeHelper *>(this)->mapNodeParGT[gnid];
}

bool GeneSpeciesTreeHelper :: AreLinsHalfSibling( const LineageCluster &lin1, const LineageCluster &lin2  ) const
{
	int idp1 = GetNodeGTPar( lin1.GetGNId() );
	if( idp1 < 0 )
	{
		return false;
	}
	int idp2 = GetNodeGTPar( lin2.GetGNId() );
	if( idp2 < 0 )
	{
		return false;
	}
	return idp1 != idp2 && GetNodeGTPar(idp1) == GetNodeGTPar(idp2);
}


int GeneSpeciesTreeHelper :: GetMRCATwoLineagesGT( int lin1, int lin2 ) const
{
	// get MRCA for two lineages (represnted by the starting node, in backward in time sense)
	YW_ASSERT_INFO( lin1 < (int)vecMRCATwoLinsGT.size() && lin2 < (int)vecMRCATwoLinsGT.size(), "range error"  );
	int mrca = vecMRCATwoLinsGT[lin1][lin2];
	YW_ASSERT_INFO( mrca >= 0, "MRCA badly initialied" );
	return mrca;
}

int GeneSpeciesTreeHelper :: GetParentLineage(int lin) const
{
	// 
	YW_ASSERT_INFO( lin < (int)vecMRCATwoLinsGT.size(), "range error"  );
	if( mapNodeParGT.find( lin ) == mapNodeParGT.end() )
	{
		// not found: meaning it is root?
		return -1;
	}
	return const_cast<GeneSpeciesTreeHelper &>(*this).mapNodeParGT[lin];
}

void GeneSpeciesTreeHelper :: FindCoalEventsBtwConfigs( const LineageConfig &linCfgSrc, const LineageConfig &linCfgDest, vector<CoalEvent> &listEvts  ) const
{
	// between two lineage configs, what has happened?
	set<int> nodeIdsSrc;
	linCfgSrc.GetSrcNodes(nodeIdsSrc);
	map<int,int> mapLinInfoSrc;
	linCfgSrc.GetSrcNodesNumLins( mapLinInfoSrc );
	set<int> nodeIdsDest;
	linCfgDest.GetSrcNodes(nodeIdsDest);
	map<int,int> mapLinInfoDest;
	linCfgDest.GetSrcNodesNumLins( mapLinInfoDest );
//cout << "FindCoalEventsBtwConfigs: nodeIdsSrc: ";
//DumpIntSet( nodeIdsSrc );
//cout << "nodeIdsDest: ";
//DumpIntSet( nodeIdsDest );

	// exam each node in src and trace backwards until reaching one node in dest
	set< pair<int,int> > nodesDone;
	FindDoneLineages( nodeIdsSrc, mapLinInfoSrc, nodeIdsDest, mapLinInfoDest, nodesDone );
	//for( set<int> :: iterator it = nodeIdsSrc.begin(); it != nodeIdsSrc.end(); ++it )
	//{
		// 
	//	FindDoneLineages( *it, nodeIdsDest, nodesDone );
	//}
//cout << "nodesDone: ";
//for(set< pair<int,int> > :: iterator itt = nodesDone.begin(); itt != nodesDone.end(); ++itt)
//{
//cout << "(" << itt->first << "," << itt->second << "),";
//}
//cout << endl;
//DumpIntSet( nodesDone );

	// Then figure out events: for each node v, add its 
	// for each done lineage, if it is not in dest cfg, then it has NUM-DEGREE - 1 
	// events associated with it; if it is in dest cfg, we need to check the
	// difference of lineage number
	for( set< pair<int,int> > :: iterator it = nodesDone.begin(); it != nodesDone.end(); ++it)
	{
		int ntid = it->first;
		int numEvt = it->second;
		// Create events
		for(int ii=0; ii<numEvt; ++ii)
		{
			CoalEvent cevt( ntid );
			listEvts.push_back( cevt );
//cout << "Find one done coalecent event at: " << ntid << endl;
		}

	}
}


void GeneSpeciesTreeHelper :: FindDoneLineages( const set<int> &nodeIdsSrc, const map<int,int> &srcLinNums, 
											   const set<int> &nodeIdsDest, const map<int,int> &destLinNums, set< pair<int,int> > &nodesDone ) const
{
	map<int,pair<int,int>  > mapDoneLins;
	nodesDone.clear();

	stack<int> nodesToUp;
	// init by the srce nodes
	for( set<int> :: iterator it = nodeIdsSrc.begin(); it != nodeIdsSrc.end(); ++it )
	{
		//int ncurr = *it;
		nodesToUp.push(*it);
	}

	// idea: trace backwards from src, and for each new node, figure out how many lineages we need to merge
	// first, process whatever we are having now
	while(  nodesToUp.empty() == false )
	{
		int ncurr = nodesToUp.top();
		nodesToUp.pop();
		int ncurrPar = GetParentLineage(ncurr);
		if( ncurrPar < 0 )
		{
			// for the root lineage, do nothing
			continue;
		}
		// if this node is in dest, skip
		if( nodeIdsDest.find(ncurr) != nodeIdsDest.end() )
		{
			continue;
		}
//cout << "FindDoneLineages: ncurr = " << ncurr << ", ncurrPar = " << ncurrPar << endl;
		int numLinsCurr = 1;
		if( srcLinNums.find( ncurr ) !=srcLinNums.end() )
		{
			numLinsCurr = const_cast<map<int,int> &>(srcLinNums)[ncurr];
		}
		int coalEvtDone = numLinsCurr -1;
		if(mapDoneLins.find(ncurrPar) == mapDoneLins.end() )
		{
			pair<int,int> pp(0, 0);
			mapDoneLins.insert( map<int,pair<int,int> > :: value_type(ncurrPar,pp) );
		}
		mapDoneLins[ncurrPar].first += numLinsCurr;
		mapDoneLins[ncurrPar].second += coalEvtDone;

		// if this node is done, we will follow up on it
		int nparNumDescs = GetGTNodeNumChildren( ncurrPar );
		if( mapDoneLins[ncurrPar].first == nparNumDescs )
		{
			// this node is done. 
			nodesToUp.push( ncurrPar );

			// add a record to it
			pair<int,int> pp1( ncurrPar, nparNumDescs - 1 - mapDoneLins[ncurrPar].second  );
			YW_ASSERT_INFO( pp1.second >=0, "num evts wrong"  );
			nodesDone.insert( pp1 );
		}
	}

}

int GeneSpeciesTreeHelper :: GetGTNodeNumChildren(int gnid) const
{
	YW_ASSERT_INFO( mapIdNode.find( gnid ) != mapIdNode.end(), "Fail to find the node" );
	return const_cast<GeneSpeciesTreeHelper *>(this)->mapIdNode[gnid]->GetChildrenNum();
}

double GeneSpeciesTreeHelper :: CalcCoalCoeffForTwoCfg( const LineageConfig &linCfgSrc, const LineageConfig &linCfgDest ) const
{
    if( fVerboseMode == true )
    {
        cout << "CalcCoalCoeffForTwoCfg: src cfg = ";
        linCfgSrc.Dump();
        cout << "CalcCoalCoeffForTwoCfg: dest cfg = ";
        linCfgDest.Dump();
    }
	// given two configurations, what is the coefficient term for the combinatorial factor
	// used in probability computation (see Degnan's paper)
	// first get the set of events
	vector<CoalEvent> listEvts;
	FindCoalEventsBtwConfigs(linCfgSrc, linCfgDest, listEvts );
	int cb = listEvts.size();
	int ub = linCfgSrc.GetTotNumLins();
	//double res = 1.0/CalcdbVal(ub, cb);
	double res = 1.0;
    //double coeff = 1.0;
    //if( GSTPIsHeuModeOn() == true )
    //{
    //    res *= CalcwbValMultiHeu( listEvts, linCfgSrc, ub, cb );
    //}
    //else
    //{
    if( IsGeneTreeMulfurcate() == false)
    {
        res *= CalcwbVal(listEvts, linCfgSrc, ub, cb);
    //cout << "CalcCoalCoeffForTwoCfg: cb = " << cb << ", ub = " << ub << ", coeff = " << res << endl;
    }
    else
    {
        res = 1.0/CalcdbVal(ub, cb);
        res *= const_cast<GeneSpeciesTreeHelper *>(this)->CalcwbValMulti( linCfgSrc, linCfgDest);
    }
    //}
	return res;
}

double GeneSpeciesTreeHelper :: CalcdbVal(int ub, int cb) const
{
//if( ub < cb + 1)
//{
//cout << "CalcdbVal: ub=" << ub << ", cb=" << cb << endl;
//}
	YW_ASSERT_INFO(  ub >= cb+1, "Can not coalescent more than the number of lineages" );

	// cahce to avoid redundnet computation
	static map< pair<int,int>, double> mapdbVals;
	pair<int,int> pp(ub,cb);
	if(mapdbVals.find( pp) != mapdbVals.end() )
	{
		return mapdbVals[pp];
	}

	// this is the number of all possible coalescents for ub lineages that will undergo cb coalescents
	double res = 1.0;
	for( int i=0; i<cb; ++i )
	{
		res *= 0.5*(ub-i)*(ub-i-1);
	}
	mapdbVals.insert( map< pair<int,int>, double> :: value_type(pp, res) );
	return res;
}

double GeneSpeciesTreeHelper :: CalcProddbValAndCoalFac(int ub, int cb) const
{
    // different from the above, calculate the term of 1/db times cb!
    YW_ASSERT_INFO(  ub >= cb+1, "Can not coalescent more than the number of lineages" );
    
	// cahce to avoid redundnet computation
	//static map< pair<int,int>, double> mapdbVals;
	//pair<int,int> pp(ub,cb);
	//if(mapdbVals.find( pp) != mapdbVals.end() )
	//{
	//	return mapdbVals[pp];
	//}
    
	// this is the number of all possible coalescents for ub lineages that will undergo cb coalescents
	double res = 1.0;
	for( int i=0; i<cb; ++i )
	{
        res *= i+1;
		res /= 0.5*(ub-i)*(ub-i-1);
	}
	//mapdbVals.insert( map< pair<int,int>, double> :: value_type(pp, res) );
	return res;
}

double GeneSpeciesTreeHelper :: CalcLogProddbValAndCoalFac(int ub, int cb) const
{
    // different from the above, calculate the term of 1/db times cb!
    YW_ASSERT_INFO(  ub >= cb+1, "Can not coalescent more than the number of lineages" );
    
	// cahce to avoid redundnet computation
	//static map< pair<int,int>, double> mapdbVals;
	//pair<int,int> pp(ub,cb);
	//if(mapdbVals.find( pp) != mapdbVals.end() )
	//{
	//	return mapdbVals[pp];
	//}
    
	// this is the number of all possible coalescents for ub lineages that will undergo cb coalescents
	double res = 0.0;
	for( int i=0; i<cb; ++i )
	{
        double resstep = i+1;
		resstep /= 0.5*(ub-i)*(ub-i-1);
        res += log(resstep);
	}
	//mapdbVals.insert( map< pair<int,int>, double> :: value_type(pp, res) );
	return res;
}

// calculate the coalescent coefficient for a set of coalescent evets
double GeneSpeciesTreeHelper :: CalcBranchFactorAt( int numDupEvts, const set<int> &sDescs, int numActiveLins, map<int,int> &mapEvtsUnder )
{
//cout << "***CalcBranchFactorAt: numDupEves = " << numDupEvts  << ", numActiveLins = " << numActiveLins << ", sDescs = ";
//DumpIntSet( sDescs );

	YW_ASSERT_INFO( numDupEvts >= 1, "Can not have zero or negative events in branch factor computation" );
	YW_ASSERT_INFO( numActiveLins >= 2, "Must have at least two lineages in coalescents" );


	// when there is no underlying evets
	if( sDescs.size() == 0 )
	{
		//return GetMultiNoDescCoeff(numDupEvts);
		return GetMultiNoDescCoeff(numActiveLins, numDupEvts);
	}


	static map< multiset<int>, double > mapNodeFactors;

	// figure out num of dup events
	//int numDupEvts = 0;

	// compute for each node, what factor we need to compute
	// get the descendent nodes desc num
	multiset<int> sDescNums;
	for( set<int> :: iterator it = sDescs.begin(); it != sDescs.end(); ++it )
	{
		//if( *it >= 0)
		//{
		// 
		YW_ASSERT_INFO( mapEvtsUnder.find(*it) != mapEvtsUnder.end(), "Fatal error in CalcBranchFactorAt" );
		//if( mapEvtsUnder.find(*it) != mapEvtsUnder.end() )
		//{
		sDescNums.insert( mapEvtsUnder[*it] );
		//}
		//else
		//{
		//	// these are non-descendent nodes
		//	sDescNums.insert( 0 );
		//}
	}
	// if there are duplicate events, then we also need to add some of these invisible lineages (with 0 events down)
	//if( (int)sDescNums.size() < numDupEvts )
	//{
	//	for(int i=0; i<  numDupEvts - (int)sDescNums.size(); ++i)
	//	{
	//		// add dummy desc events here
	//		sDescNums.insert( 0 );
	//	}
	//}

//cout << "In CalcBranchFactorAt: branch desc num = ";
//for( multiset<int> :: iterator it = sDescNums.begin(); it != sDescNums.end(); ++it )
//{
//cout << " " << *it;
//}
//cout << endl;


	// now recurisvely try all possible way of dividing the events
	double res = 1.0;
	
	// now with more than 2 descendents, need to consider all ways of partition them and then combine
	vector<int> listDescs;
	//int numTotDesc = 0;
	for( set<int> :: iterator it = sDescs.begin(); it != sDescs.end(); ++it )
	{
		listDescs.push_back( *it );
		//numTotDesc += *it;
	}
	int numTotDesc = 0;
	for( multiset<int> :: iterator it = sDescNums.begin(); it != sDescNums.end(); ++it )
	{
		numTotDesc += *it;
	}
//cout << "Num of tot desc = " << numTotDesc << endl;;
//cout << "numDupEvts = " << numDupEvts << endl;;


#if 0
	if( numActiveLins == 2 )
	//if( numDupEvts == 1 )
	//if( sDescNums.size() == 0 )
	{
		// what if we have multiple events here (although no other descendents)
		// collect num of descendent evets
		YW_ASSERT_INFO( sDescNums.size() <=2, "Can not have more than two descendents" );
		//int numDescEvts = 0;
		//for(multiset<int> :: iterator itt =  sDescNums.begin(); itt != sDescNums.end(); ++itt)
		//{
		//	numDescEvts += *itt;
		//}
cout << "CalcBranchFactorAt: simple return\n";
		return 1.0/(1.0+numTotDesc);
	}
#endif

	//------------------------------------------------------------------
	// from now on, we have at least two coalescent events

	//if( sDescNums.size() == 1 )
	//{
	//	// just return the factor directly
	//	int numNodesUnder = *( sDescNums.begin() );
	//	YW_ASSERT_INFO(numNodesUnder >=1, "Can not have zero events");
	//	return 1.0/(1.0+ numNodesUnder );
	//}
	//if( sDescNums.size() == 2 )
	//{
	//	// in this case, can not have more than one event at the parent node
	//	// just return the factor directly
	//	multiset<int> :: iterator itt =  sDescNums.begin();
	//	int numNodesUnder1 = *itt;
	//	itt++;
	//	int numNodesUnder2 = *itt;
	//	int evtTot = numNodesUnder1 + numNodesUnder2;
//cout << "Two descent: evtTot = " << evtTot << endl;
	//	// there can only be a single event at this coalescent node
	//	return 1.0/( 1.0+ evtTot);
	//}
	// if it is already computed before, just use it
	if( mapNodeFactors.find( sDescNums ) != mapNodeFactors.end() )
	{
		return mapNodeFactors[ sDescNums ];
	}



	// deal with the first (topmost) coalescent events
	//int topEventsUnder = sDescs.size()-1 + numTotDesc;
	int topEventsUnder = numDupEvts-1 + numTotDesc;
//cout << "Topmost coalescent event has this many events under: " << topEventsUnder << endl;
	YW_ASSERT_INFO( topEventsUnder >= 1, "Must be 1 or more" );
	res *= 1.0/( 1+topEventsUnder );

	YW_ASSERT_INFO(numDupEvts >=(int)listDescs.size()-1, "Must be the case" );
//	YW_ASSERT_INFO(numDupEvts>1, "numDupEvts got to be at least two here");
	int numNoDescLins = numActiveLins - listDescs.size() ;
	YW_ASSERT_INFO(numNoDescLins>=0, "numNoDescEvts got to be at least zero here");
//cout << "numNoDescLins = " << numNoDescLins << endl;

#if 0
	// now append these non-descendents (with zero coalescents in them)
	for(int i=0; i<numNoDescLins; ++i)
	{
		listDescs.push_back(-1*(i+1));
	}

	double cooeffTot = 0.0;
    for( int sz = 0; sz<=(int)listDescs.size()-1; ++sz  )
    {
		// sz: num of events on left - 1 (because we have fixed one event to this side, the last one on the list)
//cout << "dmNum = " << dmNum << ", probDS = " << probDS << endl;
        // now try next one
		//int szDescSelect = sz;
		// need to reduce if we do not have that many descendent events
		//if(  szDescSelect >= listDescs.size()-1 )
		//{
		//	szDescSelect = ;
		//}
        vector<int> posvec;
        GetFirstCombo( sz, listDescs.size()-1, posvec );

        while(true)
        {
            // to avoid double-counting, require the first position to be 0
			// now form two parts: one with the selected items (with the last item)
			// the other are those not in the partition
			set<int> stmp;
			set<int> sint1, sint2;
            PopulateSetByVec( stmp, posvec );
			YW_ASSERT_INFO( stmp.find(listDescs.size()-1) == stmp.end(), "The last item should not be in yet" );
			stmp.insert(  listDescs.size()-1  );
			for( set<int> :: iterator itg = stmp.begin(); itg != stmp.end(); ++itg)
			{
				sint1.insert( listDescs[*itg] );
			}
			// now add the last one
			sint2 = sDescs;
			SubtractSets( sint2, sint1 );
//cout << "sint1: ";
//DumpIntSet(sint1);
//cout << "sint2: ";
//DumpIntSet(sint2);

			// now calculate for this option
			double fac1 = 1.0;
			if( sint1.size() > 1 )
			{
				// have a factor for it
				fac1 = CalcBranchFactorAt(sint1.size()-1, sint1, mapEvtsUnder);
//cout << "fac1 = " << fac1 << endl;
				//res *= fac1;
			}
			double fac2 = 1.0;
			if( sint2.size() > 1 )
			{
				// have a factor for it
				fac2 = CalcBranchFactorAt(sint2.size()-1,  sint2, mapEvtsUnder);
//cout << "fac2 = " << fac2 << endl;
				//res *= fac2;
			}
			cooeffTot +=  fac1 * fac2;			;

            // get next mutation edge set
            if( GetNextCombo( sz, listDescs.size()-1, posvec ) == false )
            {
                break;
            }

		}
	}
#endif	


//#if 0
	// we  now have remaining evt-1 to deal with, 
	// we do this by first spliting events w/ descendents
	double cooeffTot = 0.0;
    for( int sz = 0; sz<=(int)listDescs.size()-1; ++sz  )
    {
		// sz: num of events on left - 1 (because we have fixed one event to this side, the last one on the list)
//cout << "dmNum = " << dmNum << ", probDS = " << probDS << endl;
        // now try next one
		//int szDescSelect = sz;
		// need to reduce if we do not have that many descendent events
		//if(  szDescSelect >= listDescs.size()-1 )
		//{
		//	szDescSelect = ;
		//}
        vector<int> posvec;
        GetFirstCombo( sz, listDescs.size()-1, posvec );

        while(true)
        {
            // to avoid double-counting, require the first position to be 0
			// now form two parts: one with the selected items (with the last item)
			// the other are those not in the partition
			set<int> stmp;
			set<int> sint1, sint2;
            PopulateSetByVec( stmp, posvec );
			YW_ASSERT_INFO( stmp.find(listDescs.size()-1) == stmp.end(), "The last item should not be in yet" );
			stmp.insert(  listDescs.size()-1  );
			for( set<int> :: iterator itg = stmp.begin(); itg != stmp.end(); ++itg)
			{
				sint1.insert( listDescs[*itg] );
			}
			// now add the last one
			sint2 = sDescs;
			SubtractSets( sint2, sint1 );
//cout << "sint1: ";
//DumpIntSet(sint1);
//cout << "sint2: ";
//DumpIntSet(sint2);

			// now, do we have choices in spliting the number of remaining hidden events?
			// if so, we should do it here
			for(int hvLeft = 0; hvLeft <=numNoDescLins; ++hvLeft)
			{
				int hvRight = numNoDescLins-hvLeft;

				double combnumLeft = 1.0;
				if( hvLeft >= 1)
				{
					combnumLeft = CalcCombNum(numNoDescLins, hvLeft);
				}
				//double combnumRight = 1.0;
				//if( hvRight >= 1)
				//{
				//	combnumRight = CalcCombNum(numNoDescLins, hvRight);
				//}
//cout << "combnumLeft = " << combnumLeft << ", hvLeft = " << hvLeft << ", hvRight = " << hvRight << endl;

				// now another round: no. of events can also be different
				// this is because there may still be un-coalesced lineages there
				for(int numEvtLeft = 0; numEvtLeft <= numDupEvts; ++numEvtLeft)
				{
					// must -1 because the split itself is an event if all lineages are done in this subtree
					int numCoalEvtAtRoot = 0;
					if(numDupEvts+1 == numNoDescLins+(int)listDescs.size())
					{
						numCoalEvtAtRoot = 1;
					}
//cout << "numCoalEvtAtRoot: " << numCoalEvtAtRoot;
					int numEvtRight = numDupEvts - numEvtLeft-numCoalEvtAtRoot;
					//int numEvtRight = numDupEvts - numEvtLeft-1;
//cout << "  numEvtLeft: " << numEvtLeft << ", numEvtRight: " << numEvtRight << endl;

					// I think the left split is always allowed. But we need to ensure right split has enough lineages
					if( (numEvtLeft == numDupEvts && hvLeft ==numNoDescLins ) ||  numEvtLeft+1 > hvLeft+(int)sint1.size() || numEvtRight+1 > hvRight+(int)sint2.size() )
					{
						// not enough to satisfy this partition, so..
						continue;
					}
//cout << "Yes, this split is OK.\n";
					//int evtDupLeft = sint1.size()+hvLeft-1;
					//int evtDupRight = sint2.size()+hvRight-1;
	//cout << "numEvtLeft = " << numEvtLeft << ", numEvtRight = " << numEvtRight << endl;
					// make sure we are not repeating
					//if(  (evtDupLeft == numDupEvts) || ((evtDupRight == numDupEvts) ))
					//{
					//	continue;
					//}
					//YW_ASSERT_INFO(  (evtDupLeft <= numDupEvts) && (evtDupRight <= numDupEvts) , "Out of bound" );

					double fac1 = 1.0;
					if( numEvtLeft >= 1 )
					{
						// have a factor for it
						fac1 = CalcBranchFactorAt(numEvtLeft, sint1, sint1.size()+hvLeft, mapEvtsUnder);
//	cout << "fac1 = " << fac1 << endl;
						//res *= fac1;
					}
					double fac2 = 1.0;
					if( numEvtRight >= 1 )
					{
						// have a factor for it
						fac2 = CalcBranchFactorAt(numEvtRight,  sint2,  sint2.size()+hvRight, mapEvtsUnder);
//	cout << "fac2 = " << fac2 << endl;
						//res *= fac2;
					}

				cooeffTot +=  fac1 * fac2 * combnumLeft;
//cout << "cooeffTot = " << cooeffTot << endl;
				}
			}
            // get next mutation edge set
            if( GetNextCombo( sz, listDescs.size()-1, posvec ) == false )
            {
                break;
            }
        }
	}
//#endif
	YW_ASSERT_INFO(cooeffTot > 0.0, "cooeffTot must be positive");
	// now multiply the entire
	res *= cooeffTot;
	// store it: YW: do not do it just yet
	//mapNodeFactors.insert( map< multiset<int>, double > :: value_type( sDescNums, res ) ) ;
//cout << "res in CalcBranchFactorAt : " << res << endl;
	// update map
	return res;
}



double GeneSpeciesTreeHelper :: CalcwbVal( const vector<CoalEvent> &listEvts, const LineageConfig &lcSrc, int ub, int cb ) const
{
//cout << "CalcwbVal: num events = " << listEvts.size() << endl;
//for(int i=0; i<(int)listEvts.size();++i)
//{
//listEvts[i].Dump();
//}
	double res = 1.0;
	int numEvts = listEvts.size();
	for(int i=1; i<=numEvts; ++i)
	{
		res *= i;
	}
	if( ub >= 0 && cb >= 0)
	{
		for( int i=0; i<cb; ++i )
		{
			res *= 2.0/((ub-i)*(ub-i-1) );
		}
	}


	// now consider the constraints imposed by ordering
	// also keep track of duplicate events
	map<int,set<int> > mapDescBides;
	map<int,int> mapEvtsUnder;		// num of events under (or exactly at, i.e. duplicate) an event
	map<int,int> mapEvtsDupCount;	// for each event, how many copies do we have
	set<int> evtNodes;
	for( int i=0; i<(int)listEvts.size(); ++i )
	{
		int nid = listEvts[i].GetID();
		evtNodes.insert(nid);

		// add a record to it
		if( mapEvtsUnder.find( nid ) == mapEvtsUnder.end() )
		{
			mapEvtsUnder.insert( map<int,int> :: value_type(nid, 1));
		}
		else
		{ 
			mapEvtsUnder[nid]++;
		}

		int pnid =  GetParentLineage(nid);
		if( pnid >= 0)
		{
			// add a record. Note pnid may not be in the range of events, but that is fine
			if( mapDescBides.find( pnid ) == mapDescBides.end()  )
			{
				set<int> sint;
				sint.insert( nid );
				mapDescBides.insert( map<int,set<int> > :: value_type(pnid, sint) );
//cout << "Add descendent (init) lineage " << nid << " to pnid = " << pnid << endl;
			}
			else
			{
				mapDescBides[pnid].insert( nid );
//cout << "Add descendent (more) lineage " << nid << " to pnid = " << pnid << endl;
			}
		}
	}
	// so far, evtunder contains exactly those for dup count
	mapEvtsDupCount = mapEvtsUnder;
	// now figure out num of descendents (direct or indirect)
	// HERE WE ASSUME NODE IDS ARE ARRANGED S.T. SMALLER NODES MEANS THEY ARE CLOSER TO LEAVES
	// YW: THIS CAN CAUSE PROBLEM IF THIS CONVENTION IS CHANGED. CAUTION!!!!!
	for( set<int> :: iterator it = evtNodes.begin(); it != evtNodes.end(); ++it )
	{
		int nid = *it;
		YW_ASSERT_INFO( mapEvtsUnder.find( nid ) != mapEvtsUnder.end(), "Fail" );
		// update its
		//int linUnder = 1;		// 1: means the current node
		//if( mapDescBides.find(nid) != mapDescBides.end() )
		//{
		//	linUnder = mapDescBides[nid].size();
		//}
		// now update its par
		int pnid =  GetParentLineage(nid);
		if( pnid >= 0)
		{
			// add a record. Note pnid may not be in the range of events, but that is fine
			if( mapEvtsUnder.find( pnid ) != mapEvtsUnder.end()  )
			{
				//mapNodesUnder.insert( map<int,int> :: value_type(pnid, linUnder) );
				mapEvtsUnder[pnid] += mapEvtsUnder[nid];
			}
			//else
			//{
			//	mapNodesUnder[pnid] += linUnder;
			//}
		}
	}
#if 0
cout << "mapEvtsUnder = ";
for( map<int,int> :: iterator itt = mapEvtsUnder.begin(); itt != mapEvtsUnder.end();++itt )
{
cout << "Node [" << itt->first << ": " << itt->second << "], ";
}
cout << endl;
cout << "mapEvtsDupCount = ";
for( map<int,int> :: iterator itt = mapEvtsDupCount.begin(); itt != mapEvtsDupCount.end();++itt )
{
cout << "mapEvtsDupCount: Node [" << itt->first << ": " << itt->second << "], ";
}
cout << endl;
cout << "evtNodes: ";
DumpIntSet( evtNodes );
#endif

	// find out involved lineages in cfgsrc
	vector<int> listNoDupEvts;
	PopulateVecBySet( listNoDupEvts, evtNodes );
	vector<int> listLinsNumForEvt;
	FindNumLinsForEvtsCfg(listNoDupEvts, mapEvtsDupCount, lcSrc, listLinsNumForEvt);

	// for each node, compute a factor
	//for( set<int> :: iterator it = evtNodes.begin(); it != evtNodes.end(); ++it )
	for( int ievt = 0; ievt <(int)listNoDupEvts.size(); ++ievt )
	{
		int nid = listNoDupEvts[ievt];
		//int nid = *it;
		// do nothing if it does not have descendent
		if( mapDescBides.find(nid) != mapDescBides.end() )
		{
			YW_ASSERT_INFO( mapEvtsDupCount.find( nid ) != mapEvtsDupCount.end(), "Fail" );

			// caculate a combinatorial factor based on the num of descendents of each branches
			double fac = const_cast<GeneSpeciesTreeHelper *>(this)->CalcBranchFactorAt( mapEvtsDupCount[nid], mapDescBides[nid], 
				listLinsNumForEvt[ievt], mapEvtsUnder  );
//cout << "**fac = " << fac << ", for node " << nid << " with set of descendents: ";
//DumpIntSet(mapDescBides[nid]);
			res *= fac;
		}
		else if( mapEvtsDupCount[nid] > 1 || listLinsNumForEvt[ievt] >= 3 )
		{
			set<int> sempty;
			double fac = const_cast<GeneSpeciesTreeHelper *>(this)->CalcBranchFactorAt( mapEvtsDupCount[nid], sempty, 
				listLinsNumForEvt[ievt], mapEvtsUnder  );
//cout << "**(nodesc)fac = " << fac << ", for node " << nid << "with dup evts = " << mapEvtsDupCount[nid]   << endl;
			res *= fac;
		}

	}


//cout << "*************Final results = " << res << endl;
	return res;
}

double GeneSpeciesTreeHelper :: CalcwbValMultiHeu( const vector<CoalEvent> &listEvts, const LineageConfig &lcSrc, int ub, int cb ) const
{
cout << "In CalcwbValMultiHeu: listEvts: ";
for(int i=0;i<(int)listEvts.size(); ++i)
{
listEvts[i].Dump();
}
cout << "lcSrc: ";
lcSrc.Dump();
cout << ", ub= " << ub << ", cb= " << cb << endl;
    // heuristic version of the coefficient computation
    // approach: use a conservative estimate of the topological constraints
    double res = 1.0;
	int numEvts = listEvts.size();
	for(int i=1; i<=numEvts; ++i)
	{
		res *= i;
	}
	if( ub >= 0 && cb >= 0)
	{
		for( int i=0; i<cb; ++i )
		{
			res *= 2.0/((ub-i)*(ub-i-1) );
		}
	}
    
cout << "res so far: " << res << endl;
	// now consider the constraints imposed by ordering
	// also keep track of duplicate events
	map<int,set<int> > mapDescBides;
	map<int,int> mapEvtsUnder;		// num of events under (or exactly at, i.e. duplicate) an event
	map<int,int> mapEvtsDupCount;	// for each event, how many copies do we have
	set<int> evtNodes;
	for( int i=0; i<(int)listEvts.size(); ++i )
	{
		int nid = listEvts[i].GetID();
		evtNodes.insert(nid);
        
		// add a record to it
		if( mapEvtsUnder.find( nid ) == mapEvtsUnder.end() )
		{
			mapEvtsUnder.insert( map<int,int> :: value_type(nid, 1));
		}
		else
		{
			mapEvtsUnder[nid]++;
		}
        
		int pnid =  GetParentLineage(nid);
		if( pnid >= 0)
		{
			// add a record. Note pnid may not be in the range of events, but that is fine
			if( mapDescBides.find( pnid ) == mapDescBides.end()  )
			{
				set<int> sint;
				sint.insert( nid );
				mapDescBides.insert( map<int,set<int> > :: value_type(pnid, sint) );
                //cout << "Add descendent (init) lineage " << nid << " to pnid = " << pnid << endl;
			}
			else
			{
				mapDescBides[pnid].insert( nid );
                //cout << "Add descendent (more) lineage " << nid << " to pnid = " << pnid << endl;
			}
		}
	}
	// so far, evtunder contains exactly those for dup count
	mapEvtsDupCount = mapEvtsUnder;
	// now figure out num of descendents (direct or indirect)
	// HERE WE ASSUME NODE IDS ARE ARRANGED S.T. SMALLER NODES MEANS THEY ARE CLOSER TO LEAVES
	// YW: THIS CAN CAUSE PROBLEM IF THIS CONVENTION IS CHANGED. CAUTION!!!!!
	for( set<int> :: iterator it = evtNodes.begin(); it != evtNodes.end(); ++it )
	{
		int nid = *it;
		YW_ASSERT_INFO( mapEvtsUnder.find( nid ) != mapEvtsUnder.end(), "Fail" );
		// update its
		//int linUnder = 1;		// 1: means the current node
		//if( mapDescBides.find(nid) != mapDescBides.end() )
		//{
		//	linUnder = mapDescBides[nid].size();
		//}
		// now update its par
		int pnid =  GetParentLineage(nid);
		if( pnid >= 0)
		{
			// add a record. Note pnid may not be in the range of events, but that is fine
			if( mapEvtsUnder.find( pnid ) != mapEvtsUnder.end()  )
			{
				//mapNodesUnder.insert( map<int,int> :: value_type(pnid, linUnder) );
				mapEvtsUnder[pnid] += mapEvtsUnder[nid];
			}
			//else
			//{
			//	mapNodesUnder[pnid] += linUnder;
			//}
		}
	}
//#if 0
    cout << "mapEvtsUnder = ";
    for( map<int,int> :: iterator itt = mapEvtsUnder.begin(); itt != mapEvtsUnder.end();++itt )
    {
        cout << "Node [" << itt->first << ": " << itt->second << "], ";
    }
    cout << endl;
    cout << "mapEvtsDupCount = ";
    for( map<int,int> :: iterator itt = mapEvtsDupCount.begin(); itt != mapEvtsDupCount.end();++itt )
    {
        cout << "mapEvtsDupCount: Node [" << itt->first << ": " << itt->second << "], ";
    }
    cout << endl;
    cout << "evtNodes: ";
    DumpIntSet( evtNodes );
//#endif
    
	// find out involved lineages in cfgsrc
	vector<int> listNoDupEvts;
	PopulateVecBySet( listNoDupEvts, evtNodes );
	vector<int> listLinsNumForEvt;
	FindNumLinsForEvtsCfg(listNoDupEvts, mapEvtsDupCount, lcSrc, listLinsNumForEvt);
    
	// for each node, simply use the factor that is the number of events under it
	//for( set<int> :: iterator it = evtNodes.begin(); it != evtNodes.end(); ++it )
	for( int ievt = 0; ievt <(int)listNoDupEvts.size(); ++ievt )
	{
		int nid = listNoDupEvts[ievt];
        double fac = 1.0;
        if(  mapEvtsUnder.find(nid) != mapEvtsUnder.end() )
        {
            YW_ASSERT_INFO( mapEvtsDupCount.find( nid ) != mapEvtsDupCount.end(), "Fail" );
            int numLinsUnder = mapEvtsUnder[nid];
            numLinsUnder -= mapEvtsDupCount[nid];
            YW_ASSERT_INFO(numLinsUnder >= 0, "Fail222");
            fac = 1.0/(1.0+numLinsUnder);
        }
cout << "For event " << nid << ": fac = " << fac << endl;
        res *= fac;
	}
    
    
    //cout << "*************Final results = " << res << endl;
	return res;
}

// a utility code that is useful
double GeneSpeciesTreeHelper :: CalcMultiBranchComboNum( const vector< pair<int,double> > &listCounts ) const
{
//cout << "CalcMultiBranchComboNum: ";
//for(int i=0; i<(int)listCounts.size(); ++i)
//{
//cout << listCounts[i].first << ", " << listCounts[i].second << endl;
//}
	// format: <#evts,#ways>
	// finally combine into the total number
	double res = 1.0;
	int numEvtsStep = 0;

	for(int jj=0; jj<(int)listCounts.size(); ++jj)
	{
		numEvtsStep += listCounts[jj].first;
		int nFactor = CalcNumNChooseK( numEvtsStep, listCounts[jj].first );
		res *= listCounts[jj].second*nFactor;

		// while in the middle, if needed, also get the factor
		//if( jj == (int)(listCounts.size()/2) && ub >= 0 && cb >= 0 )
		//{
		//	res *= 1.0/CalcdbVal(ub, cb);
		//}
	}
	return res;
}


double GeneSpeciesTreeHelper :: CalcwbValMulti( const LineageConfig &linCfgSrc, const LineageConfig &linCfgDest )
{
	// to speedup, we keep a cache of pre-computed wbval, and if we have already computed it before, then we just look it up
	multiset<int> msSrc, msDest;
	linCfgSrc.GetParLins(msSrc);
	linCfgDest.GetParLins(msDest);
	pair<multiset<int>, multiset<int> > ppSrcDest(msSrc, msDest);
static int numQuickRef = 0;
	if(cacheWbValMulti.find(ppSrcDest) != cacheWbValMulti.end() )
	{
if( (++numQuickRef % 1000 ) == 0)
{
//cout << "Quick access of CalcwbValMulti num = " << numQuickRef << endl;
}
		return cacheWbValMulti[ppSrcDest];
	}
//cout << "CalcwbValMulti: src cfg = ";
//linCfgSrc.Dump();
//cout << "CalcwbValMulti: dest cfg = ";
//linCfgDest.Dump();
	// for trees with mulfurcating nodes
	// idea: bottom up; at each internal node, consider all ways of coalesces
	// for a binary coalescence with N1 valid ways and N2 valid ways and each with k1 and k2 coalescent (plus the
	// new coalescent event, then the number of ways for the combined subtree is equal to C(k1, k1+k2)*N1*N2
	// when there are different combinations, we just add them up wrt different combinations
	set<int> nodeIdsSrc;
	linCfgSrc.GetSrcNodes(nodeIdsSrc);
	map<int,int> mapLinInfoSrc;
	linCfgSrc.GetSrcNodesNumLins( mapLinInfoSrc );
	set<int> nodeIdsDest;
	linCfgDest.GetSrcNodes(nodeIdsDest);
	map<int,int> mapLinInfoDest;
	linCfgDest.GetSrcNodesNumLins( mapLinInfoDest );
//cout << "nodeIdsDest = ";
//DumpIntSet( nodeIdsDest );

	// what nodes are currently active, and for each node, how many underlying events and number
	// of ways to arrange them: <nodeid, <#evts and #ways> >

	map<int,pair<int,int>  > mapDoneLins;
	map<int,set<int> > mapDestLinsInfo;			// what are the parents of destination lineages
	map<int, set<int> > mapDoneLinsParInfo;		// which parent's descendents are done
	map<int,NODE_COUNT_INFO > mapDoneLinsCountInfo;	// for each parent lineage, what descendents are beneath them?

	// what nodes have been absorbed
	set<int> setNodesCounted;

	// stuff to keep track during backward traversal

	stack<int> nodesToUp;
	// init by the srce nodes
	for( set<int> :: iterator it = nodeIdsSrc.begin(); it != nodeIdsSrc.end(); ++it )
	{
		//int ncurr = *it;
		//EVT_TRACE_INFO cfgInit;
		//cfgInit.nodeId = *it;
		//cfgInit.listDescendents.clear();
		nodesToUp.push( *it );

		// consider each leaf node done
		NODE_COUNT_INFO ncInfo;
		ncInfo.numEvtsBelow = 0;
		ncInfo.numWays = 1;
		mapDoneLinsCountInfo.insert( map<int,NODE_COUNT_INFO > :: value_type(*it, ncInfo) );
	}

	// idea: trace backwards from src, and for each new node, figure out how many lineages we need to merge
	// first, process whatever we are having now
	while(  nodesToUp.empty() == false )
	{
		//EVT_TRACE_INFO nodeInfoCur = nodesToUp.top();
		int ncurr = nodesToUp.top();
		nodesToUp.pop();
		int ncurrPar = GetParentLineage(ncurr);

		// record it
		if( nodeIdsDest.find(ncurr) == nodeIdsDest.end()  )
		{
			if( mapDoneLinsParInfo.find( ncurrPar) == mapDoneLinsParInfo.end() )
			{
				set<int> sdummy;
				mapDoneLinsParInfo.insert( map<int, set<int> > :: value_type(ncurrPar, sdummy) );
			}
			mapDoneLinsParInfo[ncurrPar].insert( ncurr );
		}

		if( ncurrPar < 0 )
		{
			// for the root lineage, do nothing
			continue;
		}
		// if this node is in dest, skip
		if( nodeIdsDest.find(ncurr) != nodeIdsDest.end() )
		{
			// also remember this info
			if( mapDestLinsInfo.find(ncurrPar) == mapDestLinsInfo.end() )
			{
				set<int> sempty;
				mapDestLinsInfo.insert( map<int,set<int> > :: value_type(ncurrPar, sempty) );
			}
			mapDestLinsInfo[ncurrPar].insert(ncurr);

			// hey, what does this mean?? Not sure. 2/12/11
			continue;
		}
//cout << "CalcwbValMulti: ncurr = " << ncurr << ", ncurrPar = " << ncurrPar << endl;
		int numLinsCurr = 1;
		if( mapLinInfoSrc.find( ncurr ) !=mapLinInfoSrc.end() )
		{
			numLinsCurr = mapLinInfoSrc[ncurr];
		}
		int coalEvtDone = numLinsCurr -1;
		if(mapDoneLins.find(ncurrPar) == mapDoneLins.end() )
		{
			pair<int,int> pp(0, 0);
			mapDoneLins.insert( map<int,pair<int,int> > :: value_type(ncurrPar,pp) );
		}
		if(mapDoneLinsCountInfo.find(ncurrPar) == mapDoneLinsCountInfo.end() )
		{
			NODE_COUNT_INFO plist;
			mapDoneLinsCountInfo.insert( map<int,NODE_COUNT_INFO> :: value_type(ncurrPar,plist) );
		}

		mapDoneLins[ncurrPar].first += numLinsCurr;
		mapDoneLins[ncurrPar].second += coalEvtDone;
		// rememer the child lineages
		mapDoneLinsCountInfo[ncurrPar].lineagesDesc.push_back( ncurr );
//cout << "for mapDoneLinsCountInfo: add " << ncurr <<"; after adding lineageDesc = ";
//DumpIntVec(mapDoneLinsCountInfo[ncurrPar].lineagesDesc);
		// if this node is done, we will follow up on it
		int nparNumDescs = GetGTNodeNumChildren( ncurrPar );
//cout << "Num of children of this par node: " << nparNumDescs << ", mapDoneLins[ncurrPar].first = " << mapDoneLins[ncurrPar].first  << endl;
		if( mapDoneLins[ncurrPar].first == nparNumDescs )
		{
			// when a parent node is finished, mark off all its children as counted
			for(int tt=0; tt<(int)mapDoneLinsCountInfo[ncurrPar].lineagesDesc.size(); ++tt)
			{
				int nodeIdChild = mapDoneLinsCountInfo[ncurrPar].lineagesDesc[tt];
				YW_ASSERT_INFO( mapDoneLinsCountInfo.find(nodeIdChild) != mapDoneLinsCountInfo.end(), "Fail to find children" );
				setNodesCounted.insert(nodeIdChild);
			}

			// this node is done. 
			nodesToUp.push( ncurrPar );

			// add a record to it
			//pair<int,int> pp1( ncurrPar, nparNumDescs - 1 - mapDoneLins[ncurrPar].second  );
			//YW_ASSERT_INFO( pp1.second >=0, "num evts wrong"  );
			//nodesDone.insert( pp1 );

			// now count how many events are there
			int numEvtsThisNode = nparNumDescs - 1 - mapDoneLins[ncurrPar].second;
			YW_ASSERT_INFO(numEvtsThisNode >=0, "Can not be negative");
			int numEvtsTot = numEvtsThisNode;
			vector< pair<int,double> > listDescEvtWaysCounts;
			for( int tt=0; tt<(int)mapDoneLinsCountInfo[ncurrPar].lineagesDesc.size(); ++tt )
			{
				int nodeIdChild = mapDoneLinsCountInfo[ncurrPar].lineagesDesc[tt];
				YW_ASSERT_INFO( mapDoneLinsCountInfo.find(nodeIdChild) != mapDoneLinsCountInfo.end(), "Fail to find children" );
				numEvtsTot += mapDoneLinsCountInfo[nodeIdChild].numEvtsBelow;
				pair<int,double> pp( mapDoneLinsCountInfo[nodeIdChild].numEvtsBelow, mapDoneLinsCountInfo[nodeIdChild].numWays );
				listDescEvtWaysCounts.push_back( pp );
			}
//cout << "numEvtsTot = " << numEvtsTot << endl;
			YW_ASSERT_INFO( mapDoneLinsCountInfo.find(ncurrPar) != mapDoneLinsCountInfo.end(), "Fail1" );
			mapDoneLinsCountInfo[ncurrPar].numEvtsBelow = numEvtsTot;

			// now need to update on the number of ways. Need to do this through a recurisvely procedure (implemented in DP)
			mapDoneLinsCountInfo[ncurrPar].numWays = CountParFromChildCoal(listDescEvtWaysCounts, 1);
//cout << "apDoneLinsCountInfo[ncurrPar].numWays = " << mapDoneLinsCountInfo[ncurrPar].numWays << endl;
		}
	}

	// find nodes that are in the destination (and their number of events at that)
	// pair: <parent node id, list of dest edges to that node>
	// this is like clustering destination lineages based on their parents
//cout << "Again, nodeIdsDest = ";
//DumpIntSet( nodeIdsDest );
	map<int, set< int> > setNodesDestInfo;
	for( set<int> :: iterator itt = nodeIdsDest.begin(); itt != nodeIdsDest.end(); ++itt )
	{
		int parlin = GetParentLineage(*itt);
		if( setNodesDestInfo.find( parlin ) == setNodesDestInfo.end() )
		{
			set<int> sdummy;
			setNodesDestInfo.insert(map<int, set< int> >:: value_type(parlin, sdummy) );
		}
		setNodesDestInfo[parlin].insert(*itt);
//cout << "setNodesDestInfo: add an entry " << parlin << ", " << *itt;
//cout << ", update: ";
//DumpIntSet( setNodesDestInfo[parlin] );
	}

	// need to wrap up: there may two coalescable but not yet coalesced and there may be several clusters
	// for example, there are two lineages should coalesce but there are three lineages in total
	// so check each record in the mapped nodes, and collect all those whose nodes are not yet finished
	vector< pair<int,double> > listDescCountInfo;	// <#events below, #ways below and here>
	for(map<int, set< int> > :: iterator itt = setNodesDestInfo.begin(); itt !=setNodesDestInfo.end(); ++itt)
	{
		// 
		YW_ASSERT_INFO(itt->second.size() >= 1, "Must have at least one");
//cout << "node = " << itt->first << ", itt->second: ";
//DumpIntSet(itt->second);
		// multiple lineages, so let us figure out how many events at that un-finished cluster
		int paridProc = itt->first;
		int numRemLins = itt->second.size();
//cout << "paridProc = " << paridProc << ", numRemLins = " << numRemLins << endl;
		int numDestLinsHere = 0;
		if( mapDestLinsInfo.find(paridProc) != mapDestLinsInfo.end() )
		{
			numDestLinsHere = mapDestLinsInfo[paridProc].size();
//cout << "Add size: " << mapDestLinsInfo[paridProc].size() << endl;
		}

		if( itt->second.size() == 1 && (mapDoneLinsParInfo.find(paridProc) == mapDoneLinsParInfo.end() ) )
		{
			// single ones, just record it
			pair<int,double> pp;
			int linDesc = *(itt->second.begin() );
			YW_ASSERT_INFO( mapDoneLinsCountInfo.find( linDesc ) != mapDoneLinsCountInfo.end(), "Fail to find the lineage"  );
			// in this case we just have the # of events as it has a single lineage
			//int numTotalSrcLins = numDestLinsHere + ;
			pp.first =  mapDoneLinsCountInfo[linDesc].numEvtsBelow;
			pp.second = mapDoneLinsCountInfo[linDesc].numWays;
			// it is possible several linages may trace to it, so check it
			listDescCountInfo.push_back( pp );
		}
		else
		{
			// if this par is not done yet, we do nothing since that node is not reached yet
			//if( mapDoneLinsCountInfo.find( paridProc ) != mapDoneLinsCountInfo.end() )
			//{
				//YW_ASSERT_INFO( (int)mapDoneLinsCountInfo[paridProc].lineagesDesc.size() >= numRemLins, "Not enough lineages1" );
				vector< pair<int,double> > listCCInfo;
				//else if( paridProc <0 )
				//{
				//	// this is the root lineage, so add one
				//	numDestLinsHere ++;
				//}
				int numEvtSubtot =  numDestLinsHere - numRemLins; 
				if( mapDoneLinsCountInfo.find( paridProc ) != mapDoneLinsCountInfo.end() )
				{
					numEvtSubtot+= (int)mapDoneLinsCountInfo[paridProc].lineagesDesc.size();
//cout << "done lineagesDesc.size(): " << mapDoneLinsCountInfo[paridProc].lineagesDesc.size() << endl;
				}
				if( mapDoneLinsCountInfo.find(paridProc) == mapDoneLinsCountInfo.end() && 
					mapDoneLinsParInfo.find( paridProc ) != mapDoneLinsParInfo.end() )
				{
					// if the node is not done (otherwise it is already included) and there are unfinished lineages
					numEvtSubtot+= (int)mapDoneLinsParInfo[paridProc].size();
//cout << "done mapDoneLinsParInfo.size(): " << mapDoneLinsParInfo[paridProc].size() << endl;
				}

				YW_ASSERT_INFO( numEvtSubtot >= 0, "Not enough lineages1" );
				vector<int> listDescToProcHere;
				if( mapDoneLinsCountInfo.find( paridProc ) != mapDoneLinsCountInfo.end() )
				{
					listDescToProcHere = mapDoneLinsCountInfo[paridProc].lineagesDesc;
//cout << "lineagesDesc in map: ";
//DumpIntVec(mapDoneLinsCountInfo[paridProc].lineagesDesc);
				//YW_ASSERT_INFO( mapDoneLinsCountInfo.find( paridProc ) != mapDoneLinsCountInfo.end(), "Fail to find the parent lineage"  );
				}
				if( numDestLinsHere > 0 )
				{
					for( set<int> :: iterator itjjj = mapDestLinsInfo[paridProc].begin(); itjjj != mapDestLinsInfo[paridProc].end(); ++itjjj  )
					{
						int idDesc = *itjjj;
						listDescToProcHere.push_back( idDesc );
					}
				}
				if( mapDoneLinsCountInfo.find(paridProc) == mapDoneLinsCountInfo.end() && 
					mapDoneLinsParInfo.find( paridProc ) != mapDoneLinsParInfo.end()  )
				{
					for( set<int> :: iterator itjjj = mapDoneLinsParInfo[paridProc].begin(); itjjj != mapDoneLinsParInfo[paridProc].end(); ++itjjj  )
					{
						int idDesc = *itjjj;
						listDescToProcHere.push_back( idDesc );
					}
				}

				for(int jjj=0; jjj<(int)listDescToProcHere.size(); ++jjj)
				//for(int jjj=0; jjj<(int)mapDoneLinsCountInfo[paridProc].lineagesDesc.size(); ++jjj)
				{
					int idDesc = listDescToProcHere[jjj];
					//int idDesc = mapDoneLinsCountInfo[paridProc].lineagesDesc[jjj];
					YW_ASSERT_INFO(mapDoneLinsCountInfo.find( idDesc) != mapDoneLinsCountInfo.end(), "Fail to find idDesc");
					pair<int,double> ppjj;
					ppjj.first = mapDoneLinsCountInfo[idDesc].numEvtsBelow;
					ppjj.second = mapDoneLinsCountInfo[idDesc].numWays;
					listCCInfo.push_back( ppjj);

					numEvtSubtot += mapDoneLinsCountInfo[idDesc].numEvtsBelow;
				}
				//if( numDestLinsHere > 0 )
				//{
				//	for( set<int> :: itjjj = mapDestLinsInfo[paridProc].begin(); itjjj != mapDestLinsInfo[paridProc].end(); ++itjjj  )
				//	{
				//		int idDesc = *itjjj;
				//		YW_ASSERT_INFO(mapDoneLinsCountInfo.find( idDesc) != mapDoneLinsCountInfo.end(), "Fail to find idDesc");
				//		pair<int,int> ppjj;
				//		ppjj.first = mapDoneLinsCountInfo[idDesc].numEvtsBelow;
				//		ppjj.second = mapDoneLinsCountInfo[idDesc].numWays;
				//		listCCInfo.push_back( ppjj);
				//		numEvtSubtot += mapDoneLinsCountInfo[idDesc].numEvtsBelow;
				//	}
				//}
				double numWaysNew = CountParFromChildCoal( listCCInfo, numRemLins );
				pair<int,double> pp2;
				pp2.first = numEvtSubtot;
				pp2.second = numWaysNew;
				listDescCountInfo.push_back( pp2 );
			}
		//}
	}
	// finally combine into the total number
	//int res = 1;
	//int numEvtsStep = 0;
	//for(int jj=0; jj<(int)listDescCountInfo.size(); ++jj)
	//{
	//	numEvtsStep += listDescCountInfo[jj].first;
	//	int nFactor = CalcNumNChooseK( numEvtsStep, listDescCountInfo[jj].first );
	//	res *= listDescCountInfo[jj].second*nFacCurr;
	//}
//cout << "listDescCountInfo size: " << listDescCountInfo.size() << endl;
	double res = CalcMultiBranchComboNum(listDescCountInfo );
//cout << "***CalcwbValMulti: final result is " << res << endl;
	const int MAX_WBVAL_CACHE_SIZE = 1000000;
	if( (int)cacheWbValMulti.size() < MAX_WBVAL_CACHE_SIZE )
	{
		cacheWbValMulti.insert( map< pair<multiset<int>, multiset<int> >, double> :: value_type( ppSrcDest, res) );
	}

	return res;
}

double GeneSpeciesTreeHelper :: CountParFromChildCoal( const vector< pair<int,double> > &listChildCounts, int numRemLins  ) const
{
    if( fVerboseMode == true )
    {
        cout << "CountParFromChildCoal, numRemLins = " << numRemLins  << ", listChildCounts: " << endl;
        for(int i=0; i<(int) listChildCounts.size();++i )
        {
            cout << "[" << listChildCounts[i].first << ", " << listChildCounts[i].second << ")\n";
        }
    }
	// count the number of ways of arranging these lineages under a single parent (like mulfurcating
	// tree node) and coalesce until numRemLins remaining lineages are left
	// use DP, starting pairwise and then work out all possible subsets
	//if( (int)listChildCounts.size() < numRemLins )
	//{
	//	// degenerate case: a single root, then just 1 choice
	//	return 1;
	//}

	YW_ASSERT_INFO( (int)listChildCounts.size() >= numRemLins, "Not enough lineages" );
	int szTable = 0x1 << listChildCounts.size();
//cout << "szTable = " << szTable << endl;
	// DP entry: <#evts below, #ways>
	vector<  pair<int, double> > tableDP(szTable);
	for(int i=0; i<szTable;++i)
	{
		tableDP[i].first = 0;
		tableDP[i].second = 0;
	}
	int numLinsTot = listChildCounts.size();
	for(int i=0; i<numLinsTot; ++i)
	{
		set<int> vecSingle;
		vecSingle.insert( i );
		int posSingle = ConvIntSetToPosition( vecSingle );
//cout << "posSingle = " << posSingle << endl;
		// conv to integer
		tableDP[posSingle] = listChildCounts[i];
	}
	// now try all possible subsets
	int numLinsTotExp = numLinsTot - numRemLins +1;
	for(int szSS = 2; szSS<=numLinsTotExp; ++szSS)
	{
//cout << "Processing szSS = " << szSS << endl;
		// consider all subset of size szSS
		vector<int> posvec;
		GetFirstCombo(szSS, numLinsTot,  posvec);
		while(true)
		{
			// treat each lineage to be the LAST one to coalesce and add up
			int nevts = posvec.size()-1;
			for(int iii=0; iii<(int)posvec.size(); ++iii)
			{
				nevts += listChildCounts[ posvec[iii] ].first;
			}
			set<int> setvecpos;
			PopulateSetByVec( setvecpos, posvec);
			double nways = 0.0;
			if( szSS == 2)
			{
				// just add it up
				// note event num is added by 1 since it consider the one on that node
				//nevts = listChildCounts[ posvec[0] ].first + listChildCounts[ posvec[1] ].first+1;
				nways = listChildCounts[ posvec[0] ].second * listChildCounts[ posvec[1] ].second
					* CalcNumNChooseK(listChildCounts[ posvec[0] ].first + listChildCounts[ posvec[1] ].first, 
					listChildCounts[ posvec[0] ].first);
			}
			else
			{
				// consider all ways of partitioning this set
				// to remove redundency, treat the last item as always in the first part
				for(int ttt =0; ttt<szSS-1; ++ttt)
				{
					// treat this lineage as the last one to coalesce
					// create a new vector
					vector<int> posvec2;
					GetFirstCombo(ttt, szSS-1,  posvec2);
					while(true)
					{
						set<int> setvecNew;
						setvecNew.insert( posvec[szSS-1]  );
						for( int iijj=0; iijj<(int)posvec2.size(); ++iijj  )
						{
							setvecNew.insert( posvec[ posvec2[iijj] ] );
						}
//cout << "ttt = " << ttt << ", setvecNew = ";
//DumpIntSet( setvecNew);
						set<int> setvecNewComp = setvecpos;
						SubtractSets(setvecNewComp, setvecNew);
						YW_ASSERT_INFO( setvecNewComp.size()+setvecNew.size() == setvecpos.size(), "Fail");

						if( setvecNew.size() > 0 && setvecNewComp.size() > 0 )
						{
							//setvecNew.erase( posvec[ttt] );
							int posChild1 = ConvIntSetToPosition( setvecNew );
							int posChild2 = ConvIntSetToPosition( setvecNewComp );
							//YW_ASSERT_INFO( tableDP.find(posChild1) != tableDP.end() && tableDP.find(posChild2) != tableDP.end(), 
							//	"positions not here" );
							nways += tableDP[ posChild1 ].second * tableDP[ posChild2 ].second
								* CalcNumNChooseK(tableDP[ posChild1 ].first + tableDP[ posChild2 ].first, tableDP[ posChild1 ].first);
//cout << "posChild1 = " << posChild1 << ", posChild2=" << posChild2 << ", nways now is " << nways << endl;
						}
						if( GetNextCombo(ttt, szSS-1,  posvec2) == false )
						{
							break;
						}
					}
				}
			}
			// assign it
			int posConv = ConvIntSetToPosition( setvecpos );
			tableDP[posConv].first = nevts;
			tableDP[posConv].second = nways;
//cout << "evtnum = " << nevts << ", nways = " << nways << ", for subsets ";
//DumpIntSet( setvecpos );
			if( GetNextCombo(szSS, numLinsTot,  posvec) == false )
			{
				break;
			}
		}
	}
	// now get the final results by combining numbers in the table
	// enumerate all partition of the lineages into nRemLins
	//vector<int> listChildLins;
	//for(int i=0; i<(int)listChildCounts.size(); ++i)
	//{
	//	listChildLins.push_back(i);
	//}
//cout << "Now compute the result of CountParFromChildCoal\n";
	// for now, just return the last cell
	// NEED TO CHANGE LATER
	//return tableDP[ szTable-1 ].second;

//#if 0
	//unsigned psize = numRemLins;
	//unsigned totSize = listChildCounts.size();
	//partition::iterator_k it(totSize, psize);
//cout << "here\n";
	// can not believe it? it is so hard but finally we are almost there...
	vector< vector<int> > parts;
	InitSubsetPartitionEnum(listChildCounts.size(), numRemLins, parts);

	double resFinal = 0;
	while(true)
	{
        //std::auto_ptr<std::vector<std::vector<int> > >
        //  part = it[listChildLins];
		//vector<unsigned> vecElemAssigned = *it;
// dump out
//cout << "Enumerate parts: \n";
//for(int i=0; i<(int)parts.size(); ++i)
//{
//DumpIntVec( parts[i] );
//}
		// partitions
		vector<set<int> > listPartitions;
		for(int i=0; i<(int)parts.size(); ++i)
		{
			YW_ASSERT_INFO( parts[i].size() > 0, "Can not be empty1" );
			//YW_ASSERT_INFO( parts[i] < psize, "Fail in enum" );
			//listPartitions[vecElemAssigned[i] ].insert(i);
			set<int> tset;
			PopulateSetByVec( tset, parts[i]);
			listPartitions.push_back( tset );
		}
		// assemble the list
		vector< pair<int,double> > listCountsForRes;
		for( int jjj=0; jjj<(int)listPartitions.size();++jjj )
		{
			YW_ASSERT_INFO( (int)listPartitions[jjj].size() <= numLinsTotExp, "Over the limit" );
			int ssId = ConvIntSetToPosition( listPartitions[jjj] );
			pair<int,double> ppjjj = tableDP[ssId];
			listCountsForRes.push_back( ppjjj );
		}
		resFinal += CalcMultiBranchComboNum( listCountsForRes );

		//
		if( GetNextSubsetPartitionEnum(listChildCounts.size(), numRemLins, parts) == false)
		{
			break;
		}
	}
    if( fVerboseMode == true )
    {
        cout << "number of counted value: " << resFinal << endl;
    }
	return resFinal;
//#endif
}


bool GeneSpeciesTreeHelper :: IsNodeAncestralToGT(int anc, int des) const
{
	// is anc ancestral to des (i.e. there is a path from anc to des)?
	YW_ASSERT_INFO( anc <(int)vecNodeAncestralGT.size() && des <(int)vecNodeAncestralGT.size(), "Out of range" );
	return vecNodeAncestralGT[anc][des];
}

int GeneSpeciesTreeHelper :: GetRootIdGT() 
{
	// find the one without parent
	for( map<int, TreeNode*> :: iterator it = mapIdNode.begin(); it != mapIdNode.end(); ++it)
	{
		if( it->second->GetParent() == NULL)
		{
			return it->first;
		}
	}
	YW_ASSERT_INFO(false, "Fail to find root");
	return -1;	// not found, bad
}

double GeneSpeciesTreeHelper :: CalcCoalProbBranch( int u, int v, double len )
{
	// have we computed it already?
	pair<int, int> puv(u,v);
	pair< pair<int,int>, double> ppp(puv,len);
	if( cacheBranchProb.find( ppp) != cacheBranchProb.end()  )
	{
//cout << "Cache hit\n";
		// done
		return cacheBranchProb[ppp];
	}

	double res = UtilsCalcCoalProbBranch(u, v, len);


	// save to cache
	cacheBranchProb.insert( map< pair<pair<int,int>,double>, double > :: value_type(ppp, res) ) ; 

//cout << ": prob = " << res << endl;
	return res;
}

void GeneSpeciesTreeHelper :: AddToCache( const LineageConfig &lin, const set<LineageConfig> &setAncestralCfgs)
{
	//
	const int MAX_CACHE_SIZE = 10000;
	if( (int)cacheLCAncs.size() > MAX_CACHE_SIZE )
	{
//cout << "Clear out cache...\n";
		cacheLCAncs.clear();
	}
	// add to it
	cacheLCAncs.insert( map<LineageConfig, set< LineageConfig > > :: value_type(lin, setAncestralCfgs) );
}


bool GeneSpeciesTreeHelper :: AreSrcDestCfgsCompatible( const LineageConfig &linCfgSrc, const LineageConfig &linCfgDest )
{
	// if src's size is smaller than dest, done 
	if( linCfgSrc.GetTotNumLins() < linCfgDest.GetTotNumLins() )
	{
		// 
		return false;
	}
	// compatible if each lineage in dest is ancestral to some lin in orig
	return linCfgDest.IsAncestralTo( linCfgSrc, *this );
}


// new functions: hopefully faster
void GeneSpeciesTreeHelper :: FastCollectAncConfig( int nodeId, const set<LineageConfig,LTLinCfg> &setOrigCfgs, vector< AncLinConfig *> &listAncCfgsNew )
{
	double lenOrig = treeSpecies.GetEdgeLen(nodeId);
	FastCollectAncConfig(lenOrig, setOrigCfgs, listAncCfgsNew );
}

void GeneSpeciesTreeHelper :: FastCollectAncConfig( double lenOrig, const set<LineageConfig,LTLinCfg> &setOrigCfgs, vector< AncLinConfig *> &listAncCfgsNew )
{
//#if 0
	//double lenOrig = treeSpecies.GetEdgeLen(nodeId);

//YW_ASSERT_INFO(false, "should not come here");
	// exam each orig config and find each way of extending
	// create a list if it is new
	map<string, AncLinConfig*> mapRepCfgs;		// stores representative configs
	for( set<LineageConfig,LTLinCfg> ::iterator it = setOrigCfgs.begin(); it != setOrigCfgs.end(); ++it )
	{
		// try all possible ways of coalescent
		FastCollACSingleCfg(lenOrig, *it, listAncCfgsNew, mapRepCfgs );
	}
//#endif

}

void GeneSpeciesTreeHelper :: FastCollACSingleCfg( double lenBranch, const LineageConfig &origCfgs, vector< AncLinConfig *> &listAncCfgsNew, 
												  map<string, AncLinConfig*> &mapExistACs )
{
//cout << "origCfgs: ";
//origCfgs.Dump();
	// loop over all possible iterations
	stack< pair<int,AncLinConfig *> > stackCfgToProc;
	//stack<AncLinConfig *> stackCfgToProc;
	// init AC: no events at all
	vector<CoalEvent> evtsEmpty;
	AncLinConfig *pacNew = new AncLinConfig(&origCfgs, origCfgs, evtsEmpty);
	pair<int, AncLinConfig *> pp(1, pacNew);
	stackCfgToProc.push( pp );
	//stackCfgToProc.push( pp );

	// local cache to ensure no AC is found duplicate
	set<LineageConfig> cacheLocalExistACs ;

	// this cache is local: it stores the already processed ancestral cfgs from the current cfg
	// the purpose is to avoid over-counting since CalcMultiVal already considers the situation
	// of multiple ways of coalescing. So we simply disallow symetric config here
	// but it seems there are still some small problems....
	//set<string> setStepFoundCfgStrs;
	set< multiset<int> > setStepFoundCfgParIds;

	while(stackCfgToProc.empty() == false )
	{
		YW_ASSERT_INFO( stackCfgToProc.top().second != NULL, "We do not use NULL in referring to AC" );
		AncLinConfig *pcfgCur = stackCfgToProc.top().second;
		int numRepeat = stackCfgToProc.top().first;
		stackCfgToProc.pop();
		// deal with probabilities...
		pcfgCur->UpdateProb( lenBranch, this );
		if( numRepeat > 1)
		{
			pcfgCur->SetProb( pcfgCur->GetProb() * numRepeat );
		}
//cout << "FastCollACSingleCfg: pcfgCur: ";
//pcfgCur->Dump();
		// is this new?
		//set<int> setId1;
		//pcfgCur->GetLinCfg().GetSrcNodes(setId1);
		//string sLabel1 = FindGNShapeLabel(setId1);
		string sLabel1 = FindGNShapeLabel(pcfgCur->GetLinCfg() );

		// also get par id
		multiset<int> msParIds;
		GetLinCfgParLinsSet(pcfgCur->GetLinCfg(), msParIds);

		//static int uniqeNum = 1;
		//char buf[100];
		//sprintf( buf, "%d", uniqeNum++);
		//sLabel1 += buf;
//cout << "sLabel1 = " << sLabel1 << endl;
//cout << "mapExistACs: ";
//for( map<string,AncLinConfig*> :: iterator itttt = mapExistACs.begin(); itttt != mapExistACs.end(); ++itttt  )
//{
//cout << itttt->first << ", ";
//}
//cout << endl;
		bool fDeleteAfter = false;

//#if 0		// for now, do not consolidate states
		// YW: disable config merging 
		if( mapExistACs.find( sLabel1) != mapExistACs.end() )
		//if( false)
		{
			// 2/18/11: do nothing if mulfurcating trees since it is already included
			if( IsGeneTreeMulfurcate() == false  || setStepFoundCfgParIds.find(msParIds) == setStepFoundCfgParIds.end() )
			// if finding exact the same entry, append it since this will be due to different configs
			// and not due to multiplicity of non-binary nodes (4/6/11)
			//if( cacheLocalExistACs.find(pcfgCur->GetLinCfg() ) != cacheLocalExistACs.end() )
			{
				// not new, so we just add the prob to it
				AncLinConfig *pOrigAC = mapExistACs[sLabel1];
				YW_ASSERT_INFO( pOrigAC != NULL, "Can not be NULL" );
//cout << "not new\n";
//cout << "Before appending, orig AC become: ";
//pOrigAC->Dump();
				pOrigAC->AppendAC( *pcfgCur );
//cout << "After appending, AC become: ";
//pOrigAC->Dump();
			}
			// freeup pcfgcur since it is done
			// YW: this seems to cause memory leak: TBD
			fDeleteAfter = true;
			//delete pcfgCur;
//cout << "Done with delete\n";
		}
		else
//#endif
		{
			// new, so add to stack and act as repretattive
//cout << "new\n";
			// and create a new record for this
			// add to stack
			listAncCfgsNew.push_back( pcfgCur );
if( pcfgCur->GetProb() < 0.0)
{
cout << "pcfgCur is: ";
pcfgCur->Dump();
YW_ASSERT_INFO(false, "Ancestral config has negative prob");
}
			//AncLinConfig *pacNew = &listAncCfgsNew[ listAncCfgsNew.size()-1 ];
			//AncLinConfig *pacNew = &cfgCur;
			//stackCfgToProc.push(cfgCur);
			mapExistACs.insert( map<string, AncLinConfig*> :: value_type( sLabel1, pcfgCur ) ) ;
		}

		//if( fDeleteAfter == false ||  IsGeneTreeMulfurcate() == false )
		//{
		
		// if new, we then continue to explore its descendent
		// what ACs reachable from the current cfg. If duplicate is found, simply accumulate prob here
		set< set<LineageCluster> > setCoalLins;
		pcfgCur->GetLinCfg().FindCoalescableLins( setCoalLins );
//cout << "setCoalLins size: " << setCoalLins.size() << endl;
//cout << "setCoalLins: \n";
//for(set< set<LineageCluster> >:: iterator ittg1= setCoalLins.begin(); ittg1 !=setCoalLins.end();++ittg1 )
//{
//cout << "###DUMP one set of lineages to coalesce\n";
//for( set<LineageCluster> :: iterator ittg2 = ittg1->begin(); ittg2 !=ittg1->end(); ++ittg2)
//{
//ittg2->Dump();
//}
//}
		// add record
		setStepFoundCfgParIds.insert(msParIds);


	//cout << "Work on coal lins...\n";
		// for the coalesable sequences (two or more)
		for( set< set<LineageCluster> > :: iterator it = setCoalLins.begin(); it != setCoalLins.end(); ++it )
		{
//cout << "cluster size = " << it->size() << endl;
			// 
			if( it->size() <= 1 )
			{
				continue;
			}
			// now consider each cluster with two or more sequences and merge them pairwise
			LineageConfig lcNew;
			set<LineageCluster> :: iterator it2 = it->begin();
			YW_ASSERT_INFO( it2 != it->end(), "Fail" );
			set<LineageCluster> :: iterator it3 = it2;
			it3++;
			YW_ASSERT_INFO( it3 != it->end(), "Fail" );
			//if( it->size() == 2 )
			//{
			// in this case, remove and combine
			pcfgCur->GetLinCfg().MergeLCInto(*this, *it2, *it3, lcNew );
			//}
			//else
			//{
			//	YW_ASSERT_INFO( it->size() > 2, "Fail: must have three or more lineages" );
			//	// otherwise we would simply remove one lineage, say *it2
			//	pcfgCur->GetLinCfg().RemoveLineage( *this, *it2, lcNew );
			//}
			// is this new? If never seen before inside this function call, do it
			if( cacheLocalExistACs.find( lcNew ) != cacheLocalExistACs.end() )
			{
//YW_ASSERT_INFO(false, "Should not be here");
				continue;
			}
			cacheLocalExistACs.insert( lcNew );

			// create a new AC
//#if 0
			vector<CoalEvent> evtsCoal;
			pcfgCur->GetCoalEvts(evtsCoal);
			AddCoalEvt( evtsCoal, *it2, *it3 );
			//pacNew->SetCoalEvts( evtsCoal );
			AncLinConfig *pacNew = new AncLinConfig(&origCfgs, lcNew, evtsCoal);
//cout << "Creating a new Anc Lin: ";
//pacNew->Dump();
			// YW: 10/28/10, set copy to be 1 and handle mulfurcation elsewhere
			//int numDups = (it->size() -1)*(it->size())/2;
			int numDups = 1;
			YW_ASSERT_INFO( numDups >= 1, "Invalid numDups" );
			pair<int, AncLinConfig *> pp(numDups, pacNew);
			stackCfgToProc.push(pp);
//#endif

		}
		//}

		// delete now
		if( fDeleteAfter == true )
		{
			delete pcfgCur;
			pcfgCur = NULL;
		}
		//}

		// if too many ancestral AC, stop
		if( (int)listAncCfgsNew.size() > maxLinAncCfgNum )
		{
			break;
		}
	}

//cout << "***After this config, list of ancestral configs so far***\n";
//for(  int i=0; i<(int)listAncCfgsNew.size(); ++i )
//{
//listAncCfgsNew[i]->Dump();;
//}
//cout << "Done: FastCollACSingleCfg\n";
}



void GeneSpeciesTreeHelper :: AddCoalEvt( vector<CoalEvent> &evtsCoal, const LineageCluster &lin1, const LineageCluster &lin2 ) const
{
	// for now, get MRCA of two lins and create an event
	// CAUTION: does not quite work for non-binary tree. Need change!!!!
	int nid1 = lin1.GetGNId();
	int nid2 = lin2.GetGNId();
	int ntid = GetMRCATwoLineagesGT( nid1, nid2 );
	CoalEvent cevt( ntid );
	evtsCoal.push_back( cevt );
}



string GeneSpeciesTreeHelper :: FindGNShapeLabel( const LineageConfig &cfgCurr )
{
//#if 0
	// use the label of nodes (not internal id) so as to merge nodes
	set<int> setId1;
	cfgCurr.GetSrcNodes(setId1);

	// now setup
    // 08/12/15
	string sLabel = treeGene.GetShapeLabel(setId1, false );
    string resLbl = GetDistinctString(sLabel);
//cout << "sLabel: " << sLabel << ", resLbl: " << resLbl << endl;
//cout << " for configuraiton: ";
//cfgCurr.Dump();
	return resLbl;
//#endif

}


int GeneSpeciesTreeHelper :: GetSpeciesForGTLin(int lin) const
{
	// for a GT lineage, what species it comes from
	YW_ASSERT_INFO(mapGTLinToSpecies.find( lin ) != mapGTLinToSpecies.end(), "This lineage has no information" );
	return const_cast<GeneSpeciesTreeHelper *>(this)->mapGTLinToSpecies[lin];
}

//int GeneSpeciesTreeHelper :: GetDescendentsNumGT(int nodeGT) const
//{
//	return GetGTNodeNumChildren( nodeGT );
//	if(mapNodesDescendentNumGT.find( nodeGT ) == mapNodesDescendentNumGT.end() )
//	{
//		return 0;
//	}
//	return const_cast<GeneSpeciesTreeHelper *>(this)->mapNodesDescendentNumGT[nodeGT];
//}

void GeneSpeciesTreeHelper :: GetGTNodeDescendents(int nodeGT, set<int> &listDescs) const
{
	// get descendents
	listDescs.clear();
	// ???
	if( mapNodesDescendentsGT.find( nodeGT ) != mapNodesDescendentsGT.end()  )
	{
		listDescs = const_cast<GeneSpeciesTreeHelper *>(this)->mapNodesDescendentsGT[nodeGT];
	}
}

void GeneSpeciesTreeHelper :: InitDefaultGTLinMap()
{
	mapGTLinToSpecies.clear();
	// by default, lineage i in GT belongs to species i
	for(int i=0; i<treeGene.GetNumLeaves(); ++i)
	{
		mapGTLinToSpecies.insert( map<int,int> :: value_type(i,i) );
	}
}

void GeneSpeciesTreeHelper :: SetupSpeciesIdForGTLin( const vector< pair<int,int> > &listRecords)
{
	// the list contains <gene tree id, species id>
	mapGTLinToSpecies.clear();
	for( int i=0; i<(int)listRecords.size(); ++i  )
	{
		mapGTLinToSpecies.insert( map<int,int> :: value_type(listRecords[i].first, listRecords[i].second) );
	}
}

void GeneSpeciesTreeHelper :: FindNumLinsForEvtsCfg( const vector<int> &listEvts, map<int,int> &mapEvtDupNums, 
													const LineageConfig &linCfgSrc, vector<int> &listLinsNumForEvt) const
{
#if 0
cout << "listEvts: ";
DumpIntVec( listEvts );
cout << "mapEvtDupNums: \n";
for(map<int,int>:: iterator it=mapEvtDupNums.begin(); it!=mapEvtDupNums.end();++it)
{
cout << "[" << it->first <<", " << it->second << "]\n";
}
cout << "linCfgSrc = ";
linCfgSrc.Dump();
#endif
	// for each evt for a given config, the number of lineages that may contributes to the event
	// in binary tree, it is always 2. But for mulfurcating tree, this may vary: depends on how many lineages are involved
	//YW_ASSERT_INFO( listEvts.size() == listEvtDupNums.size(), "Must match size" );
	listLinsNumForEvt.clear();
	//listLinsNumForEvt.resize(listEvts.size());
	map<int,int> mapEvtLinNums;
	for(int i=0; i<(int)listEvts.size();++i)
	{
		mapEvtLinNums.insert(map<int,int> :: value_type(listEvts[i], 0) );
		//listLinsNumForEvt[i] = 0;
	}
	// we loop through each lineage of cfg and update the value
	// be CAREFUL: internal lineages should also be treated
	set< LineageCluster > setLCsRes;
	linCfgSrc.GetLCs( setLCsRes );
	for(set< LineageCluster > :: iterator it = setLCsRes.begin(); it!= setLCsRes.end(); ++it)
	{
		// 
		int idPar = it->GetGNParId();
		if( mapEvtLinNums.find( idPar ) != mapEvtLinNums.end() )
		{
			mapEvtLinNums[idPar] ++;
		}
	}
#if 0
cout << "mapEvtLinNums: \n";
for(map<int,int>:: iterator it=mapEvtLinNums.begin(); it!=mapEvtLinNums.end();++it)
{
cout << "[" << it->first <<", " << it->second << "]\n";
}
#endif
    
	// the number for internal number is decided by numdupevts
	for(int i=0; i<(int)listEvts.size(); ++i)
	{
		int idevt = listEvts[i];		// also the id of the node in GT
        
#if 0
if(mapEvtLinNums.find( idevt ) == mapEvtLinNums.end())
{
cout << "idevt: " << idevt << endl;
}
if( mapEvtDupNums[idevt]+1 > mapEvtLinNums[idevt])
{
cout << "idevt2: " << idevt << endl;
}
#endif
		YW_ASSERT_INFO( mapEvtDupNums.find( idevt ) != mapEvtDupNums.end(), "Evt multicitlies: not specified");
		YW_ASSERT_INFO( mapEvtLinNums.find( idevt ) != mapEvtLinNums.end(), "Must has it");
		YW_ASSERT_INFO( mapEvtDupNums[idevt]+1 <= mapEvtLinNums[idevt] , "Must has it2");
		bool fDoneInLC = false;
		if( mapEvtDupNums[idevt]+1 == mapEvtLinNums[idevt]  )
		{
			fDoneInLC = true;
		}
		//if( GetNumLinsAtNode( idevt ) == mapEvtDupNums[idevt]+1  )
		if( fDoneInLC )
		{
			// this lineage is done
			YW_ASSERT_INFO( mapNodeParGT.find(idevt) != mapNodeParGT.end(), "Fail to find gtree node in par map" );
			int ppid = const_cast<GeneSpeciesTreeHelper *>(this)->mapNodeParGT[idevt];
			if( ppid >= 0  && mapEvtLinNums.find( ppid) != mapEvtLinNums.end() )
			{
				// 
				mapEvtLinNums[ppid]++;
			}
		}
	}
//cout << "AFTER: mapEvtLinNums: \n";
//for(map<int,int>:: iterator it=mapEvtLinNums.begin(); it!=mapEvtLinNums.end();++it)
//{
//cout << "[" << it->first <<", " << it->second << "]\n";
//}

	// finally set up a value
	listLinsNumForEvt.resize(listEvts.size());
	for(int i=0; i<(int)listEvts.size();++i)
	{
		int numLins = mapEvtLinNums[listEvts[i] ];
		YW_ASSERT_INFO( numLins >= 2, "Must have at least two lineages" );
		listLinsNumForEvt[i] = numLins;
	}
}

int GeneSpeciesTreeHelper :: GetNumLinsAtNode(int gnid) const
{
	YW_ASSERT_INFO( mapIdNode.find( gnid) != mapIdNode.end(), "Fail to find the gene tree node" );
	return const_cast<GeneSpeciesTreeHelper *>(this)->mapIdNode[gnid]->GetChildrenNum();
}

void GeneSpeciesTreeHelper :: GetNodesForIds(const set<int> &nids, set<TreeNode *> &setNodes) const
{
    //
    GeneSpeciesTreeHelper *pthis = const_cast<GeneSpeciesTreeHelper *>(this);
    setNodes.clear();
    for(set<int> :: iterator it = nids.begin(); it != nids.end(); ++it)
    {
        YW_ASSERT_INFO( mapIdNode.find(*it) != mapIdNode.end(), "Fail to find");
        TreeNode *pn = pthis->mapIdNode[ *it ];
        setNodes.insert(pn);
    }
}

void GeneSpeciesTreeHelper :: GetLinCfgParLinsSet(const LineageConfig &linCfgSrc, multiset<int> &sint) const
{
	// get the par ids set
	sint.clear();

	//
	linCfgSrc.GetParLins(sint);
}

////////////////////////////////////////////////////////////////////////////////////////////////


LinCfgSoreIterator :: LinCfgSoreIterator( int id, LineageConfigSore &ls ) : nodeId(id), lcStore(ls)
{
	// initialize it
	Init();
}

void LinCfgSoreIterator :: Init()
{
	YW_ASSERT_INFO( lcStore.statesDepot.find( nodeId) != lcStore.statesDepot.end(), "Fail to find it"  );
	it = lcStore.statesDepot[nodeId].begin();
}

bool LinCfgSoreIterator :: IsDone()
{
	return it == lcStore.statesDepot[nodeId].end();
}

void LinCfgSoreIterator :: Next()
{
	it ++;
}



////////////////////////////////////////////////////////////////////////////////////////////////

LineageConfigSore :: LineageConfigSore(MarginalTree &ts, PhylogenyTreeBasic &tg, GeneSpeciesTreeHelper &gst) : treeSpecies(ts), treeGene(tg), gstHelper(gst)
{
}

// add a record of LinCfg: for the purpose of computing the probablity, need to know which
// LinCfg it derives from. For leaves, we also support a plain version: with no probility
void LineageConfigSore :: AddLinCfg( int nodeId, LineageConfig &lcNew )
{
//#if 0
	// for leaves, set prob to be 1.0
	//lcNew.SetProb(1.0);
//cout << "AddLinCfg: node id = " << nodeId <<", cfg to add: ";
//lcNew.Dump();

	// in this case, this must be the only state
	//YW_ASSERT_INFO( statesDepot.find( nodeId ) == statesDepot.end(), "Fail: leaf must only have one LinCfg" );
	if(  statesDepot.find( nodeId ) == statesDepot.end()  )
	{
		set<LineageConfig,LTLinCfg> setLins;
		setLins.insert(lcNew);
		statesDepot.insert( map<int, set< LineageConfig,LTLinCfg > > :: value_type(nodeId, setLins)  );
	}
	else
	{
        // YW: 08/10/15. CHANGE: if the cfg is already in the store, append it
        if( statesDepot[nodeId].find(lcNew) != statesDepot[nodeId].end()  )
        {
            LineageConfig lcTot = *statesDepot[nodeId].find( lcNew );
            lcTot.SetProb( lcTot.GetProb()+lcNew.GetProb());
            statesDepot[nodeId].erase( lcNew );
            statesDepot[nodeId].insert( lcTot );
        }
        else
        {
            // Otherwise, just add it
            // can not be the same, already there
            //YW_ASSERT_INFO( statesDepot[nodeId].find( lcNew ) == statesDepot[nodeId].end(), "Can not insert duplicate lcNew"  );
            statesDepot[nodeId].insert( lcNew );
        }
        
		// now prune it
		Prune(nodeId);
	}
//#endif

}

void LineageConfigSore :: Prune(int nodeId)
{
	// in case there are too many cfgs, prune it
	if( (int)statesDepot[nodeId].size() >  maxLinCfgNum )
	{
		// prune half of it
		set< LineageConfig,LTLinCfg > setCfgsKept;

		// find the median of prob
		vector<double> listProbsCur;
		set< LineageConfig,LTLinCfg > &listLinCfgsMap = statesDepot[nodeId];
		for(  set< LineageConfig,LTLinCfg > :: iterator itt = listLinCfgsMap.begin(); itt != listLinCfgsMap.end(); ++itt )
		{
			listProbsCur.push_back( itt->GetProb() );
		}
		double medianProb = FindMedian(listProbsCur);
		// scan it again
		for(  set< LineageConfig,LTLinCfg > :: iterator itt = listLinCfgsMap.begin(); itt != listLinCfgsMap.end(); ++itt )
		{
			if( itt->GetProb() > medianProb)
			{
				setCfgsKept.insert( *itt );
			}
		}

		// update list
		statesDepot[nodeId] = setCfgsKept;
	}

}

#if 0
LineageConfig *LineageConfigSore :: FindLinCfgRep(int nodeId, const LineageConfig &linCfg) 
{
	// find a lin cfg that is represntative of this linCfg
	LineageConfig *pcfg = NULL;
	// check to see if any current stored items
	if(  statesDepot.find( nodeId ) != statesDepot.end()  )
	{
		int numLins = linCfg.GetTotNumLins();
		//set<int> setId1;
		//linCfg.GetSrcNodes(setId1);
		//string sLabel1 = treeGene.GetShapeLabel(setId1);
		//string sLabel1 = FindGNShapeLabel(setId1);
		string sLabel1 = gstHelper.FindGNShapeLabel(linCfg);
		set< LineageConfig,LTLinCfg > &listLinCfgsMap = statesDepot[nodeId];
		for(  set< LineageConfig,LTLinCfg > :: iterator itt = listLinCfgsMap.begin(); itt != listLinCfgsMap.end(); ++itt )
		{
			// has to have the same lin num
			if( numLins == itt->GetTotNumLins()  )
			{
				// are they the same
				//set<int> setId2;
				//itt->GetSrcNodes(setId2);
				//string sLabel2 = treeGene.GetShapeLabel(setId2);
				//string sLabel2 = FindGNShapeLabel(setId2);
				string sLabel2 = gstHelper.FindGNShapeLabel(*itt);
				if( sLabel1 == sLabel2 )
				{
					return const_cast<LineageConfig *>(&(*itt));
				}
			}
		}
	}

	return pcfg;

}
#endif



int LineageConfigSore :: GetNumLinCfgsAt(int nodeId)
{
	// 
	if( statesDepot.find( nodeId ) == statesDepot.end() )
	{
		return 0;
	}
	else
	{
		return statesDepot[nodeId].size();
	}
}

void LineageConfigSore :: DumpConfigsAt(int nodeId)
{
	if( statesDepot.find( nodeId ) == statesDepot.end() )
	{
		return;
	}
	cout << "Configs at node " << nodeId << endl;
	for( set< LineageConfig,LTLinCfg > :: iterator it = statesDepot[nodeId].begin(); it != statesDepot[nodeId].end(); ++it )
	{
		it->Dump();
	}
}

void LineageConfigSore :: SetLinCfgProbAt(int nodeId, LineageConfig &lcNew, double probNew)
{
	YW_ASSERT_INFO(statesDepot.find( nodeId ) != statesDepot.end(), "Can not find nodeId" );
	//YW_ASSERT_INFO(statesDepot[nodeId].find( lcNew ) != statesDepot[nodeId].end(), "Can not find lcNew at nodeId" );
	// YW: changed at 080510, to support pruning
	if( statesDepot[nodeId].find( lcNew ) != statesDepot[nodeId].end() )
	{
		set< LineageConfig,LTLinCfg >:: iterator it = statesDepot[nodeId].find(lcNew);
		const_cast< LineageConfig * >(&(*it))->SetProb( probNew );
	}
}

double LineageConfigSore :: CalcTotProbAt(int nodeId)
{
	// this should be the root
	// here each config is going to coalesce into the final (infinite) branch out of species tree root
	YW_ASSERT_INFO( statesDepot.find( nodeId ) != statesDepot.end(), "Bad root node2"  );
	double res = 0.0;
	for( set< LineageConfig,LTLinCfg > :: iterator it = statesDepot[nodeId].begin(); it != statesDepot[nodeId].end(); ++it  )
	{
		// If this config has more than one branch, then need to compute one additional coefficient
		if( it->GetTotNumLins() == 1 )
		{
			res += it->GetProb();
		}
		else
		{
			// need to compute the final combinatorial factor
			LineageConfig finalLC;
			LineageCluster flc;
			int rtGT = gstHelper.GetRootIdGT();
			flc.Init(rtGT, 1, -1);
			finalLC.AddLC( flc );
			double coeff = gstHelper.CalcCoalCoeffForTwoCfg(*it, finalLC);
//cout << "CalcTotProbAt: coeff = " << coeff << endl;
			res += it->GetProb()*coeff;
			if( coeff < 0.0  )
			{
				//cout << "coeff = " << coeff << endl;
				//YW_ASSERT_INFO(false, "Coefficient wrong");
                fSTELLSProbComputeWarningFlag = true;
                coeff = 0.0;    // set prob to be 0.0 for now
			}
			if( it->GetProb() < 0.0 )
			{
				//cout << "prob of the AC is wrong: ";
				//it->Dump();
				//YW_ASSERT_INFO(false, "AC prob wrong");
                fSTELLSProbComputeWarningFlag = true;
                LineageConfig *pcfgUse = const_cast<LineageConfig *>(&(*it));
                pcfgUse->SetProb(0.0);
			}
		}
	}
	return res;
}




////////////////////////////////////////////////////////////////////////////////////////////////

// main class for computing a given gene tree on species tree probability

// need two trees
GeneSpeciesTreeProb :: GeneSpeciesTreeProb(MarginalTree &ts, PhylogenyTreeBasic &tg) : treeSpecies(ts), treeGene(tg)
{
	treeSpecies.BuildDescendantInfo();
    
    // get the subtree that we will use
    set<int> species;
    GetLeavesSpeciesFromGT( treeGene, species );
    CreateSubtreeFromLeaves( treeSpecies, species, treeSpeciesUse, this->mapNewNodeToOldNode );
    
    pgstHelper = new GeneSpeciesTreeHelper(treeSpeciesUse, treeGene);
    pstoreLinCfgs = new LineageConfigSore(ts, tg, *pgstHelper);
}

GeneSpeciesTreeProb :: ~GeneSpeciesTreeProb()
{
    delete pgstHelper;
    delete pstoreLinCfgs;
}


double GeneSpeciesTreeProb :: CalcProb()
{
    if(fVerboseMode == true)
    {
        string strNW;
        treeGene.ConsNewick(strNW);
        cout << "CalcProb: species tree: " << treeSpecies.GetNewick() << ", gene tree: " << strNW << endl;
        cout << "gene tree: ";
        treeGene.Dump();
    }
    
#if 0
	if( gstHelper.IsGeneTreeMulfurcate() == true)
	{
		cout << "Caution: the gene tree is not binary...\n";
	}
#endif
	// first make sure starting from fresh
	Reset();

	// test drive the functions
	// we scan the nodes of species tree,
	//map<int, set< LineageConfig > > statesDepot;
	for(int iv = 0; iv < treeSpeciesUse.GetTotNodesNum(); ++iv)
	{
        if(fVerboseMode == true)
        {
            cout << "iv = " << iv << endl;
        }
		//set<LineageConfig> setLinsEmpty;
		//statesDepot.insert( map<int, set< LineageConfig > > :: value_type(iv, setLinsEmpty)  );

		// is it a leaf? If so, collect all the corresponding lineages for this species
		if( treeSpeciesUse.IsLeaf( iv ) == true )
		{
//cout << "Is leaf\n";
			// get label for it
			int lbl = treeSpeciesUse.GetLabel(iv );
			set<int> geneAlleles;
			pgstHelper->GetGeneAllelesForSpecies(lbl, geneAlleles);
//cout << "geneAlleles: ";
//DumpIntSet( geneAlleles );
			LineageConfig linConfig;
			pgstHelper->FormLinConfig( iv, geneAlleles, linConfig );
			linConfig.SetProb(1.0);		// leaves are prob 1.0

			// add a record for it
			//set<LineageConfig> setLins;
			//setLins.insert( linConfig );
			//statesDepot.insert( map<int, set< LineageConfig > > :: value_type(iv, setLins)  );
			//statesDepot[iv].insert( linConfig );
			pstoreLinCfgs->AddLinCfg(iv, linConfig);
//linConfig.Dump();
		}
		else
		{
//cout << "Non-leaf\n";
			//int numProc = 0;

			// otherwise, we examine its two descendents to see what kinds of states they have
			// we can delete them afterwards (TBD)
			// first find the two descendents of iv
			int nodeLeft = treeSpeciesUse.GetLeftDescendant(iv);
			int nodeRight = treeSpeciesUse.GetRightDescendant(iv);
            if( fVerboseMode == true )
            {
                cout << "LeftNode = " << nodeLeft << ", rightNode=" << nodeRight << endl;
                cout << "The number of config at left = " << pstoreLinCfgs->GetNumLinCfgsAt(nodeLeft) << endl;
                cout << "The number of config at right = " << pstoreLinCfgs->GetNumLinCfgsAt(nodeRight) << endl;
            }
			//YW_ASSERT_INFO( statesDepot.find( nodeLeft) != statesDepot.end(), "Wrong1" );
			//YW_ASSERT_INFO( statesDepot.find( nodeRight) != statesDepot.end(), "Wrong2" );
//cout << "Num of statesDepot in left = " << statesDepot[nodeLeft].size() << endl;
//cout << "Num of statesDepot in right = " << statesDepot[nodeRight].size() << endl;


			// collect ancestral cfgs on left/right
			vector<AncLinConfig *> configsNewLeft;
			set< LineageConfig,LTLinCfg > cfgListOldLeft = pstoreLinCfgs->GetLCSetAt(nodeLeft);
			pgstHelper->FastCollectAncConfig( nodeLeft, cfgListOldLeft, configsNewLeft  );
//cout << "configsNewLeft.size = " << configsNewLeft.size() << endl;
			vector<AncLinConfig *> configsNewRight;
			set< LineageConfig,LTLinCfg > cfgListOldRight = pstoreLinCfgs->GetLCSetAt(nodeRight);
			pgstHelper->FastCollectAncConfig( nodeRight, cfgListOldRight, configsNewRight  );
//cout << "configsNewRight.size = " << configsNewRight.size() << endl;
			// now let us merge them pairwise to form new cfg
			// append and insert into state depot
			for( int ii=0; ii< (int)configsNewLeft.size();  ++ii )
			{
//cout << "Here\n";
//cout << "Left ances: ";
//configsNewLeft[ii]->Dump();
				for(  int jj = 0; jj<  (int)configsNewRight.size(); ++jj )
				{
//cout << "Right ances: ";
//configsNewRight[jj]->Dump();
					// form a new one
					LineageConfig lcNew;
					pgstHelper->FormLinConfig(iv,configsNewLeft[ii]->GetLinCfg(),  configsNewRight[jj]->GetLinCfg(), lcNew);
					// setup prob
					lcNew.SetProb(  configsNewLeft[ii]->GetProb()*configsNewRight[jj]->GetProb() );
					if( lcNew.GetProb() < 0.0)
					{
						cout << "The left config prob is: " << configsNewLeft[ii]->GetProb() << ", and the right config prob is: " << configsNewRight[jj]->GetProb() << endl;
						YW_ASSERT_INFO(false, "AC: prob is negative");
					}
					//statesDepot[iv].insert( lcNew );
					pstoreLinCfgs->AddLinCfg(iv, lcNew );

					//numProc++;
					//if( (numProc % 10000 ) == 0 )
					//{
					//	cout << "Processing " << numProc << "...\n";
					//}
//cout << "After FormLinConfig, find a new LC: ";
//lcNew.Dump();
				}
			}
			// free mem
			for( int ii=0; ii< (int)configsNewLeft.size();  ++ii )
			{
				delete configsNewLeft[ii];
			}
			for(  int jj = 0; jj<  (int)configsNewRight.size(); ++jj )
			{
				delete configsNewRight[jj];
			}

//cout << "The total number of config at this node = " << storeLinCfgs.GetNumLinCfgsAt(iv) << endl;
		}
		if( fVerboseMode == true)
		{
			cout << "*************At species tree node " << iv << ", list of configurations: \n";
			pstoreLinCfgs->DumpConfigsAt(iv);
		}
	}

	int rootId = treeSpeciesUse.GetTotNodesNum()-1;
	double totProb = pstoreLinCfgs->CalcTotProbAt(rootId);
	//cout << "Total probability = " << setprecision(12)<< totProb << endl;

//#if 0
	if( fVerboseMode == true)
	{
		int maxNumCfgs = 0;
		int totNumCfgs = 0;
		for(int nid=0; nid<=rootId; ++nid)
		{
		int n1 = pstoreLinCfgs->GetNumLinCfgsAt( nid );
		totNumCfgs+= n1;
		if(maxNumCfgs < n1 )
		{
			maxNumCfgs = n1;
		}
		//cout << "At node " << nid << ", number of configurations: ";
		//cout << storeLinCfgs.GetNumLinCfgsAt( nid ) << endl;
		}
		cout << "Maximum number of config at this computation: " << maxNumCfgs << endl;
		cout << "TOTAL number of config at this computation: " << totNumCfgs << endl;
	}
//#endif
	return totProb;
}


double GeneSpeciesTreeProb :: TestNewBranchLenAt(int threadId, int branch, double lenNew, map<int,set< LineageConfig,LTLinCfg > > &origLinCfgs, bool fSetBrLen )
{
    // simply set the ebranch length
    double brOld = treeSpecies.GetEdgeLen( branch );
    
    if( fSetBrLen == true )
    {
        SetSTBranchLen(branch, lenNew);
    }
    double res = CalcLogProb();
    
    // set the old branch length back
    if( fSetBrLen == true )
    {
        SetSTBranchLen( branch, brOld );
    }
    
    
    return res;
}

void GeneSpeciesTreeProb :: SetLinCfgs(map<int,set< LineageConfig,LTLinCfg > > &linCfgsToSet)
{
//cout << "Set cfgs for nodes: ";
//for( map<int,set< LineageConfig,LTLinCfg > >::iterator it= linCfgsToSet.begin(); it !=linCfgsToSet.end(); ++it  )
//{
//cout << "node: " << it->first << endl;
//}
	// useful when we want to un-do the effect of the previous tweaking of branch length
	for( map<int,set< LineageConfig,LTLinCfg > > :: iterator itt = linCfgsToSet.begin(); itt != linCfgsToSet.end(); ++itt )
	{
		pstoreLinCfgs->SetLCSetAt(itt->first, itt->second);
	}
}

double GeneSpeciesTreeProb :: GetSTBrachLen(int br) const
{
    // use the one that we use
    return treeSpeciesUse.GetEdgeLen(br);
}

void GeneSpeciesTreeProb :: SetSTBranchLen(int br, double brLen)
{
	treeSpecies.SetBranchLen(br, brLen);
    
    // now should set the branch length of the used species tree
    UpdateBranchLenInSubtree( treeSpecies, this->mapNewNodeToOldNode, treeSpeciesUse );
    
	pgstHelper->ClearBranchLenCache();
}


