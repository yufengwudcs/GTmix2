//
//  ApproxGeneTreeProb2.cpp
//  
//
//  Created by Yufeng Wu on 9/26/24.
//
//

#include "ApproxGeneTreeProb2.h"
#include "UtilsNumerical.h"
#include "Utils4.h"
#include "pthread.h"
#include <cmath>
#include <queue>
#include <ctime>
#include <algorithm>
#include <thread>
#include <mutex>
//#include "PolynomialUtils.hpp"
//#include "FFTUtils.h"

//extern int GetNumThreads();

static bool fVerboseMode = false;
static bool fVerboseModeExtern = false;
static bool ffastSTELLSProbComputeWarningFlag = false;
//static int numBranchLenDigit = 6;
static int numBranchLenDigit = 4;
static int numBranchLenDigitPreCalc = 10001;

void AGTPSetVerbose(bool f)
{
    fVerboseModeExtern = f;
    //fVerboseMode = f;
}

#define LOG_CFG_PROB 1

// snap edge lengths
//static const int DEF_SIZE_BR_LEN_GRID_AGTP = 20;
static const double LIST_DEF_BR_LEN_GRID_VALS_AGTP[] = {0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 3.0, 5.0};
//static const double LIST_DEF_BR_LEN_GRID_VALS_AGTP[DEF_SIZE_BR_LEN_GRID_AGTP] = {0.0001, 0.0005, 0.001,0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.3, 0.5};
//static const double LIST_DEF_BR_LEN_GRID_VALS_AGTP[] = {0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
//static const double LIST_DEF_BR_LEN_GRID_VALS_AGTP[] = {0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0};
//static const double LIST_DEF_BR_LEN_GRID_VALS_AGTP[] = {0.0001, 0.0003, 0.0005, 0.00075, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 1.0};
//static const double LIST_DEF_BR_LEN_GRID_VALS_AGTP[] = {0.0001, 0.0002, 0.0004, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 1.0};
static vector<double> listDefSizesAGTP(LIST_DEF_BR_LEN_GRID_VALS_AGTP, LIST_DEF_BR_LEN_GRID_VALS_AGTP+sizeof(LIST_DEF_BR_LEN_GRID_VALS_AGTP)/sizeof(LIST_DEF_BR_LEN_GRID_VALS_AGTP[0]));

///////////////////////////////////////////////////////////////////////////////////////
// supporting multi-threading
const int MAX_NUM_FAC = 100;

ApproxGTPCache & ApproxGTPCache :: Instance()
{
    //
    static ApproxGTPCache inst;
    return inst;
}

ApproxGTPCache :: ApproxGTPCache()
{
    //
    Init();
}
ApproxGTPCache :: ~ApproxGTPCache()
{
#if 0
    //
    for(int i=0; i<(int)listMapdbValsThreadSafe.size(); ++i )
    {
        delete listMapdbValsThreadSafe[i];
    }
    for(int i=0; i<(int)listMapwbValsThreadSafe.size(); ++i)
    {
        delete listMapwbValsThreadSafe[i];
    }
    for(int i=0; i<(int)listCacheBranchProbThreadSafe.size(); ++i)
    {
        delete listCacheBranchProbThreadSafe[i];
    }
#endif
}

void ApproxGTPCache :: DumpStats() const
{
    //cout << "ApproxGTPCache: tot num of dVals calculations: " << numTotdVals << ", cached number: " << numCacheddVals << ", numTotwbVals: " << numTotwbVals << ", numCachedwbVals: " << numCachedwbVals << endl;
}

void ApproxGTPCache :: Init()
{
/*
    numTotdVals = 0;
    numCacheddVals = 0;
    numTotwbVals = 0;
    numCachedwbVals = 0;
    
    //
    int numThreads = GetNumThreads();
//cout << "ApproxGTPCache Init: Number of threads to use: " << numThreads << endl;
    for(int i=0; i<numThreads; ++i)
    {
        map< pair<int,int>, double> *p1 = new map< pair<int,int>, double>;
        listMapdbValsThreadSafe.push_back( p1 );
        map< pair<int,int>, double> *p2 = new map< pair<int,int>, double>;
        listMapwbValsThreadSafe.push_back( p2 );
        map< pair<pair<int,int>,int>, double > *p3 = new map< pair<pair<int,int>,int>, double >;
        listCacheBranchProbThreadSafe.push_back(p3);
    }
 */
    
    /*
    listFactorials.resize(MAX_NUM_FAC, 1.0);
    for(int i=1; i<MAX_NUM_FAC; ++i)
    {
        listFactorials[i] = listFactorials[i-1]*i;
    }*/
    
    listProdChoose2.resize(MAX_NUM_FAC);
    for(int ub=1; ub<MAX_NUM_FAC; ++ub)
    {
        listProdChoose2[ub].resize( ub );
        listProdChoose2[ub][0] = 1.0;
        if( ub >= 2 )
        {
            for(int cb = 1; cb<ub; ++cb)
            {
                listProdChoose2[ub][cb] = listProdChoose2[ub][cb-1] * 0.5 * (ub-cb+1) * (ub-cb);
            }
        }
    }
}

/*
map< pair<int,int>, double> & ApproxGTPCache :: GetdbValCache(int tid)
{
    //
    YW_ASSERT_INFO( tid<(int)listMapdbValsThreadSafe.size(), "Overflow1" );
    ++numTotdVals;
    return *listMapdbValsThreadSafe[tid];
}

map< pair<int,int>, double> & ApproxGTPCache :: GetwbValCache(int tid)
{
    //
    YW_ASSERT_INFO( tid<(int)listMapwbValsThreadSafe.size(), "Overflow2" );
    ++numTotwbVals;
    return *listMapwbValsThreadSafe[tid];
}

map< pair<pair<int,int>,int>, double > &  ApproxGTPCache :: GetBrProbCache(int tid)
{
    //
    YW_ASSERT_INFO( tid<(int)listCacheBranchProbThreadSafe.size(), "Overflow3" );
    return *listCacheBranchProbThreadSafe[tid];
}
*/

/*
double ApproxGTPCache :: GetFac(int k) const
{
    if( k >=  MAX_NUM_FAC)
    {
        return MAX_DOUBLE_VAL;
    }
    return listFactorials[k];
} */

double ApproxGTPCache :: GetProdChoose2(int ub, int cb) const
{
    //
    if( ub >=  MAX_NUM_FAC)
    {
        return MAX_DOUBLE_VAL;
    }
    return listProdChoose2[ub][cb];
}

///////////////////////////////////////////////////////////////////////////////////////
// pre-compute some commonly used coalescent terms
const int MAX_PRE_HAP_NUM = 20;
const int MAX_LRU_CACHE_SIZE = 10000;

ApproxGTPCache2 & ApproxGTPCache2 :: Instance()
{
    static ApproxGTPCache2 inst;
    return inst;
}
void ApproxGTPCache2 :: InitNumThreads(int nt)
{
    listCachedCoalProb.clear();
    listCachedCoalProb.resize(nt);
    
    // clear out all LRU caches too
    for(unsigned int i=0; i<listCachedCoalProbLRU.size(); ++i)
    {
        delete listCachedCoalProbLRU[i];
    }
    listCachedCoalProbLRU.clear();
    for(int i=0; i<nt; ++i)
    {
        auto x = new lru_cache2< pair<pair<int,int>,int>, double, KeyHasher >;
        x->setmaxsize(MAX_LRU_CACHE_SIZE);
        listCachedCoalProbLRU.push_back( x );
    }
}
void ApproxGTPCache2 :: Init(int maxNumLins)
{
    listdbVals.resize(maxNumLins+1);
    //
    for(int ub=1; ub<=maxNumLins; ++ub)
    {
        listdbVals[ub].resize( ub+1 );
        for(int cb=0; cb<ub; ++cb)
        {
            listdbVals[ub][cb] = ApproxGTPCache :: Instance().GetProdChoose2(ub, cb);
        }
    }
    
    
#if 0
    listwbVals.resize(maxNumLins+1);
    //
    for(int ub=1; ub<=maxNumLins; ++ub)
    {
        listwbVals[ub].resize(ub+1);
        for(int cb=0; cb<=ub; ++cb)
        {
            double res = 1.0;
#if 0   // another approximation, making each coalescent history only has one way to be consistent with the gene tree topology (i.e., behave like a caterpillar tree
            int numEvts = cb;
            
            //res = ApproxGTPCache::Instance().GetFac(numEvts);
            
            if( ub >= 0 && cb >= 0)
            {
                res *= 1.0/ApproxGTPCache :: Instance().GetProdChoose2(ub, cb);
            }
#endif
            
            // now include the fake branch factors
            // simplifying assumption: all cb events involves ub-cb+1 lineages
            //int vb = ub-cb;
            //int numLinsNotUsed = vb-1;
            //res *= CalcBranchFactorAtFake( 1, ub-numLinsNotUsed, ub-numLinsNotUsed );
            //res *= 1.0/ApproxGTPCache::Instance().GetFac(cb);
            listwbVals[ub][cb] = res;
#if 0
            double res = 1.0;
            int numEvts = cb;
            res = ApproxGTPCache::Instance().GetFac(numEvts);
            
            if( ub >= 0 && cb >= 0)
            {
                res *= 1.0/ApproxGTPCache :: Instance().GetProdChoose2(ub, cb);
            }
            
            // now include the fake branch factors
            // simplifying assumption: all cb events involves ub-cb+1 lineages
            int vb = ub-cb;
            int numLinsNotUsed = vb-1;
            //res *= CalcBranchFactorAtFake( 1, ub-numLinsNotUsed, ub-numLinsNotUsed );
            res *= 1.0/ApproxGTPCache::Instance().GetFac(cb);
            listwbVals[ub][cb] = res;
#endif
        }
    }
#endif

#if 0
#define AGTP_USE_CACHED_BRANCH_PROB 1
#ifdef AGTP_USE_CACHED_BRANCH_PROB
    // now calculate branch prob
    
    //if( maxNumLins <= MAX_PRE_HAP_NUM )
    //{
    //listCacheBranchProb.resize( GetMaxBrProbPreValu()+1 );
    listCacheBranchProb.resize( MAX_PRE_HAP_NUM+1 );
    for(unsigned int u=0; u<listCacheBranchProb.size(); ++u)
    {
        listCacheBranchProb[u].resize(u+1);
        for(unsigned int v=1; v<=u; ++v )
        {
            listCacheBranchProb[u][v].resize( numBranchLenDigitPreCalc );
            for(unsigned int k=0; k<listCacheBranchProb[u][v].size(); ++k)
            {
                double len = ConvGridToDist(k);
                listCacheBranchProb[u][v][k] = UtilsCalcCoalProbBranch(u, v, len);
            }

        }
    }
    //}
#endif
#endif
}
double ApproxGTPCache2 :: GetdbValCache(int ub, int cb)
{
    return listdbVals[ub][cb];
}
double ApproxGTPCache2 :: GetwbValCache(int ub, int cb)
{
    return 1.0;
    //return listwbVals[ub][cb];
}

/*
double ApproxGTPCache2 :: GetBrProbCache(int u, int v, double brLen)
{
    int grid = ConvDistToGrid(brLen);
    return listCacheBranchProb[u][v][grid];
}*/

int ApproxGTPCache2 :: GetMaxBrProbPreValu() const
{
    return 8;
}

ApproxGTPCache2 :: ApproxGTPCache2()
{
    // by default, only one thread
    InitNumThreads(1);
}
ApproxGTPCache2 :: ~ApproxGTPCache2()
{
    for(unsigned int i=0; i<listCachedCoalProbLRU.size(); ++i)
    {
        delete listCachedCoalProbLRU[i];
    }
}

int ApproxGTPCache2 ::  ConvDistToGrid(double dist) const
{
    int res = (int)(dist*(numBranchLenDigitPreCalc-1)+0.5);
    //int res = (int)(dist*(pow(10.0, numBranchLenDigit))+0.5);
    if( res >=numBranchLenDigitPreCalc)
    {
        res = numBranchLenDigitPreCalc-1;
    }
    return res;
}
double ApproxGTPCache2 :: ConvGridToDist(int grid) const
{
    return grid * (1.0/(numBranchLenDigitPreCalc-1));
    //return grid * (1.0/pow(10.0, numBranchLenDigit));
}
/*
double ApproxGTPCache2 :: GetCoalProbBtwCfgLins( int u, int v, double brLen )
{
    // use cached prob if it is has been pre-computed
    if( u <= MAX_PRE_HAP_NUM )
    {
        return GetBrProbCache(u, v, brLen);
    }
    
    if( u < v || v == 0)
    {
        // impossible
        return MAX_NEG_DOUBLE_VAL;
    }
    pair<pair<int,int>,int> pp(std::make_pair(u,v), ConvDistToGrid(brLen));
    auto it = listCachedCoalProb.find(pp);
    if( it != listCachedCoalProb.end() )
    {
        return it->second;
    }
    // claculate it first
    int cb = u-v;
    double res = 1.0;
    res = 1.0/GetdbValCache(u, cb);
//cout << "** coeff after dbVal: " << res << endl;
//#endif
    res *= GetwbValCache( u, cb);
//cout << "coeff after wbVal: " << res << endl;
    
    // get branch if still worth it (value not too small)
    if( res >= 1.0e-300)
    {
//cout << "u:" << u << ",v:" << v << ", brLen:" << brLen << endl;
        res *= UtilsCalcCoalProbBranch2(u, v, brLen);
    }
    listCachedCoalProb[pp] = res;
//cout << "coeff after br prob: " << res << endl;
    return res;
}*/
double ApproxGTPCache2 :: GetCoalProbBtwCfgLins2( int u, int v, double brLen, int tid )
{
    // don't use cache
    if( u < v || v == 0)
    {
        // impossible
        return MAX_NEG_DOUBLE_VAL;
    }
    
    int grLen = ConvDistToGrid(brLen);
    double brLenUse = ConvGridToDist(grLen);
    
    // no-lock because each thread access different cache
    //static std::mutex mut;
    pair<pair<int,int>,int> pp(std::make_pair(u,v), grLen);
    double res = 1.0;
    //mut.lock();
    auto it = listCachedCoalProb[tid].find(pp);
    if( it != listCachedCoalProb[tid].end() )
    {
        res = it->second;
        //mut.unlock();
        return res;
    }
    //mut.unlock();
    
    // claculate it first
    int cb = u-v;
    //double res = 1.0;
    double dval = ApproxGTPCache :: Instance().GetProdChoose2(u, cb);
    res = 1.0/dval;
//cout << "** coeff after dbVal: " << res << endl;
//#endif
    //res *= GetwbValCache( u, cb);
//cout << "coeff after wbVal: " << res << endl;
    
    // get branch if still worth it (value not too small)
    if( res >= 1.0e-300)
    {
//cout << "u:" << u << ",v:" << v << ", brLen:" << brLen << endl;
        //res *= UtilsCalcCoalProbBranch2(u, v, brLen);
        res *= UtilsCalcCoalProbBranch2(u, v, brLenUse);
    }
//cout << "coeff after br prob: " << res << endl;
    //mut.lock();
    listCachedCoalProb[tid][pp] = res;
    //mut.unlock();
    return res;
}
double ApproxGTPCache2 :: GetCoalProbBtwCfgLins2ConvertedLen( int u, int v, double brLen, int grLen, int tid )
{
    // don't use cache
    if( u < v || v == 0)
    {
        // impossible
        return MAX_NEG_DOUBLE_VAL;
    }
    
    double res = 1.0;
    
    
    // no-lock since each thread access its own cache
    //static std::mutex mut;
    pair<pair<int,int>,int> pp(std::make_pair(u,v), grLen);
    //mut.lock();
    auto it = listCachedCoalProb[tid].find(pp);
    if( it != listCachedCoalProb[tid].end() )
    {
        res = it->second;
        //mut.unlock();
        return res;
    }
    //mut.unlock();
    
    
    // claculate it first
    int cb = u-v;
    //double res = 1.0;
    double dval = ApproxGTPCache :: Instance().GetProdChoose2(u, cb);
    res = 1.0/dval;
//cout << "** coeff after dbVal: " << res << endl;
//#endif
    //res *= GetwbValCache( u, cb);
//cout << "coeff after wbVal: " << res << endl;
    
    // get branch if still worth it (value not too small)
    if( res >= 1.0e-300)
    {
//cout << "u:" << u << ",v:" << v << ", brLen:" << brLen << endl;
        //res *= UtilsCalcCoalProbBranch2(u, v, brLen);
        res *= UtilsCalcCoalProbBranch2(u, v, brLen);
    }
//cout << "coeff after br prob: " << res << endl;
    //mut.lock();
//#if 0
    listCachedCoalProb[tid][pp] = res;
//#endif
    //mut.unlock();
    return res;
}

double ApproxGTPCache2 :: GetCoalProbBtwCfgLins2ConvertedLenLRU( int u, int v, double brLen, int grLen, int tid )
{
    // don't use cache
    if( u < v || v == 0)
    {
        // impossible
        return MAX_NEG_DOUBLE_VAL;
    }
    
    double res = 1.0;
    
    // no-lock since each thread access its own cache
    //static std::mutex mut;
    pair<pair<int,int>,int> pp(std::make_pair(u,v), grLen);
    //mut.lock();
    bool f = listCachedCoalProbLRU[tid]->exists(pp, res);
    if( f )
    {
        //res = listCachedCoalProbLRU[tid]->get(pp);
        //mut.unlock();
        return res;
    }
    //mut.unlock();
    
    // claculate it first
    int cb = u-v;
    //double res = 1.0;
    double dval = ApproxGTPCache :: Instance().GetProdChoose2(u, cb);
    res = 1.0/dval;
//cout << "** coeff after dbVal: " << res << endl;
//#endif
    //res *= GetwbValCache( u, cb);
//cout << "coeff after wbVal: " << res << endl;
    
    // get branch if still worth it (value not too small)
    if( res >= 1.0e-300)
    {
//cout << "u:" << u << ",v:" << v << ", brLen:" << brLen << endl;
        //res *= UtilsCalcCoalProbBranch2(u, v, brLen);
        res *= UtilsCalcCoalProbBranch2(u, v, brLen);
    }
//cout << "coeff after br prob: " << res << endl;
    //mut.lock();
//#if 0
    listCachedCoalProbLRU[tid]->put(pp, res);
//#endif
    //mut.unlock();
    return res;
}

void ApproxGTPCache2 :: DumpStats() const
{
    int nt=0;
    for(unsigned int i=0; i<listdbVals.size(); ++i)
    {
        nt += listdbVals[i].size();
    }
    cout << "ApproxGTPCache2: Num of listdbVals entries: " << nt << endl;
    
    /*
    nt = 0;
    for(unsigned int i=0; i<listCacheBranchProb.size(); ++i)
    {
        for(unsigned int j=0; j<listCacheBranchProb[i].size(); ++j)
        {
            nt += listCacheBranchProb[i][j].size();
        }
    }
    cout << "ApproxGTPCache2: Num of listCacheBranchProb entries: " << nt << endl;  */
    
    nt = 0;
    for(unsigned int i=0; i<listCachedCoalProb.size(); ++i)
    {
        nt += listCachedCoalProb[i].size();
    }
    cout << "ApproxGTPCache2: Num of listCachedCoalProb entries: " << nt << endl;
}

///////////////////////////////////////////////////////////////////////////////////////
// for benchmark run time


ApproxGTPStats & ApproxGTPStats :: Instance()
{
    static ApproxGTPStats inst;
    return inst;
}

void ApproxGTPStats :: DumpStats() const
{
    //
    cout << "************************ Execution Statistics **************************\n";
    cout << "      Maximum size of configuration list: " << maxSzACList << endl;
    cout << "      Total number of gene probability computation (for all given input trees): " << listProbCalcTimeEnd.size() << endl;
    cout << "      Total number of branch optimization: " << brOptNums << endl;
    double tmTot = 0.0;
    for(int i=0; i<(int)listProbCalcTimeEnd.size(); ++i )
    {
        tmTot += (listProbCalcTimeEnd[i] - listProbCalcTimeStart[i])/(double)CLOCKS_PER_SEC;
    }
    if( listProbCalcTimeEnd.size() > 0 )
    {
        cout << "     Average time for computing the gene tree probability is " << tmTot/listProbCalcTimeEnd.size() << " seconds\n";
    }
}

void ApproxGTPStats :: RecordProbComputeStart()
{
    //
    //auto duration = now.time_since_epoch();
    //auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    listProbCalcTimeStart.push_back( GetCurrentTimeTick() );
}

void ApproxGTPStats :: RecordProbComputeEnd()
{
    //
    //auto duration = now.time_since_epoch();
    //auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    listProbCalcTimeEnd.push_back( GetCurrentTimeTick() );
}

void ApproxGTPStats :: RecordMaxACSize(int szACList)
{
    if( szACList > maxSzACList )
    {
        maxSzACList = szACList;
    }
}

void ApproxGTPStats :: RecordBranchLenOpt()
{
    ++brOptNums;
}

ApproxGTPStats :: ApproxGTPStats() : maxSzACList(0), brOptNums(0)
{
    //
}


//////////////////////////////////////////////////////////////////////////////////
// For compact helper

CompactApproxGeneTreeProbHelper :: CompactApproxGeneTreeProbHelper( MarginalTree &treeSpeciesIn, PhylogenyTreeBasic &treeGeneIn ) : treeSpecies(treeSpeciesIn), treeGene(treeGeneIn)
{
    ConsClusterInfo();
}

void CompactApproxGeneTreeProbHelper :: ConsClusterInfo()
{
    // now traverse the tree and find the MRCA
    // simple approach: traverse the tree; for each node,
    // find its left and right and then set their MRCA to the current one
    map<int, TreeNode*> mapIdNode;
    PhylogenyTreeIterator itorTree(treeGene);
    itorTree.Init();
    //treeGene.InitPostorderWalk();
    while(itorTree.IsDone() == false )
    {
        TreeNode *pn = itorTree.GetCurrNode( ) ;
        itorTree.Next();
        if( pn == NULL )
        {
            break;      // done with all nodes
        }
        int idSelf = pn->GetID();
        
        // keep track of id-->node info
        //cout << "In Gene tree traversal, find a node idSelf: " << idSelf << endl;
        mapIdNode.insert( map<int, TreeNode*> :: value_type( idSelf, pn )  );
    }
    
    // YW: 05/18/16: do we really need this?
    this->treeSpecies.BuildDescendantInfo();
    
#if 0
cout << "ConsClusterInfo: ";
cout << "Current species tree: " << treeSpecies.GetNewick() << endl;
treeSpecies.Dump();
cout << "Current gene tree: ";
string strGT;
treeGene.ConsNewick(strGT);
cout << strGT << endl;
treeGene.Dump();
#endif
    
// print out some info
//for(int sbIndex = 0; sbIndex < treeSpecies.GetTotNodesNum(); ++sbIndex )
//{
//cout << "ST node: " << sbIndex << ", left chilld: " << treeSpecies.GetLeftDescendant(sbIndex) << ", right child: " << treeSpecies.GetRightDescendant(sbIndex) << endl;
//}

    // collect gene alleles for taxon
    listTaxonGids.resize( treeSpecies.GetNumLeaves() );
    for(int t = 0; t<(int)listTaxonGids.size(); ++t)
    {
        //cout << "GetGeneAllelesForSpecies: taxon = " << taxon << endl;
        // get all alleles in GT for a gene
        char idbuf[100];
        snprintf(idbuf, 100, "%d", t);
        string idstr = idbuf;
        // search for each leaf, to find match
        treeGene.GetLeavesIdsWithLabel(  idstr, listTaxonGids[t] );
    }

    //
    listSTLinClusterList.clear();
    //listSTLinClusterRootMap.clear();
    
    // construct cluster info for species tree branch
    // approach: consider each species tree branch; find their taxa; get the corresponding
    // gene lineages (taxa); find the clusters of them
    vector<set<int> > leafNodeLabels;
    treeSpecies.ConsDecedentLeavesInfoLabels( leafNodeLabels );
    
    for( int sbIndex = 0; sbIndex < treeSpecies.GetTotNodesNum(); ++sbIndex )
    {
        //
        set<int> &setTaxaBelow = leafNodeLabels[sbIndex];
        
        set<int> geneAlleles;
        GetGeneAllelesForSpeciesSet( setTaxaBelow, geneAlleles);
        
        // find all nodes
        set<TreeNode *> setLeaves;
        GetNodesForIds(geneAlleles, mapIdNode, setLeaves);
        
        // find all the subtrees clustered
        set< set<TreeNode *> > setSubtreeClades;
        treeGene.FindCladeOfSubsetLeavesExact( setLeaves, setSubtreeClades );
        set< set<TreeNode *> > setSubtreeCladesMaximal;
        PhylogenyTreeBasic :: GroupLeavesToSubtrees( setLeaves,  setSubtreeClades, setSubtreeCladesMaximal);
        
        // collect all roots
        vector<ApproxGeneTreeProbCluster> vecClusters;
        for( set<set<TreeNode *> > :: iterator it = setSubtreeCladesMaximal.begin(); it != setSubtreeCladesMaximal.end(); ++it )
        {
            TreeNode *pp = treeGene.GetSubtreeRootForLeaves( *it );
            // create subtrees
            ApproxGeneTreeProbCluster clus( pp );
            vecClusters.push_back(clus);
        }
        listSTLinClusterList.push_back( vecClusters );
//cout << "For species tree branch: " << sbIndex << ", number of clusters: " << vecClusters.size() << endl;
    }
    
#if 0
    // also update the other
    listSTLinClusterRootMap.clear();
    listSTLinClusterRootMap.resize( listSTLinClusterList.size() );
    for(int i=0; i<(int)listSTLinClusterRootMap.size(); ++i)
    {
        for(int jj=0; jj<(int)listSTLinClusterList[i].size(); ++jj )
        {
            TreeNode *prootc = listSTLinClusterList[i][jj].GetClusterRoot();
            listSTLinClusterRootMap[i].insert( map<TreeNode *,int> :: value_type( prootc, jj ) );
        }
    }
#endif

    
    // find the set of lineage numbers under each ST node
    listSTLinClusterSizeList.resize( listSTLinClusterList.size() );
    for(int i=0; i<(int)this->listSTLinClusterList.size(); ++i)
    {
        for(int j=0; j<(int)listSTLinClusterList[i].size(); ++j)
        {
            TreeNode *pn = listSTLinClusterList[i][j].GetClusterRoot();
            set<TreeNode *> setLeaves;
            pn->GetAllLeavesUnder( setLeaves );
            listSTLinClusterSizeList[i].push_back( setLeaves.size() );
        }
    }
    

#if 0
cout << "List of lower bound of lineages: \n";
for(int i=0; i<(int)listSTLinClusterSizeListMin.size(); ++i)
{
cout << "ST ndoe: " << i << ": list of cluster min-sizes: ";
DumpIntVec( listSTLinClusterSizeListMin[i] );
}
#endif
}

void CompactApproxGeneTreeProbHelper :: GetGeneAllelesForSpeciesSet( const set<int> &setTaxa, set<int> &geneAlleles)
{
    //
    geneAlleles.clear();
    for( set<int> :: const_iterator it = setTaxa.begin(); it != setTaxa.end(); ++it)
    {
        //set<int> ss;
        //GetGeneAllelesForSpecies( *it, ss );
        //UnionSets( geneAlleles, ss );
        UnionSets( geneAlleles, GetGeneAllelesForSpecies( *it) );
    }
}

const set<int> & CompactApproxGeneTreeProbHelper :: GetGeneAllelesForSpecies( int taxon)
{
    return listTaxonGids[taxon];
}

void CompactApproxGeneTreeProbHelper :: GetNodesForIds(const set<int> &nids, const map<int, TreeNode*> &mapIdNode, set<TreeNode *> &setNodes) const
{
    //
    //ApproxGeneTreeProbHelper *pthis = const_cast<ApproxGeneTreeProbHelper *>(this);
    setNodes.clear();
    for(set<int> :: iterator it = nids.begin(); it != nids.end(); ++it)
    {
        map<int, TreeNode*> :: const_iterator it2 = mapIdNode.find(*it);
        //YW_ASSERT_INFO(it2 != mapIdNode.end(), "Fail to find");
        TreeNode *pn = it2->second;
        setNodes.insert(pn);
    }
}


//////////////////////////////////////////////////////////////////////////////////
// Heuristic prob calculation

ApproxGeneTreeProbHeu2 :: ApproxGeneTreeProbHeu2(MarginalTree &treeSpeciesIn, const vector<PhylogenyTreeBasic *> &listGeneTreePtrsIn, int tidIn) : treeSpecies(treeSpeciesIn), listGeneTreePtrs(listGeneTreePtrsIn), threadIndex(tidIn)
{
    //strST = treeSpecies.GetNewickSorted(true);
    
    Init();
//cout << "ApproxGeneTreeProbHeu2::Done init\n";
    // snap branch lengths of ST
    //SnapSpeciesTreeLens();
    
    // Init prob calculation
    InitProbs();
//cout << "ApproxGeneTreeProbHeu2::Done init probs\n";
}

ApproxGeneTreeProbHeu2 :: ~ApproxGeneTreeProbHeu2()
{
    for(auto x: listGTPtHelpers)
    {
        delete x;
    }
    listGTPtHelpers.clear();
}

void ApproxGeneTreeProbHeu2 :: CalcLogProb(vector<double> &listLogProbs)
{
//cout << "STApproxGeneTreeProb2: CalcLogProb:\n";
    // standarize br lengths
    //SnapSpeciesTreeLens();
    
    // compute all the gene trees in parallele since this will be better for cache performance (?)
//cout << "CalcLogProb:\n";
    listLogProbs.resize(GetNumGTs());
    
//cout << "Num of root cfgs (lower): " << setCfgsRootLower.size() << endl;
    // calculate prob based on this
    // this should be the root
    // here each config is going to coalesce into the final (infinite) branch out of species tree root
    //YW_ASSERT_INFO( setCfgsRootLower.size() > 0, "Root cfgs: not constructed"  );
    for(unsigned int tr=0; tr<GetNumGTs(); ++tr)
    {
        double logprob = CalcLogProbForTree(tr);
        listLogProbs[tr] = logprob;
    }
    
//cout << "ApproxGeneTreeProbHeu2 :: CalcLogProb: species tree: " << treeSpecies.GetNewickSorted(true) << ", listLogProbs: ";
//DumpDoubleVec(listLogProbs);
}
/*
static void UtilAGUpdateSTBrLenMT( ApproxGeneTreeProbHeu2 *pProbCalc2, int numTotTrees, vector<bool> *pvecProcessed, int nodeST, MarginalTree *ptreeSpecies )
{
    int indTree = 0;
    static std::mutex mut;
    int numSTNodes = ptreeSpecies->GetTotNodesNum();
    while(indTree < numTotTrees)
    {
        // find a node that is not being processed
        mut.lock();
        while( indTree < numTotTrees && (*pvecProcessed)[indTree] == true )
        {
            ++indTree;
        }
        if( indTree < numTotTrees )
        {
            (*pvecProcessed)[indTree] = true;
        }
        mut.unlock();
        if( indTree >= numTotTrees )
        {
            break;
        }
        
        // real work
        int nodeCur = nodeST;
        while( ptreeSpecies->GetRoot() != nodeCur )
        {
            int nodePar = ptreeSpecies->GetParent(nodeCur);
            
            // construct upper DP first for the current node
            pProbCalc2->CalcUpperProbAt(indTree, nodeCur);
            
            // now construct the lower DP for the parent of the current node
            pProbCalc2->ConsLowerProbAt(indTree, nodePar);
            
            nodeCur = nodePar;
        }
 
    }
}


void ApproxGeneTreeProbHeu2 :: UpdateSTBrLenMTImp(int nodeST)
{
    int numThreadsUse = numThreads;
    if( (int)listGeneTreePtrs.size() < numThreadsUse )
    {
        numThreadsUse = listGeneTreePtrs.size();
    }
    vector<bool> vecProcessed(listGeneTreePtrs.size());
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        vecProcessed[i] = false;
    }
    // vector for a node being processed
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        thread *pthr = new thread(UtilAGUpdateSTBrLenMT, this, listGeneTreePtrs.size(), &vecProcessed, nodeST, &(this->treeSpecies) );
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
*/
void ApproxGeneTreeProbHeu2 :: UpdateSTBrLen(int nodeST, double brLenNew)
{
    //YW_ASSERT_INFO( nodeST >=0 && nodeST <treeSpecies.GetTotNodesNum()-1, "nodeST: overflow" );
    
    // change branch length
    treeSpecies.SetBranchLen(nodeST, brLenNew);
    
    // calculate prob
    //CalcLogProb(listLogProbs);
    //InitProbs();
//#if 0
    
    /*
    if( numThreads > 1 )
    {
//cout << "Multithreading: update brLen\n";
        UpdateSTBrLenMTImp(nodeST);
    }
    else
    {  */
        // re-init probs for those parts that are affected
        // now deal with each internal node
        int nodeCur = nodeST;
        while( treeSpecies.GetRoot() != nodeCur )
        {
            int nodePar = treeSpecies.GetParent(nodeCur);
            
            // construct upper DP first for the current node
            for(int tr=0; tr<GetNumGTs(); ++tr)
            {
                CalcUpperProbAt(tr, nodeCur);
                
                // now construct the lower DP for the parent of the current node
                ConsLowerProbAt(tr, nodePar);
            }
            
            nodeCur = nodePar;
        }
    //}
//#endif
}

void ApproxGeneTreeProbHeu2 :: OnUpdateSTBrLen()
{
//cout << "SGSTApproxGeneTreeProb : OnUpdateSTBrLen" << endl;
    //SnapSpeciesTreeLens();
    //InitProbs();
}

string ApproxGeneTreeProbHeu2 :: GetSpeciesTreeNW() const
{
    //return strST;
    return treeSpecies.GetNewickSorted(true);
}

/*
void ApproxGeneTreeProbHeu2 :: SnapSpeciesTreeLens()
{
    for(int node=0; node<treeSpecies.GetTotNodesNum()-1; ++node)
    {
        SnapSpeciesTreeLenAt(node);
    }
}

void ApproxGeneTreeProbHeu2 :: SnapSpeciesTreeLenAt(int node)
{
    double len = treeSpecies.GetEdgeLen(node);
    int pos = GetClosestTo2( listDefSizesAGTP, len );
    double lenConv = listDefSizesAGTP[pos];
    treeSpecies.SetBranchLen(node, lenConv);
}  */

// utility code
void ApproxGeneTreeProbHeu2 :: Init()
{
    // init helper
    for(unsigned int i=0; i<listGeneTreePtrs.size(); ++i)
    {
        CompactApproxGeneTreeProbHelper *pHelper = new CompactApproxGeneTreeProbHelper(treeSpecies, *listGeneTreePtrs[i]);
        listGTPtHelpers.push_back(pHelper);
    }
    
    // init max lins
    listMaxLinGTLinNumAt.resize(listGeneTreePtrs.size());
    for(unsigned int i=0; i<listGeneTreePtrs.size(); ++i)
    {
        for(int nv = 0; nv<treeSpecies.GetTotNodesNum(); ++nv)
        {
            listMaxLinGTLinNumAt[i].push_back( GetMaxLinGTLinNumAt(i, nv) );
        }
    }
}
/*
static void UtilAGTPHeu2InitProbMT( ApproxGeneTreeProbHeu2 *pProbCalc2, int numTotTrees, vector<bool> *pvecProcessed,  MarginalTree *ptreeSpecies, vector<vector<vector<double> > > *ptblDPLower, vector<vector<vector<double> > > *ptblDPUpper  )
{
    int indTree = 0;
    static std::mutex mut;
    int numSTNodes = ptreeSpecies->GetTotNodesNum();
    while(indTree < numTotTrees)
    {
        // find a node that is not being processed
        mut.lock();
        while( indTree < numTotTrees && (*pvecProcessed)[indTree] == true )
        {
            ++indTree;
        }
        if( indTree < numTotTrees )
        {
            (*pvecProcessed)[indTree] = true;
        }
        mut.unlock();
        if( indTree >= numTotTrees )
        {
            break;
        }
        
        (*ptblDPLower)[indTree].resize( numSTNodes );
        (*ptblDPUpper)[indTree].resize( numSTNodes );
        for(int nodeST=0; nodeST<numSTNodes; ++nodeST)
        {
            //int minLinNum = GetMinLinGTLinNumAt(tr, nodeST);
            //int maxLinNum = GetMaxLinGTLinNumAt(tr, nodeST);
            (*ptblDPLower)[indTree][nodeST].resize( pProbCalc2->GetNumDPCellsAt(indTree, nodeST) );
            // init to very small num
            for(unsigned int i=0; i<(*ptblDPLower)[indTree][nodeST].size(); ++i)
            {
                (*ptblDPLower)[indTree][nodeST][i] = MAX_NEG_DOUBLE_VAL;
            }
            
            // same thing for DP upper
            (*ptblDPUpper)[indTree][nodeST].resize( (*ptblDPLower)[indTree][nodeST].size() );
            // init to very small num
            for(unsigned int i=0; i<(*ptblDPUpper)[indTree][nodeST].size(); ++i)
            {
                (*ptblDPUpper)[indTree][nodeST][i] = MAX_NEG_DOUBLE_VAL;
            }
        }
        // init for each leaf
        for(int nodeST=0; nodeST<pProbCalc2->GetNumTaxa(); ++nodeST)
        {
            // set the max to be prob-1
            (*ptblDPLower)[indTree][nodeST][ (int)(*ptblDPLower)[indTree][nodeST].size()-1 ] = 0.0;
        }
        // now deal with each internal node
        for(int nodeST=pProbCalc2->GetNumTaxa(); nodeST<numSTNodes; ++nodeST)
        {
            int childLeft = ptreeSpecies->GetLeftDescendant(nodeST);
            int childRight = ptreeSpecies->GetRightDescendant(nodeST);
            
            // construct upper DP first for two children
            pProbCalc2->CalcUpperProbAt(indTree, childLeft);
            // then right
            pProbCalc2->CalcUpperProbAt(indTree, childRight);
            
            // now construct the lower DP for the current node
            //int numLinsMin = GetMinLinGTLinNumAt(indexGT, nodeST);
            // consider all possible split
            pProbCalc2->ConsLowerProbAt(indTree, nodeST);
        }
    }
}

void ApproxGeneTreeProbHeu2 :: InitDPTblMT()
{
    //
    int numThreadsUse = numThreads;
    if( (int)listGeneTreePtrs.size() < numThreadsUse )
    {
        numThreadsUse = listGeneTreePtrs.size();
    }
    vector<bool> vecProcessed(listGeneTreePtrs.size());
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        vecProcessed[i] = false;
    }
    // vector for a node being processed
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        thread *pthr = new thread(UtilAGTPHeu2InitProbMT, this, listGeneTreePtrs.size(), &vecProcessed, &(this->treeSpecies), &(this->tblDPLower), &(this->tblDPUpper) );
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
*/

void ApproxGeneTreeProbHeu2 :: InitProbs()
{
    // init prob calculation
    tblDPLower.resize( GetNumGTs() );
    tblDPUpper.resize( GetNumGTs() );
    
    /*
    if( numThreads > 1 )
    {
//cout << "Multithreading: init\n";
        InitDPTblMT();
    }
    else
    {  */
        // prepare DP tables for all gene trees
        int numSTNodes = treeSpecies.GetTotNodesNum();
        for(int tr=0; tr<GetNumGTs(); ++tr)
        {
            tblDPLower[tr].resize( numSTNodes );
            tblDPUpper[tr].resize( numSTNodes );
            for(int nodeST=0; nodeST<numSTNodes; ++nodeST)
            {
                //int minLinNum = GetMinLinGTLinNumAt(tr, nodeST);
                //int maxLinNum = GetMaxLinGTLinNumAt(tr, nodeST);
                tblDPLower[tr][nodeST].resize( GetNumDPCellsAt(tr, nodeST) );
                // init to very small num
                for(unsigned int i=0; i<tblDPLower[tr][nodeST].size(); ++i)
                {
                    tblDPLower[tr][nodeST][i] = MAX_NEG_DOUBLE_VAL;
                }
                
                // same thing for DP upper
                tblDPUpper[tr][nodeST].resize(tblDPLower[tr][nodeST].size() );
                // init to very small num
                for(unsigned int i=0; i<tblDPUpper[tr][nodeST].size(); ++i)
                {
                    tblDPUpper[tr][nodeST][i] = MAX_NEG_DOUBLE_VAL;
                }
            }
            // init for each leaf
            for(int nodeST=0; nodeST<GetNumTaxa(); ++nodeST)
            {
                // set the max to be prob-1
                tblDPLower[tr][nodeST][ (int)tblDPLower[tr][nodeST].size()-1 ] = 0.0;
            }
            // now deal with each internal node
            for(int nodeST=GetNumTaxa(); nodeST<numSTNodes; ++nodeST)
            {
                // construct upper DP first for two children: left first
                int childLeft = treeSpecies.GetLeftDescendant(nodeST);
                CalcUpperProbAt(tr, childLeft);
                
                // then right
                int childRight = treeSpecies.GetRightDescendant(nodeST);
                CalcUpperProbAt(tr, childRight);
                
                // now construct the lower DP for the current node
                //int numLinsMin = GetMinLinGTLinNumAt(indexGT, nodeST);
                // consider all possible split
                ConsLowerProbAt(tr, nodeST);
            }
        }
    //}
    
//cout << "now init listLogCoeffs...\n";
    
    // init coal coefficients
    listLogCoeffs.resize( GetNumGTs() );
    int mroot = treeSpecies.GetRoot();
    for(int tr = 0; tr<GetNumGTs(); ++tr)
    {
        listLogCoeffs[tr].resize( tblDPLower[tr][mroot].size() );
        for( int nl=0; nl<tblDPLower[tr][mroot].size(); ++nl )
        {
            int idConv = ConvDPCellIndexToLinNum(tr, mroot,nl);
//cout << "idConv:" << idConv << ", nl: " << nl << endl;
            double coeff = CalcCoalCoeff(idConv, 1);
            if( coeff <= 1.0e-300)
            {
                coeff = 1.0e-300;
            }
            listLogCoeffs[tr][nl] = log( coeff );
        }
    }
//cout << "InitProbs: done\n";
}

double ApproxGeneTreeProbHeu2 :: CalcLogProbForTree(int indexGT)
{
    // DP for calculating prob of the gene tree
    int mroot = treeSpecies.GetRoot();

    // now summerize at the root
    double res = MAX_NEG_DOUBLE_VAL;
    vector<double> listLogPs;
    for( int nl=0; nl<tblDPLower[indexGT][mroot].size(); ++nl )
    {
#if 0
        double coeff = CalcCoalCoeff(ConvDPCellIndexToLinNum(indexGT, mroot,nl), 1);
        //res = GetLogSumOfTwo(res,  log(coeff) + tblDPLower[indexGT][mroot][nl] );
        listLogPs.push_back(log(coeff) + tblDPLower[indexGT][mroot][nl] );
#endif
//cout << "coeff: " << listLogCoeffs[indexGT][nl] << ", lower: " << tblDPLower[indexGT][mroot][nl] << endl;
        listLogPs.push_back( listLogCoeffs[indexGT][nl] + tblDPLower[indexGT][mroot][nl] );
    }
//cout << "ApproxGeneTreeProbHeu2 :: CalcLogProbForTree: indexGT" << indexGT << ", listLogPs: ";
//DumpDoubleVec(listLogPs);
    res = GetLogSumOfLogs2( listLogPs );
    return res;
}

int ApproxGeneTreeProbHeu2 :: GetMaxLinGTLinNumAt(int indexGT, int nodeST) const
{
    const vector<int> &listCnts = listGTPtHelpers[indexGT]->GetCfgBoundsAt(nodeST);
    int res = 0;
    for(unsigned int i=0; i<listCnts.size(); ++i)
    {
        res += listCnts[i];
    }
    return res;
}

int ApproxGeneTreeProbHeu2 :: GetMinLinGTLinNumAt(int indexGT, int nodeST) const
{
    return listGTPtHelpers[indexGT]->GetNumAncClusters(nodeST);
}

// prob of u lineages coalescing into v along a branch of length brLen
double ApproxGeneTreeProbHeu2 :: CalcCoalProbFor(double brLen, int u, int v)
{
    //return ApproxGTPCache2 :: Instance().GetCoalProbBtwCfgLins(u, v, brLen);
    return ApproxGTPCache2 :: Instance().GetCoalProbBtwCfgLins2(u, v, brLen, threadIndex);
}
double ApproxGeneTreeProbHeu2 :: CalcCoalProbForLenGrid(int grLen, double brLen, int u, int v)
{
    // grLen: grid of length
    //return ApproxGTPCache2 :: Instance().GetCoalProbBtwCfgLins(u, v, brLen);
    //return ApproxGTPCache2 :: Instance().GetCoalProbBtwCfgLins2ConvertedLen(u, v, brLen, grLen, threadIndex);
    return ApproxGTPCache2 :: Instance().GetCoalProbBtwCfgLins2ConvertedLenLRU(u, v, brLen, grLen, threadIndex);
}

void ApproxGeneTreeProbHeu2 :: CalcUpperProbAt(int indexGT, int nodeST)
{
    // re-init upper at this node
    //int nlMin = GetMinLinGTLinNumAt(indexGT, nodeST);
    //for(unsigned int i=0; i<tblDPUpper[indexGT][nodeST].size(); ++i)
    //{
    //    tblDPUpper[indexGT][nodeST][i] = MAX_NEG_DOUBLE_VAL;
    //}
    
    double lenBr = treeSpecies.GetEdgeLen(nodeST);
    
    // find length grid
    int grLen = ApproxGTPCache2 :: Instance().ConvDistToGrid(lenBr);
    double lenBrUse = ApproxGTPCache2 :: Instance().ConvGridToDist(grLen);
    
    //int lenBrInt = RoundToInt(lenBr*LEN_FAC_APTP2);
    int nlMin = GetMinLinGTLinNumAt(indexGT, nodeST);
    //int nlConv = ConvDPCellIndexToLinNum(indexGT, nodeST, 0);
    int nlConv = nlMin;
    for(unsigned int nl=0; nl<tblDPUpper[indexGT][nodeST].size(); ++nl)
    {
        // enum all src lin nums on that side
        //int nlLowConv = ConvDPCellIndexToLinNum(indexGT, nodeST, nl);
        int nlLowConv = nl + nlMin;
        vector<double> listLogPs;
        //tblDPUpper[indexGT][nodeST][nl] = -1.0e100;
        for( unsigned int nlLow = nl; nlLow < tblDPLower[indexGT][nodeST].size(); ++nlLow )
        {
#if 0
            //double probCoal = CalcCoalProbFor( lenBr, nlLowConv, nlConv );
            auto it22 = mapCacheCoalProbs.find(lenBrInt);
            double logprobCoal = 0.0;
            bool fCache = false;
            if( it22 != mapCacheCoalProbs.end() )
            {
                auto it23 = it22->second.find(nlLowConv);
                if( it23 != it22->second.end() )
                {
                    auto it24 = it23->second.find(nlConv);
                    if( it24 != it23->second.end() )
                    {
                    //    logprobCoal = mapCacheCoalProbs[lenBrInt][nlLowConv][nlConv];
                        logprobCoal = it24->second;
                        fCache = true;
                    }
                }
            }
            if( fCache == false )
            {
                logprobCoal = log(CalcCoalProbFor( lenBr, nlLowConv, nlConv ));
                mapCacheCoalProbs[lenBrInt][nlLowConv][nlConv] = logprobCoal;
            }
#endif
            //double probCoalStep = CalcCoalProbFor( lenBr, nlLowConv, nlConv );
            double probCoalStep = CalcCoalProbForLenGrid(grLen, lenBrUse, nlLowConv, nlConv );
            //double probCoalStep = GetCoalProbBtwCfgLins2ConvertedLen( nlLowConv, nlConv, lenBrUse, grLen);
            //double probCoalStep = 0.000000001;
            if( probCoalStep < 1.0e-300)
            {
                probCoalStep = 1.0e-300;
            }
            double logprobCoal = log(probCoalStep);
//cout << "logprobCoal: " << logprobCoal << ", tblLower:" << tblDPLower[indexGT][nodeST][nlLow] << endl;
            double logprob = logprobCoal + tblDPLower[indexGT][nodeST][nlLow];
            //tblDPUpper[indexGT][nodeST][nl] = GetLogSumOfTwo( tblDPUpper[indexGT][nodeST][nl], logprob );
            listLogPs.push_back(logprob);
            ++nlLowConv;
        }
//cout << "ConsUpperProbAt:listLogPs: nl=" << nl << ", indexGT:" << indexGT << " " << ", nodeST: " << nodeST << " ";
//DumpDoubleVec(listLogPs);
        tblDPUpper[indexGT][nodeST][nl] = GetLogSumOfLogs(listLogPs);
        ++nlConv;
    }
}

void ApproxGeneTreeProbHeu2 :: ConsLowerProbAt(int indexGT, int nodeST)
{
    // re-init lower at this node
#if 0
    for(unsigned int i=0; i<tblDPLower[indexGT][nodeST].size(); ++i)
    {
        tblDPLower[indexGT][nodeST][i] = MAX_NEG_DOUBLE_VAL;
    }
#endif
#if 0
    // for testing
    int childLeft = treeSpecies.GetLeftDescendant(nodeST);
    int childRight = treeSpecies.GetRightDescendant(nodeST);
    int nlLeftConv = ConvDPCellIndexToLinNum(indexGT, childLeft, 0);
    vector<vector<double> > listLogPs(tblDPLower[indexGT][nodeST].size());
    //int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, ConvDPCellIndexToLinNum(indexGT, childLeft, 0) + ConvDPCellIndexToLinNum(indexGT, childRight, 0) );
    int nlTotConv = ConvDPCellIndexToLinNum(indexGT, nodeST, 0 );
    for(int nlTot=0; nlTot <(int)tblDPUpper[indexGT][nodeST].size(); ++nlTot )
    {
        //if( tblDPUpper[indexGT][childLeft][nlLeft] <= MAX_NEG_DOUBLE_VAL/2 )
        //{
        //    continue;
        //}
        int nlLeftConv = ConvDPCellIndexToLinNum(indexGT, childLeft, 0);
        for(int nlLeft = 0; nlLeft <(int)tblDPUpper[indexGT][childLeft].size();++nlLeft )
        {
            //int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, nlLeftConv +nlRightConv );
            //YW_ASSERT_INFO( nlTot <(int)tblDPLower[indexGT][nodeST].size(), "Overflow111" );
            int nlRightConv = nlTotConv-nlLeftConv;
            int nlRight = ConvLinNumToDPCellIndex(indexGT, childRight, nlRightConv );
            if( nlRight >=0 && nlRight < (int)tblDPUpper[indexGT][childRight].size() )
            {
                double logprob = tblDPUpper[indexGT][childLeft][nlLeft] + tblDPUpper[indexGT][childRight][nlRight];
                //tblDPLower[indexGT][nodeST][nlTot] = GetLogSumOfTwo(tblDPLower[indexGT][nodeST][nlTot], logprob);
                listLogPs[nlTot].push_back(logprob);
            }
            ++nlLeftConv;
        }
        ++nlTotConv;
    }
    // now finalize prob
    for(unsigned int i=0; i<listLogPs.size(); ++i)
    {
//cout << "ConsLowerProbAt:listLogPs[i]: i=" << i << ", indexGT:" << indexGT <<" " << ", nodeST: " << nodeST << " ";
//DumpDoubleVec(listLogPs[i]);
        //tblDPLower[indexGT][nodeST][i] = GetLogSumOfLogs2( listLogPs[i] );
        tblDPLower[indexGT][nodeST][i] = GetLogSumOfLogs3( listLogPs[i] );
    }
#endif

//#if 0
    
    int childLeft = treeSpecies.GetLeftDescendant(nodeST);
    int childRight = treeSpecies.GetRightDescendant(nodeST);
    int nlLeftConv = ConvDPCellIndexToLinNum(indexGT, childLeft, 0);
    
#if 0  // debug
    // new code based on convolution starts here...
    int nlRightConv = ConvDPCellIndexToLinNum(indexGT, childRight, 0);
//#if 0
    cout << "Upper left: ";
    DumpDoubleVec(tblDPUpper[indexGT][childLeft]);
    cout << "Upper right: ";
    DumpDoubleVec(tblDPUpper[indexGT][childRight]);
//#endif
    // use FFT-based convolution to speedup
    vector<double> vecLogProd;
    UtilConvoluteTwoLogVecs(tblDPUpper[indexGT][childLeft], tblDPUpper[indexGT][childRight], vecLogProd);
    
    // do a conversion
cout << "vecLogProd: ";
DumpDoubleVec(vecLogProd);
    int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, nlLeftConv + nlRightConv );
cout << "nlLeftConv: " << nlLeftConv << ", nlRightConv: " << nlRightConv << ", nlTot: " << nlTot << ", tblDPLower[indexGT][nodeST].size():" << tblDPLower[indexGT][nodeST].size() << endl;
    for(unsigned int i=0; i<vecLogProd.size(); ++i)
    {
        //YW_ASSERT_INFO( nlTot + i < (int)tblDPLower[indexGT][nodeST].size(), "Overflow444" );
        tblDPLower[indexGT][nodeST][nlTot+i] = vecLogProd[i];
    }
    
#endif
    
//#if 0
//cout << "At nodeST:" << nodeST << ", for tree " << indexGT << ", convolution size: left: " << (int)tblDPUpper[indexGT][childLeft].size() << ", right: " << (int)tblDPUpper[indexGT][childRight].size() << endl;
    //vector<vector<double> > listLogPs(tblDPLower[indexGT][nodeST].size());
    //int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, ConvDPCellIndexToLinNum(indexGT, childLeft, 0) + ConvDPCellIndexToLinNum(indexGT, childRight, 0) );

    int nlLeftMin = GetMinLinGTLinNumAt(indexGT, childLeft);
    //int nlLeftMax  = GetMaxLinGTLinNumAt(indexGT, childLeft);
    int nlRightMin = GetMinLinGTLinNumAt(indexGT, childRight);
    //int nlRightMax  = GetMaxLinGTLinNumAt(indexGT, childRight);
    int nlRightMax  = listMaxLinGTLinNumAt[indexGT][childRight];
//cout << "nlLeftMiin:" << GetMinLinGTLinNumAt(indexGT, childLeft) << ", nlLeftMax:" << GetMaxLinGTLinNumAt(indexGT, childLeft) << ", nlRightMin:" << nlRightMin << ", nlRightMax:" << nlRightMax << ", tot min:" << GetMinLinGTLinNumAt(indexGT, nodeST) << ", tot max:" << GetMaxLinGTLinNumAt(indexGT, nodeST) << endl;
    
    int nlTotConv = ConvDPCellIndexToLinNum(indexGT, nodeST, 0 );
    for(int nlTot=0; nlTot <(int)tblDPUpper[indexGT][nodeST].size(); ++nlTot )
    {
//cout << "nlTot:" << nlTot << endl;
        if(nlTotConv < nlLeftMin + nlRightMin )
        {
            //cout << "Skipped\n";
            ++nlTotConv;
            continue;
        }
        
        vector<double> listLogPs;
        //if( tblDPUpper[indexGT][childLeft][nlLeft] <= MAX_NEG_DOUBLE_VAL/2 )
        //{
        //    continue;
        //}
        //int nLeft0 = std::max(0, ConvLinNumToDPCellIndex(indexGT, childLeft, nlTotConv - nlRightMax) );
        int nLeft0 = std::max( (nlTotConv - nlRightMax) - nlLeftMin, 0 );
        //int nlLeftConv = ConvDPCellIndexToLinNum(indexGT, childLeft, nLeft0);
        int nlLeftConv = nlLeftMin + nLeft0;
//cout << "nLeft0:" << nLeft0 << ", left table size: " <<  (int)tblDPUpper[indexGT][childLeft].size() << ", nlTotConv:" << nlTotConv << ", nlLeftConv: " << nlLeftConv << endl;
        for(int nlLeft = nLeft0; nlLeft <(int)tblDPUpper[indexGT][childLeft].size() && nlTotConv-nlLeftConv >= nlRightMin; ++nlLeft )
        {
            //int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, nlLeftConv + nlRightConv );
            //YW_ASSERT_INFO( nlTot <(int)tblDPLower[indexGT][nodeST].size(), "Overflow111" );
            int nlRightConv = nlTotConv-nlLeftConv;
            //int nlRight = ConvLinNumToDPCellIndex(indexGT, childRight, nlRightConv );
            int nlRight = nlRightConv - nlRightMin;
            
//cout << "nlRight:" << nlRight << ", nlRightConv:" << nlRightConv  << ", right bound:" << (int)tblDPUpper[indexGT][childRight].size() << ", nlLeft:" << nlLeft << ", nlLeftConv:" << nlLeftConv << ", nlTot:" << nlTot << ", nlTotConv:" << nlTotConv << endl;
            
            //YW_ASSERT_INFO(nlRight >=0 && nlRight < (int)tblDPUpper[indexGT][childRight].size(), "WRONG111");
            //if( nlRight >=0 && nlRight < (int)tblDPUpper[indexGT][childRight].size() )
            //{
                double logprob = tblDPUpper[indexGT][childLeft][nlLeft] + tblDPUpper[indexGT][childRight][nlRight];
                //tblDPLower[indexGT][nodeST][nlTot] = GetLogSumOfTwo(tblDPLower[indexGT][nodeST][nlTot], logprob);
                //listLogPs[nlTot].push_back(logprob);
                listLogPs.push_back(logprob);
            //}
            //if( nlRight < 0 )
            //{
                // no need to continue
            //    break;
            //}
            ++nlLeftConv;
        }
        ++nlTotConv;
        
        tblDPLower[indexGT][nodeST][nlTot] = GetLogSumOfLogs3( listLogPs );
    }
//cout << "tblDPLower at this node:";
//DumpDoubleVec(tblDPLower[indexGT][nodeST]);
    
#if 0
    for(int nlLeft = 0; nlLeft <(int)tblDPUpper[indexGT][childLeft].size(); ++nlLeft )
    {
        //if( tblDPUpper[indexGT][childLeft][nlLeft] <= MAX_NEG_DOUBLE_VAL/2 )
        //{
        //    continue;
        //}
        int nlRightConv = ConvDPCellIndexToLinNum(indexGT, childRight, 0);
        int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, nlLeftConv + nlRightConv );
        for(int nlRight = 0; nlRight <(int)tblDPUpper[indexGT][childRight].size(); ++nlRight )
        {
            //int nlTot = ConvLinNumToDPCellIndex(indexGT, nodeST, nlLeftConv + nlRightConv );
            //YW_ASSERT_INFO( nlTot <(int)tblDPLower[indexGT][nodeST].size(), "Overflow111" );
            double logprob = tblDPUpper[indexGT][childLeft][nlLeft] + tblDPUpper[indexGT][childRight][nlRight];
            //tblDPLower[indexGT][nodeST][nlTot] = GetLogSumOfTwo(tblDPLower[indexGT][nodeST][nlTot], logprob);
            listLogPs[nlTot].push_back(logprob);
            ++nlRightConv;
            ++nlTot;
        }
        ++nlLeftConv;
    }
#endif
    // now finalize prob
    //for(unsigned int i=0; i<listLogPs.size(); ++i)
    //{
//cout << "ConsLowerProbAt:listLogPs[i]: i=" << i << ", indexGT:" << indexGT << " " << ", nodeST: " << nodeST << " ";
//DumpDoubleVec(listLogPs[i]);
        //tblDPLower[indexGT][nodeST][i] = GetLogSumOfLogs2( listLogPs[i] );
    //    tblDPLower[indexGT][nodeST][i] = GetLogSumOfLogs3( listLogPs[i] );
    //}
//#endif
//#endif
}

double ApproxGeneTreeProbHeu2 :: CalcCoalCoeff(int ub, int vb)
{
    //YW_ASSERT_INFO( ub >= vb, "ApproxGeneTreeProbHeu2 :: CalcCoalCoeff: ub must be at least as large as vb" );
    int cb = ub-vb;
    //double res = 1.0/CalcdbVal(ub, cb);
    double res = 1.0/ApproxGTPCache2 :: Instance().GetdbValCache(ub, cb);
    //res *= ApproxGTPCache2 :: Instance().GetwbValCache( ub, cb);
    return res;
}

int ApproxGeneTreeProbHeu2 :: GetNumDPCellsAt(int indexGT, int nodeST) const
{
    //return GetMaxLinGTLinNumAt(indexGT, nodeST) - GetMinLinGTLinNumAt(indexGT, nodeST)+1;
    return this->listMaxLinGTLinNumAt[indexGT][nodeST] - GetMinLinGTLinNumAt(indexGT, nodeST)+1;
}
int ApproxGeneTreeProbHeu2 :: ConvDPCellIndexToLinNum(int indexGT, int nodeST, int nc) const
{
    return nc + GetMinLinGTLinNumAt(indexGT, nodeST);
}
int ApproxGeneTreeProbHeu2 :: ConvLinNumToDPCellIndex(int indexGT, int nodeST, int nl) const
{
    return nl - GetMinLinGTLinNumAt(indexGT, nodeST);
}

double ApproxGeneTreeProbHeu2 :: GetCoalProbBtwCfgLins2ConvertedLen( int u, int v, double brLen, int grLen )
{
    // don't use cache
    //if( u <v || v == 0 )
    //{
    //    cout << "ApproxGeneTreeProbHeu2 :: GetCoalProbBtwCfgLins2ConvertedLen: u=" << u << ", v=" << v << endl;
    //}
    YW_ASSERT_INFO( u >= v && v > 0, "FATAL ERROR2345");
    
    pair<pair<int,int>,int> pp(std::make_pair(u,v), grLen);
    double res = 1.0;
    auto it = listCachedCoalProb.find(pp);
    if( it != listCachedCoalProb.end() )
    {
        res = it->second;
        return res;
    }
    
    // claculate it first
    int cb = u-v;
    //double res = 1.0;
    double dval = ApproxGTPCache :: Instance().GetProdChoose2(u, cb);
    res = 1.0/dval;
//cout << "** coeff after dbVal: " << res << endl;
//#endif
    //res *= GetwbValCache( u, cb);
//cout << "coeff after wbVal: " << res << endl;
    
    // get branch if still worth it (value not too small)
    if( res >= 1.0e-300)
    {
//cout << "u:" << u << ",v:" << v << ", brLen:" << brLen << endl;
        //res *= UtilsCalcCoalProbBranch2(u, v, brLen);
        res *= UtilsCalcCoalProbBranch2(u, v, brLen);
    }
//cout << "coeff after br prob: " << res << endl;
    listCachedCoalProb[pp] = res;
    return res;
}


//////////////////////////////////////////////////////////////////////////////////
// Analyze common subtrees in species trees

STSubtreeDepot& STSubtreeDepot :: Instance()
{
    static STSubtreeDepot inst;
    return inst;
}

int STSubtreeDepot :: GetSTIdFor(const string &subtreeNW)
{
    auto it = mapNWtoId.find(subtreeNW);
    if( it != mapNWtoId.end())
    {
        return it->second;
    }
    ++idNextToUse;
    mapNWtoId[subtreeNW] = idNextToUse;
    return idNextToUse;
}

void STSubtreeDepot :: DumpStats() const
{
    cout << "STSubtreeDepot: number of species subtrees: " << mapNWtoId.size() << endl;
}

STSubtreeDepot :: STSubtreeDepot() : idNextToUse(0)
{
}

