//
//  AGFastLenOpt.cpp
//  
//
//  Created by Yufeng Wu on 3/18/24.
//

#include "AGFastLenOpt.hpp"
#include "AGGeneTreeProb.hpp"
#include "Utils4.h"
#include <thread>
#include <mutex>
#include <cmath>
using namespace std;

//static const int DEF_SIZE_BR_LEN_GRID = 44;
//static const double LIST_DEF_BR_LEN_GRID_VALS[DEF_SIZE_BR_LEN_GRID] = {0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0};
static const int DEF_SIZE_BR_LEN_GRID = 35;
static const double LIST_DEF_BR_LEN_GRID_VALS[DEF_SIZE_BR_LEN_GRID] = {0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 3.0, 5.0};
//static const int DEF_SIZE_BR_LEN_GRID = 20;
//static const double LIST_DEF_BR_LEN_GRID_VALS[DEF_SIZE_BR_LEN_GRID] = {0.0001, 0.0005, 0.001,0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.3, 0.5};
static vector<double> listDefSizesAGLO(LIST_DEF_BR_LEN_GRID_VALS, LIST_DEF_BR_LEN_GRID_VALS+sizeof(LIST_DEF_BR_LEN_GRID_VALS)/sizeof(LIST_DEF_BR_LEN_GRID_VALS[0]));
//static const int DEF_SIZE_MR_GRID = 13;
//static const double LIST_DEF_MR_GRID_VALS[DEF_SIZE_MR_GRID] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
static const int DEF_SIZE_MR_GRID = 11;
static const double LIST_DEF_MR_GRID_VALS[DEF_SIZE_MR_GRID] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
//static const int DEF_SIZE_MR_GRID = 7;
//static const double LIST_DEF_MR_GRID_VALS[DEF_SIZE_MR_GRID] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99};

//***********************************************************************************
// Length grid

AGGTBranchLengthGrid & AGGTBranchLengthGrid :: Instance()
{
    static AGGTBranchLengthGrid instance;
    return instance;
}

AGGTBranchLengthGrid :: AGGTBranchLengthGrid()
{
    Init();
}

double AGGTBranchLengthGrid :: SnapToGrid(double brOrig, int *posCloseOut) const
{
// no longer snap branch length
//return brOrig;
    int posClose = GetClosestTo(listBrLenGridVals, brOrig);
    //int posClose = GetClosestTo2(listBrLenGridVals, brOrig);
    YW_ASSERT_INFO( posClose >= 0, "Fail to find" );
    if( posCloseOut != NULL)
    {
        *posCloseOut = posClose;
    }
    return listBrLenGridVals[posClose];
}

double AGGTBranchLengthGrid :: ReducePercison(double brOrig) const
{
    // only keep 4 digits after demical point
    const int FAC = 10000;
    return ((double)((int)(FAC*brOrig+0.5)))/FAC;
}

double AGGTBranchLengthGrid :: SnapMRToGrid(double mr) const
{
//return mr;
    int posClose = GetClosestTo(listMRGridVals, mr);
    YW_ASSERT_INFO( posClose >= 0, "Fail to find" );
    return listMRGridVals[posClose];
}

void AGGTBranchLengthGrid :: Init()
{
    listBrLenGridVals.clear();
    for(int i=0; i<DEF_SIZE_BR_LEN_GRID; ++i)
    {
        listBrLenGridVals.push_back(LIST_DEF_BR_LEN_GRID_VALS[i]);
    }
    listMRGridVals.clear();
    for(int i=0; i<DEF_SIZE_MR_GRID; ++i)
    {
        listMRGridVals.push_back( LIST_DEF_MR_GRID_VALS[i] );
    }
}



//***********************************************************************************
// Optimize branch lengths of a given network topology wrt gene trees

//int maxNumOptRound = 1;     // test 2 rounds opt


AGFastLenOpt :: AGFastLenOpt(GenealogicalNetwork &netIn, vector<PhylogenyTreeBasic *> &listGeneTreesIn, int threadIndIn) : network(netIn), listGeneTrees(listGeneTreesIn), pBranchInNetCur(NULL), pMixNodeCur(NULL), logprobPre( MAX_NEG_DOUBLE_VAL ), maxNumOptRound(1), threadIndex(threadIndIn)
{
    //pProbCalc = new AGGeneTreeQuickCoal2(network);
    //pProbCalc2 = new AGGeneTreeProb(network, listGeneTrees);
    pProbCalc2 = new AGGeneTreeProb2(network, listGeneTrees, threadIndIn);
}

AGFastLenOpt:: ~AGFastLenOpt()
{
    //delete pProbCalc;
    delete pProbCalc2;
}

/*
void AGFastLenOpt :: SetMultithread(int numT)
{
    numThreads = numT;
    this->pProbCalc2->SetMultithread(numThreads);
} */

double AGFastLenOpt :: CalcProb()
{
    /*
    if( this->numThreads > 1 )
    {
        return CalcProbMT();
    } */
    double res = pProbCalc2->CalcProb();
    //double res = 0.0;
    //for(unsigned int i=0; i<listGeneTrees.size(); ++i)
    //{
    //    res += pProbCalc->CalcGeneTreeProbHeu( *listGeneTrees[i] );
    //}
//cout << "---- AGFastLenOpt: prob: " << res << ", for network: \n";
//this->network.DumpMargTrees(false);
    return res;
}

double AGFastLenOpt :: Optimize()
{
//return Optimize2();
    // snap all edge length
    SnapBrLens(this->network);
#if 0
cout << "Before optimize, network: \n";
//this->network.DumpMargTrees(false);
this->network.Dump();
cout << "Set of critical edges: \n";
for(auto x: setBranchesCritical)
{
x->Dump();
cout << endl;
}
#endif
    // optimize over the branch lengths by iterate through all branches
    //set<GenealogicalNetworkBranch *> setAllBranchesInNet;
    //network.GetAllBranches( setAllBranchesInNet );
    vector<GenealogicalNetworkBranch *> vecAllBranchesInNet;
    //network.GetAllBranchesOrdered(vecAllBranchesInNet);
    network.GetAllBranchesOrdered2(vecAllBranchesInNet);
    
    if( setBranchesCritical.size() > 0 )
    {
        // only use these critical edges
        vector<GenealogicalNetworkBranch *> vecAllBranchesInNet2;
        for(unsigned int i=0; i<vecAllBranchesInNet.size(); ++i)
        {
            GenealogicalNetworkBranch *pb = vecAllBranchesInNet[i];
            if( setBranchesCritical.find(pb) != setBranchesCritical.end() )
            {
                vecAllBranchesInNet2.push_back(pb);
            }
        }
        vecAllBranchesInNet = vecAllBranchesInNet2;
    }
    
    //double loglikeliBest = MAX_NEG_DOUBLE_VAL;
    double loglikeliBest = this->logprobPre;
    
//cout << "^^^^^Optimizing branch length: \n";
    int round = 0;
    while(true)
    {
//cout << "Current network: ";
//network.DumpMargTrees(false);
        //
        double loglikeliStep = MAX_NEG_DOUBLE_VAL;
        for( vector<GenealogicalNetworkBranch *>::iterator it= vecAllBranchesInNet.begin(); it != vecAllBranchesInNet.end(); ++it )
        {
            //
            pBranchInNetCur = *it;
            
            if(pBranchInNetCur == NULL)
            {
                
                // YW: don't optimize over root branch
                continue;
            }
            
            if( setBranchesSkip.find(pBranchInNetCur) != setBranchesSkip.end() )
            {
                continue;
            }
            
            
            pMixNodeCur = NULL;
            double brCurStep = pBranchInNetCur->GetLength();
            double brCurNew = 0.0;
            double brMin = GetMinBrLen(round, brCurStep);
            double brMax = GetMaxBrLen(round, brCurStep);
            double tolBr = GetTolBrLen(round);
//cout << "BRENT (length): [" << brMin << "," << brMax << "]: tol: " << tolBr << endl;
            setProcessedLens.clear();
            setProcessedMRs.clear();
            double loglikeliStepStep = -1.0*Func1DMinBrent( brMin, brCurStep, brMax, tolBr,  &brCurNew);
            //setProcessedLens.clear();
            //setProcessedMRs.clear();
//cout << "Orig len: " << brCurStep << ", new len: " << brCurNew;
            brCurNew = AGGTBranchLengthGrid::Instance().SnapToGrid(brCurNew);
//cout << "After snapping, new length: " << brCurNew << endl;
//cout << "Orig len: " << brCurStep << ", new len: " << brCurNew;
//cout << "one round of Brent: best log-likelihood: " << loglikeliStepStep << ", current best: " <<  loglikeliStep << ", new len: " << brCurNew << endl;
//cout << "   current branch: ";
//pBranchInNetCur->Dump();
//cout << endl;
            
            if( IsSignificantlyLarge( loglikeliStepStep, loglikeliStep ) == false )
            {
                // this branch does not change much
                pBranchInNetCur->SetLength( brCurStep );
                setBranchesSkip.insert( pBranchInNetCur );
            }
            else
            {
                pBranchInNetCur->SetLength( brCurNew );
                loglikeliStep = loglikeliStepStep;
            }
        }
        
//cout << "After one round of branch length search, loglikeliStep: " << loglikeliStep << endl;
        // now search over all mixing node to adjust mixing coefficients
        vector<GenealogicalNetworkNode *> listMixNodes;
        network.GetMixNodes( listMixNodes );
        for(int i=0; i<(int)listMixNodes.size(); ++i)
        {
            pMixNodeCur = listMixNodes[i];
            if( setMixNodesSkip.find( pMixNodeCur) != setMixNodesSkip.end() )
            {
                continue;
            }
            
            pBranchInNetCur = NULL;
            double mrCur = pMixNodeCur->GetMixRatio();
            double mrBest = mrCur;
            double mrMin = GetMinMR(round, mrCur);
            double mrMax = GetMaxMR(round, mrCur);
            double tolMr = GetTolMR(round);
//cout << "BRENT (mix ratio): [" << mrMin << "," << mrMax << "]: tol: " << tolMr << endl;
            double loglikeliStepStep2 = -1.0*Func1DMinBrent(mrMin, mrCur, mrMax, tolMr, &mrBest);
//cout << "Orig mr: " << mrCur << ", new mr: " << mrBest;
            mrBest = AGGTBranchLengthGrid::Instance().SnapMRToGrid(mrBest);
//cout << "After snap, new mr: " << mrBest << endl;
//cout << "Orig mr: " << mrCur << ", new mr: " << mrBest;
//cout << "   one round of Brent (for mixiing estimate): best log-likelihood: " << loglikeliStepStep2 << endl;
            
            if( IsSignificantlyLarge( loglikeliStepStep2, loglikeliStep ) == false )
            {
                // this branch does not change much
                pMixNodeCur->SetMixRatio(mrCur);
                // forbid more search if we only do one round
                //if( maxNumOptRound <= 1 )
                //{
                //    setMixNodesSkip.insert( pMixNodeCur );
                //}
            }
            else
            {
                pMixNodeCur->SetMixRatio(mrBest);
                loglikeliStep = loglikeliStepStep2;
            }
        }
        
//cout << "--round " << round << "-- At one round of branch length search: prev lokelihood: " << loglikeliBest << ", improved likelihood: " << loglikeliStep << endl;
        // stop if the likelihood value does not improve significantly
        if( IsSignificantlyLarge( loglikeliStep, loglikeliBest ) == false )
        {
            loglikeliBest = loglikeliStep;
            break;
        }
        else
        {
            loglikeliBest = loglikeliStep;
        }
        ++round;
        
        if( round >= maxNumOptRound )
        {
            break;
        }
    }
//cout << "Done length opt: loglikelihood: " << loglikeliBest << endl;
// make sure the prob matches
//double probRecalc = CalcProb();
//cout << "VERIFYING prob: " << probRecalc << endl;
    
    return loglikeliBest;
}

double AGFastLenOpt :: Optimize2()
{
    // snap all edge length
    SnapBrLens(this->network);
#if 0
cout << "Before optimize, network: \n";
//this->network.DumpMargTrees(false);
this->network.Dump();
cout << "Set of critical edges: \n";
for(auto x: setBranchesCritical)
{
x->Dump();
cout << endl;
}
#endif
    // optimize over the branch lengths by iterate through all branches
    //set<GenealogicalNetworkBranch *> setAllBranchesInNet;
    //network.GetAllBranches( setAllBranchesInNet );
    vector<GenealogicalNetworkBranch *> vecAllBranchesInNet;
    //network.GetAllBranchesOrdered(vecAllBranchesInNet);
    network.GetAllBranchesOrdered2(vecAllBranchesInNet);
    
    if( setBranchesCritical.size() > 0 )
    {
        // only use these critical edges
        vector<GenealogicalNetworkBranch *> vecAllBranchesInNet2;
        for(unsigned int i=0; i<vecAllBranchesInNet.size(); ++i)
        {
            GenealogicalNetworkBranch *pb = vecAllBranchesInNet[i];
            if( setBranchesCritical.find(pb) != setBranchesCritical.end() )
            {
                vecAllBranchesInNet2.push_back(pb);
            }
        }
        vecAllBranchesInNet = vecAllBranchesInNet2;
    }
    
    //double loglikeliBest = MAX_NEG_DOUBLE_VAL;
    double loglikeliBest = this->logprobPre;
    
//cout << "^^^^^Optimizing branch length: \n";
    int round = 0;
    while(true)
    {
//cout << "Current network: ";
//network.DumpMargTrees(false);
        //
        //double loglikeliStep = MAX_NEG_DOUBLE_VAL;
        double loglikeliStep = loglikeliBest;
        for( vector<GenealogicalNetworkBranch *>::iterator it= vecAllBranchesInNet.begin(); it != vecAllBranchesInNet.end(); ++it )
        {
            //
            pBranchInNetCur = *it;
            
            if(pBranchInNetCur == NULL)
            {
                // YW: don't optimize over root branch
                continue;
            }
            
            if( setBranchesSkip.find(pBranchInNetCur) != setBranchesSkip.end() )
            {
                continue;
            }
            
            pMixNodeCur = NULL;
            double brCurStep = pBranchInNetCur->GetLength();
            double brCurNew = pBranchInNetCur->GetLength();
            double loglikeliStepStep = OptimizeBranch(pBranchInNetCur, loglikeliStep, brCurNew);
            
            //setProcessedLens.clear();
            //setProcessedMRs.clear();
//cout << "Orig len: " << brCurStep << ", new len: " << brCurNew;
            brCurNew = AGGTBranchLengthGrid::Instance().SnapToGrid(brCurNew);
//cout << "After snapping, new length: " << brCurNew << endl;
//cout << "Orig len: " << brCurStep << ", new len: " << brCurNew;
//cout << "one round of Brent: best log-likelihood: " << loglikeliStepStep << ", current best: " <<  loglikeliStep << ", new len: " << brCurNew << endl;
//cout << "   current branch: ";
//pBranchInNetCur->Dump();
//cout << endl;
            
            if( IsSignificantlyLarge( loglikeliStepStep, loglikeliStep ) == false )
            {
                // this branch does not change much
                pBranchInNetCur->SetLength( brCurStep );
                setBranchesSkip.insert( pBranchInNetCur );
            }
            else
            {
                pBranchInNetCur->SetLength( brCurNew );
                loglikeliStep = loglikeliStepStep;
            }
        }
        // now search over all mixing node to adjust mixing coefficients
        vector<GenealogicalNetworkNode *> listMixNodes;
        network.GetMixNodes( listMixNodes );
        for(int i=0; i<(int)listMixNodes.size(); ++i)
        {
            pMixNodeCur = listMixNodes[i];
            if( setMixNodesSkip.find( pMixNodeCur) != setMixNodesSkip.end() )
            {
                continue;
            }
            
            pBranchInNetCur = NULL;
            double mrCur = pMixNodeCur->GetMixRatio();
            double mrBest = mrCur;
            double mrMin = GetMinMR(round, mrCur);
            double mrMax = GetMaxMR(round, mrCur);
            double tolMr = GetTolMR(round);
//cout << "BRENT (mix ratio): [" << mrMin << "," << mrMax << "]: tol: " << tolMr << endl;
            double loglikeliStepStep2 = -1.0*Func1DMinBrent(mrMin, mrCur, mrMax, tolMr, &mrBest);
//cout << "Orig mr: " << mrCur << ", new mr: " << mrBest;
            mrBest = AGGTBranchLengthGrid::Instance().SnapMRToGrid(mrBest);
//cout << "After snap, new mr: " << mrBest << endl;
//cout << "Orig mr: " << mrCur << ", new mr: " << mrBest;
//cout << "   one round of Brent (for mixiing estimate): best log-likelihood: " << loglikeliStepStep2 << endl;
            
            if( IsSignificantlyLarge( loglikeliStepStep2, loglikeliStep ) == false )
            {
                // this branch does not change much
                pMixNodeCur->SetMixRatio(mrCur);
                setMixNodesSkip.insert( pMixNodeCur );
            }
            else
            {
                pMixNodeCur->SetMixRatio(mrBest);
                loglikeliStep = loglikeliStepStep2;
            }
        }
        
//cout << "--round " << round << "-- At one round of branch length search: prev lokelihood: " << loglikeliBest << ", improved likelihood: " << loglikeliStep << endl;
        // stop if the likelihood value does not improve significantly
        if( IsSignificantlyLarge( loglikeliStep, loglikeliBest ) == false )
        {
            loglikeliBest = loglikeliStep;
            break;
        }
        else
        {
            loglikeliBest = loglikeliStep;
        }
        ++round;
        
        if( round >= maxNumOptRound )
        {
            break;
        }
    }
//cout << "Done length opt: loglikelihood: " << loglikeliBest << endl;
// make sure the prob matches
//double probRecalc = CalcProb();
//cout << "VERIFYING prob: " << probRecalc << endl;
    
    return loglikeliBest;
}

double AGFastLenOpt :: OptimizeBranch(GenealogicalNetworkBranch *brCurr, double logprobPre, double &brLenBest)
{
    // use golden section optimizaion
    const double invphi = (std::sqrt(5) - 1) / 2; // phi
    const double invphi2 = (3 - std::sqrt(5)) / 2;  // 1 / phi^2
    const double tol = 0.2;
    
    //
    int a = 0, b = (int)listDefSizesAGLO.size()-1;
    double h = b-a;
    int n = (int)(std::ceil(std::log(tol / h) / std::log(invphi)));
//cout << "Golden section: number of iterations: " << n << endl;
    int c = (int)(a + invphi2 * h);
    int d = (int)(a + invphi * h);
    UpdateProbComputeForBranch( brCurr, listDefSizesAGLO[c] );
    double yc = CalcProb();
    UpdateProbComputeForBranch( brCurr, listDefSizesAGLO[d] );
    double yd = CalcProb();
    bool fBetterFound = IsSignificantlyLarge(yc, logprobPre) || IsSignificantlyLarge(yd, logprobPre);
    int numFailedEvals = 0;
    const int MAX_NUM_NO_INC_TRY = 3;

    for(int i=0; i<n; ++i  )
    {
//cout << "a:" << a << ", b:" << b << ", c:" << c << ", d:" << d << ", yc:" << yc << ", yd:" << yd << ", h:" << h << ", logprobPre:" << logprobPre << endl;
        if (yc > yd)  // yc > yd to find the maximum
        {
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = (int)(a + invphi2 * h);
            //double ycOld = yc;
            UpdateProbComputeForBranch( brCurr, listDefSizesAGLO[c] );
            yc = CalcProb();
            if( IsSignificantlyLarge(yc, logprobPre) == false )
            {
                if( fBetterFound ==false &&  ++numFailedEvals > MAX_NUM_NO_INC_TRY )
                {
                    //break;
                }
            }
            else
            {
                fBetterFound = true;
            }
        }
        else
        {
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = (int)(a + invphi * h);
            UpdateProbComputeForBranch( brCurr, listDefSizesAGLO[d] );
            //double ydOld = yd;
            yd = CalcProb();
            if( IsSignificantlyLarge(yd, logprobPre) == false )
            {
                if( fBetterFound == false && ++numFailedEvals > MAX_NUM_NO_INC_TRY )
                {
                    //break;
                }
            }
            else
            {
                fBetterFound = true;
            }
        }
    }
    if( yc > yd )
    {
        brLenBest = listDefSizesAGLO[(a+d)/2];
    }
    else
    {
        brLenBest = listDefSizesAGLO[(b+c)/2];
    }
    return std::max(yc, yd);
}

double AGFastLenOpt :: EvaluateAt(double pt, void *pParam)
{
#if 0
cout << "___________Now updating branch length to: " << pt << " for branch: ";
if( pBranchInNetCur != NULL )
pBranchInNetCur->Dump();
else
{
cout << " ROOT branch";
}
cout << endl;
#endif
    //
    double lenUse = AGGTBranchLengthGrid::Instance().SnapToGrid(pt);
    // 06/13/24: no longer snap br length
    //double lenUse = pt;
    //double lenUse = AGGTBranchLengthGrid::Instance().ReducePercison(pt);
    if( pBranchInNetCur != NULL )
    {
        if( setProcessedLens.find(lenUse) != setProcessedLens.end() )
        {
            return setProcessedLens[lenUse];
        }
//cout << "EvaluateAt: lenUse:" << lenUse << endl;
        UpdateProbComputeForBranch( pBranchInNetCur, lenUse );
    }
    else
    {
        double mrUse = lenUse;
        if( setProcessedMRs.find(mrUse) != setProcessedMRs.end() )
        {
            return setProcessedMRs[mrUse];
        }
        UpdateProbComputeForMR( pMixNodeCur, mrUse );
    }
    double prob = -1.0*CalcProb();
    if( pBranchInNetCur != NULL )
    {
        setProcessedLens[lenUse] = prob;
    }
    else
    {
        setProcessedMRs[lenUse] = prob;
    }
//cout << "AGFastLenOpt :: EvaluateAt: prob: " << prob << endl;
    return prob;
}

void AGFastLenOpt :: UpdateProbComputeForBranch( GenealogicalNetworkBranch *pBrToChange, double brNew )
{
    // find the affected marginal trees
    double lenOld = pBrToChange->GetLength();
    
    if( std::fabs(brNew - lenOld) < 0.00000001 )
    {
        return;
    }
    
    pBrToChange->SetLength(brNew);
    //pprobCompute->UpdateBranch( pBrToChange );
    pProbCalc2->UpdateBranch( pBrToChange, brNew - lenOld );
}

void AGFastLenOpt :: UpdateProbComputeForMR( GenealogicalNetworkNode *pNodeMRChange, double mrNew )
{
    YW_ASSERT_INFO( pNodeMRChange != NULL, "Fail" );
    //
    pNodeMRChange->SetMixRatio( mrNew );
    pProbCalc2->UpdateWts();
}

#if 0
// use mutex to
//static void UtilCalcProbMT( AGGeneTreeQuickCoal2 * pProbCalc, int numTotTrees, vector<bool> *pvecProcessed, vector<double> *pListProb, vector<PhylogenyTreeBasic *> *plistGeneTrees )
static void UtilCalcProbMT( AGGeneTreeProb * pProbCalc2, int numTotTrees, vector<bool> *pvecProcessed, vector<double> *pListProb, vector<PhylogenyTreeBasic *> *plistGeneTrees )
{
    int indTree = 0;
    static std::mutex mut;
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
        //
        double pStep = pProbCalc->CalcGeneTreeProbHeu( *(*plistGeneTrees)[indTree] );
        
        // update prune ST result
        (*pListProb)[indTree] = pStep;
    }
}
#endif

/*
double AGFastLenOpt :: CalcProbMT()
{
//cout << "Num of threads: " << numThreads << endl;
    this->pProbCalc2->SetMultithread(numThreads);
    return this->pProbCalc2->CalcProb();
#if 0
    // do multi-threading
//cout << "Start multi-threading...\n";
    int numThreadsUse = numThreads;
    if( (int)listGeneTrees.size() < numThreadsUse )
    {
        numThreadsUse = listGeneTrees.size();
    }
    vector<bool> vecProcessed(listGeneTrees.size());
    for(int i=0; i<(int)listGeneTrees.size(); ++i)
    {
        vecProcessed[i] = false;
    }
    vector<double> listProbs(listGeneTrees.size());
    for(int k=0; k<listGeneTrees.size(); ++k)
    {
        listProbs[k] = 0.0;
    }
    // vector for a node being processed
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        thread *pthr = new thread(UtilCalcProbMT, this->pProbCalc, listGeneTrees.size(), &vecProcessed, &listProbs, &this->listGeneTrees );
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
    
    // consider all rSPR move
    double probSum = 0.0;
    
    // statt from each node
    for( int k = 0; k < (int)listProbs.size(); ++k )
    {
        probSum += listProbs[k];
//cout << "--- max Q value over ALL sites: " << probMaxStep << endl;
    }
    return probSum;
#endif
}
*/

const double MIN_BR_LEN = 0.001;
//const double MAX_BR_LEN = 1.0;
const double MAX_BR_LEN = 0.5;
//const double TOL_BR_LEN = 0.001;
const double TOL_BR_LEN_INIT = 0.01;
//const double TOL_BR_LEN_INIT = 0.001;
//const double TOL_BR_LEN_MIN = 0.0001;
const double TOL_BR_LEN_MIN = 0.001;
//const double MIN_MR = 0.01;
const double MIN_MR = 0.001;
const double MAX_MR = 0.99;
//const double MIN_MR = 0.4;
//const double MAX_MR = 0.6;
//const double TOL_MR_MIN = 0.01;
const double TOL_MR_MIN = 0.001;
//const double TOL_MR_INIT = 0.05;
//const double TOL_MR_INIT = 0.01;
const double TOL_MR_INIT = 0.001;

static double GNLOFactor(int round)
{
    const double facInit = 8.0;
    double res = facInit/( 0x1 << round );
    if( res < 2.0 )
    {
        res = 2.0;
    }
    return res;
}

double AGFastLenOpt :: GetMinBrLen(int round, double brLenCur) const
{
    if( round == 0 )
    {
        return MIN_BR_LEN;
    }
    else
    {
        double fac=GNLOFactor(round);
        double len = brLenCur/fac;
        if( len < MIN_BR_LEN)
        {
            len = MIN_BR_LEN;
        }
        return len;
    }
}
double AGFastLenOpt :: GetMaxBrLen(int round, double brLenCur) const
{
    if( round == 0 )
    {
        return MAX_BR_LEN;
    }
    else
    {
        double fac=GNLOFactor(round);
        double len = brLenCur*fac;
        if( len > MAX_BR_LEN)
        {
            len = MAX_BR_LEN;
        }
        return len;
    }
}
double AGFastLenOpt :: GetTolBrLen(int round) const
{
    if( round == 0)
    {
        return TOL_BR_LEN_INIT;
    }
    else
    {
        double tol = TOL_BR_LEN_INIT;
        tol = tol/(0x1 << round );
        if( tol < TOL_BR_LEN_MIN)
        {
            tol = TOL_BR_LEN_MIN;
        }
        return tol;
    }
}
double AGFastLenOpt :: GetMinMR(int round, double mrCur) const
{
//    return MIN_MR;
//#if 0
    if( round == 0 )
    {
        return MIN_MR;
    }
    else
    {
        double fac=GNLOFactor(round);
        double len = mrCur/fac;
        if( len < MIN_MR)
        {
            len = MIN_MR;
        }
        return len;
    }
//#endif
}
double AGFastLenOpt :: GetMaxMR(int round, double mrCur) const
{
//    return MAX_MR;
//#if 0
    if( round == 0 )
    {
        return MAX_MR;
    }
    else
    {
        double fac=GNLOFactor(round);
        double len = mrCur*fac;
        if( len > MAX_MR)
        {
            len = MAX_MR;
        }
        return len;
    }
//#endif
}
double AGFastLenOpt :: GetTolMR(int round) const
{
    //return TOL_MR_MIN;
//#if 0
    if( round == 0)
    {
        return TOL_MR_INIT;
    }
    else
    {
        double tol = TOL_MR_INIT;
        tol = tol/(0x1 << round );
        if( tol < TOL_MR_MIN)
        {
            tol = TOL_MR_MIN;
        }
        return tol;
    }
//#endif
}


void AGFastLenOpt :: SnapBrLens(GenealogicalNetwork &networkCurrent)
{
// 06/13/24: no longer snap br lengths
//return;
    //
    set<GenealogicalNetworkBranch *> setAllBranchesInNet;
    networkCurrent.GetAllBranches( setAllBranchesInNet );
    for( set<GenealogicalNetworkBranch *> :: iterator it = setAllBranchesInNet.begin(); it != setAllBranchesInNet.end(); ++it )
    {
        GenealogicalNetworkBranch *pBr = *it;
        YW_ASSERT_INFO(pBr != NULL, "Wrong223");
        double brCurr = pBr->GetLength();
        double brConv = AGGTBranchLengthGrid::Instance().SnapToGrid(brCurr);
        pBr->SetLength(brConv);
    }
    //cout << "SnapBrLens: \n";
    vector<GenealogicalNetworkNode *> listMixNodes;
    network.GetMixNodes( listMixNodes );
//cout << "Num of mixed nodes: " << listMixNodes.size() << endl;
    for(int i=0; i<(int)listMixNodes.size(); ++i)
    {
        double mr = listMixNodes[i]->GetMixRatio();
        double mrConv = AGGTBranchLengthGrid::Instance().SnapMRToGrid(mr);
        listMixNodes[i]->SetMixRatio(mrConv);
    }
}

