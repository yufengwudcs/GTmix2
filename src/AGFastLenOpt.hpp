//
//  AGFastLenOpt.hpp
//  
//
//  Created by Yufeng Wu on 3/18/24.
//

#ifndef AGFastLenOpt_hpp
#define AGFastLenOpt_hpp

#include "UtilsNumerical.h"
#include "GenealogicalNetwork.h"
//#include "AGGeneTreeQuickCoal.hpp"
#include <vector>


class PhylogenyTreeBasic;
//class AGGeneTreeProb;
class AGGeneTreeProb2;

//***********************************************************************************
// Length grid

class AGGTBranchLengthGrid
{
public:
    static AGGTBranchLengthGrid &Instance();
    double SnapToGrid(double brOrig, int *posCloseOut = NULL) const;
    double ReducePercison(double brOrig) const;
    double SnapMRToGrid(double mr) const;

private:
    AGGTBranchLengthGrid();
    void Init();
    
    vector<double> listBrLenGridVals;
    vector<double> listMRGridVals;
};

//***********************************************************************************
// Optimize branch lengths of a given network topology wrt gene trees

class AGFastLenOpt : public NumericalAlgoUtils
{
public:
    AGFastLenOpt(GenealogicalNetwork &netIn, vector<PhylogenyTreeBasic *> &listGeneTreesIn, int threadInd = 0);
    ~AGFastLenOpt();
    double CalcProb();
    double Optimize();
    virtual double EvaluateAt(double pt, void *pParam);
    //void SetMultithread(int numT);
    // this is the pre-computed log prob
    void SetLogProbPre(double b) { logprobPre = b; }
    void SetCritcalBrs(const set<GenealogicalNetworkBranch *> & setBranchesCriticalIn) { setBranchesCritical = setBranchesCriticalIn; }
    void SetMaxNumOptRounds(int nopt) { maxNumOptRound = nopt; }
    void SetThreadIndex(int t) { threadIndex = t; }
    int GetThreadIndex() const { return threadIndex; }
    
private:
    double Optimize2();
    double OptimizeBranch(GenealogicalNetworkBranch *brCurr, double logprobPre, double &brLenBest);
    //double CalcProbMT();
    void UpdateProbComputeForBranch( GenealogicalNetworkBranch *pBrToChange, double brNew );
    void UpdateProbComputeForMR( GenealogicalNetworkNode *pNodeMRChange, double mrNew );
    double GetMinBrLen(int round, double brLenCur) const;
    double GetMaxBrLen(int round, double brLenCur) const;
    double GetTolBrLen(int round) const;
    double GetMinMR(int round, double mrCur) const;
    double GetMaxMR(int round, double mrCur) const;
    double GetTolMR(int round) const;
    void SnapBrLens(GenealogicalNetwork &networkCurrent);
    
    GenealogicalNetwork &network;
    vector<PhylogenyTreeBasic *> &listGeneTrees;
    //AGGeneTreeQuickCoal2 * pProbCalc;
    //AGGeneTreeProb *pProbCalc2;
    AGGeneTreeProb2 *pProbCalc2;
    GenealogicalNetworkBranch *pBranchInNetCur;
    GenealogicalNetworkNode *pMixNodeCur;
    set<GenealogicalNetworkBranch *> setBranchesCritical;
    set<GenealogicalNetworkBranch *> setBranchesSkip;
    set<GenealogicalNetworkNode *> setMixNodesSkip;
    map<double, double> setProcessedLens;
    map<double, double> setProcessedMRs;
    //int numThreads;
    double logprobPre;
    int maxNumOptRound;
    int threadIndex;
};


#endif /* AGFastLenOpt_hpp */
