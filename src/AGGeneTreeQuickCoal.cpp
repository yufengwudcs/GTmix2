//
//  AGGeneTreeQuickCoal.cpp
//  
//
//  Created by Yufeng Wu on 2/16/24.
//
#include <algorithm>
#include <queue>
#include <stack>
#include <iostream>
#include <cmath>
#include "UtilsNumerical.h"
using namespace std;

#include "AGGeneTreeQuickCoal.hpp"
#include "PhylogenyTreeBasic.h"
#include "GenealogicalNetwork.h"

// *************************************************************************************

static double CalcCombNumAG(int n, int k)
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

// *************************************************************************************
// ways to approximate
static double UtilsCalcCoalProbBranchPoissonApproxAG( int u, int v, double len )
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

static double UtilsCalcCoalProbBranchAG( int u, int v, double len )
{
    // prob of coalsecing from u lineage to v lineages with len time
    // CAUTION: we assume u >= v > =1
//cout << "u:" << u << "v:" << v << endl;
    YW_ASSERT_INFO( u >= v && (u == 0 || v >= 1), "Here, it must be u>=v and v>=1 unless u=0" );
    
    // if u = 0, then nothing
    if( u == 0 )
    {
        return 1.0;
    }
    // YW: test approx
    const int THRES_APPROX_U_VAL = 30;
    const double THRES_APPROX_LEN = 0.0001;
    if( u > THRES_APPROX_U_VAL || len < THRES_APPROX_LEN )
    {
        // use approximation if approximation is needed
        return UtilsCalcCoalProbBranchPoissonApproxAG( u, v, len );
    }
    
//cout << "CalcCoalProbBranch: u=" << u << ", v=" << v << ", len=" << len;
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
        //if( fVerboseMode == true)
        {
            //cout << "*Caution* A minor numerical issue detected: change value " << res << " to 0.0, u=" << u << ", v=" << v << ", len=" << len << "\n";
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

#if 0
static double UtilsCalcCoalProbBranchNormalApproxAG( int u, int v, double len )
{
    // from Griffiths, 1984
    double a = 0.5*u*len;
    double b = -0.5*len;
    double c = a*b/((a+b)*exp(b)-a);
    double valMean = 2.0*c/len;
    double valVar = valMean*(c+b)*(c+b)*(1.0+  c/(c+b) - c/a -c/(a+b) -2.0*c  )/(b*b);
    double valSD = sqrt(valVar);
    
cout << "UtilsCalcCoalProbBranchNormalApprox: a=" << a << ", b=" << b << ", c=" << c << ", valMean=" << valMean << ", valVar=" << valVar << endl;
    
    // now convert to standard nomral
    double v1 = ( v - valMean )/valSD;
    
    // use CDF of normal approximation
    double v2 = (v-1.0-valMean)/valSD;
    
    return CalcApproxCDFStdNormal(v1) - CalcApproxCDFStdNormal(v2);
}
#endif


//***********************************************************************************
// Placing coalescents in the given tree onto network using parsimony

AGGeneTreeQuickCoal2 :: AGGeneTreeQuickCoal2( GenealogicalNetwork &agIn ) : ag(agIn)
{
    Init();
}

double AGGeneTreeQuickCoal2 ::  CalcGeneTreeProbHeu(PhylogenyTreeBasic &phTree)
{
//return CalcGeneTreeProbHeu2(phTree);
return CalcGeneTreeProbHeu4(phTree);
    // for each AG branch, what is the expected number of branches at bottom and top of this branch
    std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsBottom;
    std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsTop;
    
    // prob of a GT lineage not involved in coalescents from one point of AG to an ancestral position of AG
    std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > mapGTLinNotCoalFromDescToAncInAG;
    
    // return log-prob
    // perform init for the tree
    // init expected lin nums
    InitExpNumLins(phTree, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    // init GT lin not coalseced
    InitGTLinNotCoals(mapGTLinNotCoalFromDescToAncInAG, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    
    // now compute
    map< TreeNode *, map< GenealogicalNetworkBranch *, double> > tblProbGTSubtreeAt;
    
    // first, set prob 1.0 all GT leaf lineages to the bottom of all AG leaf branch
    vector<TreeNode *> setAllNodesGT;
    phTree.GetAllNodes(setAllNodesGT);
    map<int, set<TreeNode *> > mapTaxaLeaves;
    for(auto x: setAllNodesGT )
    {
        if( x->IsLeaf() )
        {
            int pop = x->GetIntLabel();
            mapTaxaLeaves[pop].insert(x);
        }
    }
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    for(auto x : setAllBranches)
    {
        if( x->IsLeafBranch() )
        {
            GenealogicalNetworkNode *pn = x->GetDestNode();
            int pop = pn->GetTaxonId();
            YW_ASSERT_INFO( mapTaxaLeaves.find(pop) != mapTaxaLeaves.end(), "Fail to find pop" );
            
            //
            for(auto y: mapTaxaLeaves[pop] )
            {
                tblProbGTSubtreeAt[y][x] = 1.0;
            }
        }
    }
    
    // cache prob of triple AG br
    map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, double> > > mapTripleAGBrs;
    
    
    // add the root branch
    setAllBranches.insert(NULL);
    
    // consider each coalescent in GT
    for(unsigned int i=0; i< setAllNodesGT.size(); ++i )
    {
        TreeNode *pnt = setAllNodesGT[i];
        if( pnt->IsLeaf() )
        {
            continue;
        }
        // need to find out which AG branches can harbor this coalescent
        set<int> setPopsUnder;
        pnt->GetAllDescendIntLbls(setPopsUnder);
        
        // now consider each AG branch
        for(auto x : setAllBranches)
        {
            GenealogicalNetworkBranch *pnagbr = x;
            set<int> setPopsUnderAG;
            if(pnagbr != NULL)
            {
                auto itt = mapPopsAtNodes.find( pnagbr->GetDestNode() );
                YW_ASSERT_INFO( itt != mapPopsAtNodes.end(), "Fail to find 666" );
                setPopsUnderAG = itt->second;
            }
            else
            {
                // use the all pops
                setPopsUnderAG = setAllPops;
            }
            
            //
            if( IsSetContainer( setPopsUnderAG, setPopsUnder ) )
            {
                // now consider its two children in GT
                YW_ASSERT_INFO( pnt->GetChildrenNum()==2, "Not binary tree" );
                TreeNode *pnc1 = pnt->GetChild(0);
                TreeNode *pnc2 = pnt->GetChild(1);
                YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc1) != tblProbGTSubtreeAt.end(), "Child1 not processed yet" );
                YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc2) != tblProbGTSubtreeAt.end(), "Child2 not processed yet" );
                
                // collect all valid child br choices
                vector<pair<GenealogicalNetworkBranch *, pair<double,double> > > listBrProbsChild2;
                for(auto x2: tblProbGTSubtreeAt[pnc2] )
                {
                    GenealogicalNetworkBranch *pnagbrc2 = x2.first;
                    double probc2 = x2.second;
                    auto itc2 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc2].find(pnagbr);
                    if(itc2 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc2].end() )
                    {
                        // invalid
                        continue;
                    }
                    listBrProbsChild2.push_back( std::make_pair( x2.first, std::make_pair(probc2, itc2->second) ) );
                }
                
                double probTot = 0.0;
                
                // now consider where to place the two
                for(auto x1: tblProbGTSubtreeAt[pnc1] )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = x1.first;
                    double probc1 = x1.second;
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    //for(auto x2 : tblProbGTSubtreeAt[pnc2] )
                    for(auto x2 : listBrProbsChild2)
                    {
                        GenealogicalNetworkBranch *pnagbrc2 = x2.first;
                        //double probc2 = x2.second;
                        double probc2 = x2.second.first;
                        //auto itc2 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc2].find(pnagbr);
                        //if( itc2 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc2].end() )
                        //{
                        //    continue;
#if 0
                            cout << "WRONG! pnagbr: ";
                            if( pnagbr != NULL )
                            {
                                pnagbr->Dump();
                            }
                            else cout << "  ROOT branch ";
                            cout << endl;
                            cout << "    pnagbrc2: ";
                            if( pnagbrc2 != NULL )
                            {
                                pnagbrc2->Dump();
                            }
                            else cout << "  ROOT branch ";
                            cout << endl;
#endif
                        //}
                        //YW_ASSERT_INFO(itc2 != mapGTLinNotCoalFromDescToAncInAG[pnagbrc2].end(), "Fail to find no-coal prob2");
                        //double probc2NoCoal = itc2->second;
                        double probc2NoCoal = x2.second.second;
                        
                        double probThisBr = 0.0;
                        bool fFound = false;
                        // see if it is already in cache
                        auto itxx1 = mapTripleAGBrs.find(pnagbr);
                        if( itxx1 != mapTripleAGBrs.end() )
                        {
                            auto itxx2 = itxx1->second.find(pnagbrc1);
                            if( itxx2 != itxx1->second.end() )
                            {
                                auto itxx3 = itxx2->second.find(pnagbrc2);
                                if( itxx3 != itxx2->second.end())
                                {
                                    fFound = true;
                                    probThisBr = itxx3->second;
                                }
                            }
                        }
                        if( fFound == false )
                        {
                            probThisBr = CalcCoalProbAtBranch( pnagbr, pnagbrc1, pnagbrc2, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                            mapTripleAGBrs[pnagbr][pnagbrc1][pnagbrc2] = probThisBr;
                        }
                        // add into
                        probTot += probc1*probc1NoCoal*probc2*probc2NoCoal*probThisBr;
                    }
                }
                
                // save it
                tblProbGTSubtreeAt[pnt][pnagbr] = probTot;
            }
        }
    }
#if 0
    cout << "CalcGeneTreeProbHeu: DUMP all probs\n";
    for(auto x : tblProbGTSubtreeAt )
    {
        cout << "GT node: ";
        set<int> ss;
        x.first->GetDescendentLabelSet(ss);
        DumpIntSet(ss);
        x.first->Dump();
        for(auto y: x.second )
        {
            cout << "prob: " << y.second << ", AG branch: ";
            if( y.first != NULL )
            {
                y.first->Dump();
            }
            else cout << "  ROOT branch ";
            cout << endl;
        }
    }
#endif
    // sum over for all GT root's prob at different AG branch
    double res = 0.0;
    for(auto x : tblProbGTSubtreeAt[phTree.GetRoot()] )
    {
        res += x.second;
    }
    return log(res);
}

double AGGeneTreeQuickCoal2 ::  CalcGeneTreeProbHeu2(PhylogenyTreeBasic &phTree)
{
return CalcGeneTreeProbHeu3(phTree);
    // for each AG branch, what is the expected number of branches at bottom and top of this branch
    std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsBottom;
    std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsTop;
    
    // prob of a GT lineage not involved in coalescents from one point of AG to an ancestral position of AG
    std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > mapGTLinNotCoalFromDescToAncInAG;
    
    // return log-prob
    // perform init for the tree
    // init expected lin nums
    InitExpNumLins(phTree, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    // init GT lin not coalseced
    InitGTLinNotCoals(mapGTLinNotCoalFromDescToAncInAG, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    
    // now compute
    map< TreeNode *, map< GenealogicalNetworkBranch *, double> > tblProbGTSubtreeAt;
    
    // first, set prob 1.0 all GT leaf lineages to the bottom of all AG leaf branch
    vector<TreeNode *> setAllNodesGT;
    phTree.GetAllNodes(setAllNodesGT);
    map<int, set<TreeNode *> > mapTaxaLeaves;
    for(auto x: setAllNodesGT )
    {
        if( x->IsLeaf() )
        {
            int pop = x->GetIntLabel();
            mapTaxaLeaves[pop].insert(x);
        }
    }
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    for(auto x : setAllBranches)
    {
        if( x->IsLeafBranch() )
        {
            GenealogicalNetworkNode *pn = x->GetDestNode();
            int pop = pn->GetTaxonId();
            YW_ASSERT_INFO( mapTaxaLeaves.find(pop) != mapTaxaLeaves.end(), "Fail to find pop" );
            
            //
            for(auto y: mapTaxaLeaves[pop] )
            {
                tblProbGTSubtreeAt[y][x] = 1.0;
            }
        }
    }
    
    // cache prob of triple AG br
    //map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, double> > > mapTripleAGBrs;
    
    // add the root branch
    setAllBranches.insert(NULL);
    
    // pre-calc prob of branches; this is to speedup computation
    map<GenealogicalNetworkBranch *, double> mapAGBrProbPre0;
    //, mapAGBrProbPre1, mapAGBrProbPre2;
    for(auto x : setAllBranches)
    {
        GenealogicalNetworkBranch *pnagbr = x;
        // first, three lineages at the same AGbr
        double p0 = CalcCoalProbAtBranch( pnagbr, pnagbr, pnagbr, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
        mapAGBrProbPre0[pnagbr] = p0;
#if 0
        GenealogicalNetworkBranch *pbrOther = NULL;
        if( pnagbr == NULL )
        {
            pbrOther = *setAllBranches.begin();
            if( pbrOther == NULL )
            {
                pbrOther = *setAllBranches.rbegin();
                YW_ASSERT_INFO( pbrOther != NULL, "Still null" );
            }
        }
        // one side is on this branch, the other isn't
        double p1 = CalcCoalProbAtBranch( pnagbr, pnagbr, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
        mapAGBrProbPre1[pnagbr] = p1;
        
        // now neither at this agbr
        double p2 = CalcCoalProbAtBranch( pnagbr, pbrOther, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
        mapAGBrProbPre2[pnagbr] = p2;
#endif
    }

    // consider each coalescent in GT
    for(unsigned int i=0; i< setAllNodesGT.size(); ++i )
    {
        TreeNode *pnt = setAllNodesGT[i];
        if( pnt->IsLeaf() )
        {
            continue;
        }
        // need to find out which AG branches can harbor this coalescent
        set<int> setPopsUnder;
        pnt->GetAllDescendIntLbls(setPopsUnder);
        
        // now consider each AG branch
        for(auto x : setAllBranches)
        {
            GenealogicalNetworkBranch *pnagbr = x;
            set<int> setPopsUnderAG;
            if(pnagbr != NULL)
            {
                auto itt = mapPopsAtNodes.find( pnagbr->GetDestNode() );
                YW_ASSERT_INFO( itt != mapPopsAtNodes.end(), "Fail to find 666" );
                setPopsUnderAG = itt->second;
            }
            else
            {
                // use the all pops
                setPopsUnderAG = setAllPops;
            }
            
            //
            if( IsSetContainer( setPopsUnderAG, setPopsUnder ) )
            {
                // now consider its two children in GT
                YW_ASSERT_INFO( pnt->GetChildrenNum()==2, "Not binary tree" );
                TreeNode *pnc1 = pnt->GetChild(0);
                TreeNode *pnc2 = pnt->GetChild(1);
                YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc1) != tblProbGTSubtreeAt.end(), "Child1 not processed yet" );
                YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc2) != tblProbGTSubtreeAt.end(), "Child2 not processed yet" );
                
                double probTot = 0.0;
                
                // first do the aspect where both children are on the same branch
                auto &tblc1 = tblProbGTSubtreeAt[pnc1];
                if(tblc1.find(pnagbr)==tblc1.end())
                {
#if 0
                    cout << "pnagbr: ";
                    if( pnagbr != NULL )
                    {
                        pnagbr->Dump();
                    }
                    else
                    {
                        cout << "  ROOT branch ";
                    }
                    cout << endl;
                    cout << "pnt: ";
                    pnt->Dump();
                    cout << "pops: ";
                    DumpIntSet(setPopsUnder);
#endif
                }
                YW_ASSERT_INFO( pnc1->IsLeaf() || tblc1.find(pnagbr)!=tblc1.end(), "Cannot find current agbr1" );
                auto &tblc2 = tblProbGTSubtreeAt[pnc2];
                YW_ASSERT_INFO( pnc2->IsLeaf() || tblc2.find(pnagbr)!=tblc2.end(), "Cannot find current agbr2" );
                double probc1this = 0.0, probc2this = 0.0;
                if( tblc1.find(pnagbr)!=tblc1.end() )
                {
                    probc1this = tblc1[pnagbr];
                }
                if( tblc2.find(pnagbr)!=tblc2.end() )
                {
                    probc2this = tblc2[pnagbr];
                }
#if 0  // not needed
                //probTot += probc1this*probc2this*CalcCoalProbAtBranch( pnagbr, pnagbr, pnagbr, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                probTot += probc1this*probc2this*mapAGBrProbPre0[pnagbr];
#endif
                // then consider one side of child is pnagbr
                //GenealogicalNetworkBranch *pbrOther = NULL;
                //if( pnagbr == NULL )
                //{
                //    pbrOther = *setAllBranches.begin();
                //    if( pbrOther == NULL )
                //    {
                //        pbrOther = *setAllBranches.rbegin();
                //        YW_ASSERT_INFO( pbrOther != NULL, "Still null" );
                //    }
                //}
                //double probC1 = CalcCoalProbAtBranch( pnagbr, pnagbr, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                //double probC1 = mapAGBrProbPre1[ pnagbr];
                double probC1 = mapAGBrProbPre0[ pnagbr];
                double probC1OneSide1 = 0.0;
                for(auto x1: tblProbGTSubtreeAt[pnc1] )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = x1.first;
#if 0           // not neeed
                    if( pnagbrc1 == pnagbr)
                    {
                        continue;
                    }
#endif
                    double probc1 = x1.second;
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    probC1OneSide1 +=  probc1 * probc1NoCoal;
                }
                //probTot += probC1OneSide1 * probc2this * probC1;
                double probC1OneSide2 = 0.0;
                for(auto x1: tblProbGTSubtreeAt[pnc2] )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = x1.first;
#if 0
                    if( pnagbrc1 == pnagbr)
                    {
                        continue;
                    }
#endif
                    double probc1 = x1.second;
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    probC1OneSide2 +=  probc1 * probc1NoCoal;
                }
                //probTot += probC1OneSide2 * probc1this * probC1;
                
                // finally neither branch are in the children
                //probTot += probC1OneSide1 * probC1OneSide2 * CalcCoalProbAtBranch( pnagbr, pbrOther, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                //probTot += probC1OneSide1 * probC1OneSide2 * mapAGBrProbPre2[ pnagbr];
                probTot += probC1OneSide1 * probC1OneSide2 * mapAGBrProbPre0[ pnagbr];
                
                // save it
                tblProbGTSubtreeAt[pnt][pnagbr] = probTot;
            }
        }
    }
#if 0
    cout << "CalcGeneTreeProbHeu: DUMP all probs\n";
    for(auto x : tblProbGTSubtreeAt )
    {
        cout << "GT node: ";
        set<int> ss;
        x.first->GetDescendentLabelSet(ss);
        DumpIntSet(ss);
        x.first->Dump();
        for(auto y: x.second )
        {
            cout << "prob: " << y.second << ", AG branch: ";
            if( y.first != NULL )
            {
                y.first->Dump();
            }
            else cout << "  ROOT branch ";
            cout << endl;
        }
    }
#endif
    // sum over for all GT root's prob at different AG branch
    double res = 0.0;
    for(auto x : tblProbGTSubtreeAt[phTree.GetRoot()] )
    {
        res += x.second;
    }
    return log(res);
}

double AGGeneTreeQuickCoal2 ::  CalcGeneTreeProbHeu3(PhylogenyTreeBasic &phTree)
{
    // for each AG branch, what is the expected number of branches at bottom and top of this branch
    std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsBottom;
    std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsTop;
    
    // prob of a GT lineage not involved in coalescents from one point of AG to an ancestral position of AG
    std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > mapGTLinNotCoalFromDescToAncInAG;
    
    // get nodes of phTree
    vector<TreeNode *> listNodesTree;
    phTree.GetAllNodes(listNodesTree);
    for(unsigned int i=0; i<listNodesTree.size(); ++i)
    {
        TreeNode *pn = listNodesTree[i];
        pn->SetUserData(i);
    }
    // also get the branches of network
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    vector<GenealogicalNetworkBranch *> vecAllBranches;
    int idUse = 0;
    for(auto x: setAllBranches)
    {
        x->SetUserData(idUse++);
        vecAllBranches.push_back(x);
    }
    // add the NULL root branch
    setAllBranches.insert(NULL);
    vecAllBranches.push_back(NULL);
    
    // return log-prob
    // perform init for the tree
    // init expected lin nums
    InitExpNumLins(phTree, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    // init GT lin not coalseced
    InitGTLinNotCoals(mapGTLinNotCoalFromDescToAncInAG, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    
    // now compute
    vector< vector<double> > tblProbGTSubtreeAt( listNodesTree.size() );
    for(unsigned int i=0; i<tblProbGTSubtreeAt.size(); ++i)
    {
        tblProbGTSubtreeAt[i].resize(setAllBranches.size());
        for(unsigned int j=0; j<setAllBranches.size(); ++j)
        {
            tblProbGTSubtreeAt[i][j] = 0.0;
        }
    }
    
    // first, set prob 1.0 all GT leaf lineages to the bottom of all AG leaf branch
    vector<TreeNode *> setAllNodesGT;
    phTree.GetAllNodes(setAllNodesGT);
    map<int, set<TreeNode *> > mapTaxaLeaves;
    for(auto x: setAllNodesGT )
    {
        if( x->IsLeaf() )
        {
            int pop = x->GetIntLabel();
            mapTaxaLeaves[pop].insert(x);
        }
    }

    for(auto x : setAllBranches)
    {
        if( x != NULL && x->IsLeafBranch() )
        {
            GenealogicalNetworkNode *pn = x->GetDestNode();
            int pop = pn->GetTaxonId();
            YW_ASSERT_INFO( mapTaxaLeaves.find(pop) != mapTaxaLeaves.end(), "Fail to find pop" );
            
            //
            for(auto y: mapTaxaLeaves[pop] )
            {
                tblProbGTSubtreeAt[y->GetUserData()][x->GetUserData()] = 1.0;
            }
        }
    }
    
    // cache prob of triple AG br
    //map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, double> > > mapTripleAGBrs;

    
    // pre-calc prob of branches; this is to speedup computation
    vector< double> mapAGBrProbPre0( setAllBranches.size() );
    //, mapAGBrProbPre1, mapAGBrProbPre2;
    //for(auto x : setAllBranches)
    for(unsigned int ii=0; ii<vecAllBranches.size(); ++ii)
    {
        GenealogicalNetworkBranch *pnagbr = vecAllBranches[ii];
        // first, three lineages at the same AGbr
        double p0 = CalcCoalProbAtBranch( pnagbr, pnagbr, pnagbr, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
        mapAGBrProbPre0[ii] = p0;
    }

    // consider each coalescent in GT
    for(unsigned int i=0; i< setAllNodesGT.size(); ++i )
    {
        TreeNode *pnt = setAllNodesGT[i];
        if( pnt->IsLeaf() )
        {
            continue;
        }
        // need to find out which AG branches can harbor this coalescent
        set<int> setPopsUnder;
        pnt->GetAllDescendIntLbls(setPopsUnder);
        
        // now consider each AG branch
        for(auto x : setAllBranches)
        {
            GenealogicalNetworkBranch *pnagbr = x;
            int pnagbrIndex = idUse;
            if( pnagbr != NULL )
            {
                pnagbrIndex = pnagbr->GetUserData();
            }
            set<int> setPopsUnderAG;
            if(pnagbr != NULL)
            {
                auto itt = mapPopsAtNodes.find( pnagbr->GetDestNode() );
                YW_ASSERT_INFO( itt != mapPopsAtNodes.end(), "Fail to find 666" );
                setPopsUnderAG = itt->second;
            }
            else
            {
                // use the all pops
                setPopsUnderAG = setAllPops;
            }
            
            //
            if( pnagbr == NULL || IsSetContainer( setPopsUnderAG, setPopsUnder ) )
            {
                // now consider its two children in GT
                YW_ASSERT_INFO( pnt->GetChildrenNum()==2, "Not binary tree" );
                TreeNode *pnc1 = pnt->GetChild(0);
                TreeNode *pnc2 = pnt->GetChild(1);
                //YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc1) != tblProbGTSubtreeAt.end(), "Child1 not processed yet" );
                //YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc2) != tblProbGTSubtreeAt.end(), "Child2 not processed yet" );
                
                double probTot = 0.0;
                
                // first do the aspect where both children are on the same branch
                auto &tblc1 = tblProbGTSubtreeAt[pnc1->GetUserData()];
                //if(tblc1.find(pnagbr)==tblc1.end())
                //{
#if 0
                    cout << "pnagbr: ";
                    if( pnagbr != NULL )
                    {
                        pnagbr->Dump();
                    }
                    else
                    {
                        cout << "  ROOT branch ";
                    }
                    cout << endl;
                    cout << "pnt: ";
                    pnt->Dump();
                    cout << "pops: ";
                    DumpIntSet(setPopsUnder);
#endif
                //}
                //YW_ASSERT_INFO( pnc1->IsLeaf() || tblc1.find(pnagbr)!=tblc1.end(), "Cannot find current agbr1" );
                auto &tblc2 = tblProbGTSubtreeAt[pnc2->GetUserData()];
                //YW_ASSERT_INFO( pnc2->IsLeaf() || tblc2.find(pnagbr)!=tblc2.end(), "Cannot find current agbr2" );
                //double probc1this = 0.0, probc2this = 0.0;
                //if( tblc1.find(pnagbr)!=tblc1.end() )
                //{
                //    probc1this = tblc1[pnagbr];
                //}
                //if( tblc2.find(pnagbr)!=tblc2.end() )
                //{
                //    probc2this = tblc2[pnagbr];
                //}
                double probc1this = tblc1[pnagbrIndex], probc2this = tblc2[pnagbrIndex];
#if 0  // not needed
                //probTot += probc1this*probc2this*CalcCoalProbAtBranch( pnagbr, pnagbr, pnagbr, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                probTot += probc1this*probc2this*mapAGBrProbPre0[pnagbr];
#endif
                // then consider one side of child is pnagbr
                //GenealogicalNetworkBranch *pbrOther = NULL;
                //if( pnagbr == NULL )
                //{
                //    pbrOther = *setAllBranches.begin();
                //    if( pbrOther == NULL )
                //    {
                //        pbrOther = *setAllBranches.rbegin();
                //        YW_ASSERT_INFO( pbrOther != NULL, "Still null" );
                //    }
                //}
                //double probC1 = CalcCoalProbAtBranch( pnagbr, pnagbr, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                //double probC1 = mapAGBrProbPre1[ pnagbr];
                double probC1 = mapAGBrProbPre0[ pnagbrIndex];
                double probC1OneSide1 = 0.0;
                for(unsigned int jj = 0; jj< tblProbGTSubtreeAt[pnc1->GetUserData()].size(); ++jj )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = vecAllBranches[jj];
#if 0           // not neeed
                    if( pnagbrc1 == pnagbr)
                    {
                        continue;
                    }
#endif
                    //double probc1 = x1.second;
                    double probc1 = tblProbGTSubtreeAt[pnc1->GetUserData()][jj];
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    probC1OneSide1 +=  probc1 * probc1NoCoal;
                }
                //probTot += probC1OneSide1 * probc2this * probC1;
                double probC1OneSide2 = 0.0;
                for(unsigned int jj = 0; jj< tblProbGTSubtreeAt[pnc2->GetUserData()].size(); ++jj )
                //for(auto x1: tblProbGTSubtreeAt[pnc2->GetUserData()] )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = vecAllBranches[jj];
                    //GenealogicalNetworkBranch *pnagbrc1 = x1.first;
#if 0
                    if( pnagbrc1 == pnagbr)
                    {
                        continue;
                    }
#endif
                    double probc1 = tblProbGTSubtreeAt[pnc2->GetUserData()][jj];
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    probC1OneSide2 +=  probc1 * probc1NoCoal;
                }
                //probTot += probC1OneSide2 * probc1this * probC1;
                
                // finally neither branch are in the children
                //probTot += probC1OneSide1 * probC1OneSide2 * CalcCoalProbAtBranch( pnagbr, pbrOther, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                //probTot += probC1OneSide1 * probC1OneSide2 * mapAGBrProbPre2[ pnagbr];
                probTot += probC1OneSide1 * probC1OneSide2 * mapAGBrProbPre0[ pnagbrIndex];
                
                // save it
                tblProbGTSubtreeAt[pnt->GetUserData()][pnagbrIndex] = probTot;
            }
        }
    }
#if 0
    cout << "CalcGeneTreeProbHeu: DUMP all probs\n";
    for(auto x : tblProbGTSubtreeAt )
    {
        cout << "GT node: ";
        set<int> ss;
        x.first->GetDescendentLabelSet(ss);
        DumpIntSet(ss);
        x.first->Dump();
        for(auto y: x.second )
        {
            cout << "prob: " << y.second << ", AG branch: ";
            if( y.first != NULL )
            {
                y.first->Dump();
            }
            else cout << "  ROOT branch ";
            cout << endl;
        }
    }
#endif
    // sum over for all GT root's prob at different AG branch
    double res = 0.0;
    for(auto x : tblProbGTSubtreeAt[phTree.GetRoot()->GetUserData()] )
    {
        //res += x.second;
        res += x;
    }
    return log(res);
}

double AGGeneTreeQuickCoal2 ::  CalcGeneTreeProbHeu4(PhylogenyTreeBasic &phTree)
{
#if 0
cout << "CalcGeneTreeProbHeu4: phTree: ";
string strNW1;
phTree.ConsNewickSorted(strNW1);
cout << strNW1 << endl;
phTree.Dump();
cout << "^^^ Network: ";
ag.Dump();
#endif
    // for each AG branch, what is the prob of having number of branches at bottom and top of this branch
    std::map< GenealogicalNetworkBranch *, std::vector< double> > mapAGBrDistNumLinsBottom;
    std::map< GenealogicalNetworkBranch *, std::vector<double> > mapAGBrDistNumLinsTop;
    
    // prob of a GT lineage not involved in coalescents from one point of AG to an ancestral position of AG
    std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > mapGTLinNotCoalFromDescToAncInAG;
    
    // get nodes of phTree
    vector<TreeNode *> listNodesTree;
    phTree.GetAllNodes(listNodesTree);
    for(unsigned int i=0; i<listNodesTree.size(); ++i)
    {
        TreeNode *pn = listNodesTree[i];
        pn->SetUserData(i);
    }
    // also get the branches of network
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    vector<GenealogicalNetworkBranch *> vecAllBranches;
    int idUse = 0;
    for(auto x: setAllBranches)
    {
        x->SetUserData(idUse++);
        vecAllBranches.push_back(x);
    }
    // add the NULL root branch
    setAllBranches.insert(NULL);
    vecAllBranches.push_back(NULL);
    
    // return log-prob
    // perform init for the tree
    // init expected lin nums
    InitDistNumLins(phTree, mapAGBrDistNumLinsBottom, mapAGBrDistNumLinsTop);
    // init GT lin not coalseced
    InitGTLinNotCoals2(mapGTLinNotCoalFromDescToAncInAG, mapAGBrDistNumLinsBottom, mapAGBrDistNumLinsTop);
    
    // now compute
    vector< vector<double> > tblProbGTSubtreeAt( listNodesTree.size() );
    for(unsigned int i=0; i<tblProbGTSubtreeAt.size(); ++i)
    {
        tblProbGTSubtreeAt[i].resize(setAllBranches.size());
        for(unsigned int j=0; j<setAllBranches.size(); ++j)
        {
            tblProbGTSubtreeAt[i][j] = 0.0;
        }
    }
    
    // first, set prob 1.0 all GT leaf lineages to the bottom of all AG leaf branch
    vector<TreeNode *> setAllNodesGT;
    phTree.GetAllNodes(setAllNodesGT);
    map<int, set<TreeNode *> > mapTaxaLeaves;
    for(auto x: setAllNodesGT )
    {
        if( x->IsLeaf() )
        {
            int pop = x->GetIntLabel();
            mapTaxaLeaves[pop].insert(x);
        }
    }

    for(auto x : setAllBranches)
    {
        if( x != NULL && x->IsLeafBranch() )
        {
            GenealogicalNetworkNode *pn = x->GetDestNode();
            int pop = pn->GetTaxonId();
            YW_ASSERT_INFO( mapTaxaLeaves.find(pop) != mapTaxaLeaves.end(), "Fail to find pop" );
            
            //
            for(auto y: mapTaxaLeaves[pop] )
            {
                tblProbGTSubtreeAt[y->GetUserData()][x->GetUserData()] = 1.0;
            }
        }
    }
    
    // cache prob of triple AG br
    //map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, map<GenealogicalNetworkBranch *, double> > > mapTripleAGBrs;

    
    // pre-calc prob of branches; this is to speedup computation
    vector< double> mapAGBrProbPre0( setAllBranches.size() );
    //, mapAGBrProbPre1, mapAGBrProbPre2;
    //for(auto x : setAllBranches)
    for(unsigned int ii=0; ii<vecAllBranches.size(); ++ii)
    {
        GenealogicalNetworkBranch *pnagbr = vecAllBranches[ii];
        // first, three lineages at the same AGbr
        double p0 = CalcCoalProbAtBranch2( pnagbr, pnagbr, pnagbr, mapAGBrDistNumLinsBottom, mapAGBrDistNumLinsTop );
        mapAGBrProbPre0[ii] = p0;
    }

    // consider each coalescent in GT
    for(unsigned int i=0; i< setAllNodesGT.size(); ++i )
    {
        TreeNode *pnt = setAllNodesGT[i];
        if( pnt->IsLeaf() )
        {
            continue;
        }
        // need to find out which AG branches can harbor this coalescent
        set<int> setPopsUnder;
        pnt->GetAllDescendIntLbls(setPopsUnder);
        
        // now consider each AG branch
        for(auto x : setAllBranches)
        {
            GenealogicalNetworkBranch *pnagbr = x;
            int pnagbrIndex = idUse;
            if( pnagbr != NULL )
            {
                pnagbrIndex = pnagbr->GetUserData();
            }
            set<int> setPopsUnderAG;
            if(pnagbr != NULL)
            {
                auto itt = mapPopsAtNodes.find( pnagbr->GetDestNode() );
                YW_ASSERT_INFO( itt != mapPopsAtNodes.end(), "Fail to find 666" );
                setPopsUnderAG = itt->second;
            }
            else
            {
                // use the all pops
                setPopsUnderAG = setAllPops;
            }
            
            //
            if( pnagbr == NULL || IsSetContainer( setPopsUnderAG, setPopsUnder ) )
            {
                // now consider its two children in GT
                YW_ASSERT_INFO( pnt->GetChildrenNum()==2, "Not binary tree" );
                TreeNode *pnc1 = pnt->GetChild(0);
                TreeNode *pnc2 = pnt->GetChild(1);
                //YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc1) != tblProbGTSubtreeAt.end(), "Child1 not processed yet" );
                //YW_ASSERT_INFO( tblProbGTSubtreeAt.find(pnc2) != tblProbGTSubtreeAt.end(), "Child2 not processed yet" );
                
                double probTot = 0.0;
                
                // first do the aspect where both children are on the same branch
                auto &tblc1 = tblProbGTSubtreeAt[pnc1->GetUserData()];
                //if(tblc1.find(pnagbr)==tblc1.end())
                //{
#if 0
                    cout << "pnagbr: ";
                    if( pnagbr != NULL )
                    {
                        pnagbr->Dump();
                    }
                    else
                    {
                        cout << "  ROOT branch ";
                    }
                    cout << endl;
                    cout << "pnt: ";
                    pnt->Dump();
                    cout << "pops: ";
                    DumpIntSet(setPopsUnder);
#endif
                //}
                //YW_ASSERT_INFO( pnc1->IsLeaf() || tblc1.find(pnagbr)!=tblc1.end(), "Cannot find current agbr1" );
                auto &tblc2 = tblProbGTSubtreeAt[pnc2->GetUserData()];
                //YW_ASSERT_INFO( pnc2->IsLeaf() || tblc2.find(pnagbr)!=tblc2.end(), "Cannot find current agbr2" );
                //double probc1this = 0.0, probc2this = 0.0;
                //if( tblc1.find(pnagbr)!=tblc1.end() )
                //{
                //    probc1this = tblc1[pnagbr];
                //}
                //if( tblc2.find(pnagbr)!=tblc2.end() )
                //{
                //    probc2this = tblc2[pnagbr];
                //}
                double probc1this = tblc1[pnagbrIndex], probc2this = tblc2[pnagbrIndex];
#if 0  // not needed
                //probTot += probc1this*probc2this*CalcCoalProbAtBranch( pnagbr, pnagbr, pnagbr, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                probTot += probc1this*probc2this*mapAGBrProbPre0[pnagbr];
#endif
                // then consider one side of child is pnagbr
                //GenealogicalNetworkBranch *pbrOther = NULL;
                //if( pnagbr == NULL )
                //{
                //    pbrOther = *setAllBranches.begin();
                //    if( pbrOther == NULL )
                //    {
                //        pbrOther = *setAllBranches.rbegin();
                //        YW_ASSERT_INFO( pbrOther != NULL, "Still null" );
                //    }
                //}
                //double probC1 = CalcCoalProbAtBranch( pnagbr, pnagbr, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                //double probC1 = mapAGBrProbPre1[ pnagbr];
                double probC1 = mapAGBrProbPre0[ pnagbrIndex];
                double probC1OneSide1 = 0.0;
                for(unsigned int jj = 0; jj< tblProbGTSubtreeAt[pnc1->GetUserData()].size(); ++jj )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = vecAllBranches[jj];
#if 0           // not neeed
                    if( pnagbrc1 == pnagbr)
                    {
                        continue;
                    }
#endif
                    //double probc1 = x1.second;
                    double probc1 = tblProbGTSubtreeAt[pnc1->GetUserData()][jj];
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    probC1OneSide1 +=  probc1 * probc1NoCoal;
                }
                //probTot += probC1OneSide1 * probc2this * probC1;
                double probC1OneSide2 = 0.0;
                for(unsigned int jj = 0; jj< tblProbGTSubtreeAt[pnc2->GetUserData()].size(); ++jj )
                //for(auto x1: tblProbGTSubtreeAt[pnc2->GetUserData()] )
                {
                    GenealogicalNetworkBranch *pnagbrc1 = vecAllBranches[jj];
                    //GenealogicalNetworkBranch *pnagbrc1 = x1.first;
#if 0
                    if( pnagbrc1 == pnagbr)
                    {
                        continue;
                    }
#endif
                    double probc1 = tblProbGTSubtreeAt[pnc2->GetUserData()][jj];
                    auto itc1 = mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].find(pnagbr);
                    if(itc1 == mapGTLinNotCoalFromDescToAncInAG[pnagbrc1].end() )
                    {
                        // invalid
                        continue;
                    }
                    double probc1NoCoal = itc1->second;
                    probC1OneSide2 +=  probc1 * probc1NoCoal;
                }
                //probTot += probC1OneSide2 * probc1this * probC1;
                
                // finally neither branch are in the children
                //probTot += probC1OneSide1 * probC1OneSide2 * CalcCoalProbAtBranch( pnagbr, pbrOther, pbrOther, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                //probTot += probC1OneSide1 * probC1OneSide2 * mapAGBrProbPre2[ pnagbr];
                probTot += probC1OneSide1 * probC1OneSide2 * mapAGBrProbPre0[ pnagbrIndex];
                
                // save it
                tblProbGTSubtreeAt[pnt->GetUserData()][pnagbrIndex] = probTot;
            }
        }
    }
#if 0
    cout << "CalcGeneTreeProbHeu: DUMP all probs\n";
    for(auto x : tblProbGTSubtreeAt )
    {
        cout << "GT node: ";
        set<int> ss;
        x.first->GetDescendentLabelSet(ss);
        DumpIntSet(ss);
        x.first->Dump();
        for(auto y: x.second )
        {
            cout << "prob: " << y.second << ", AG branch: ";
            if( y.first != NULL )
            {
                y.first->Dump();
            }
            else cout << "  ROOT branch ";
            cout << endl;
        }
    }
#endif
#if 0
    cout << "CalcGeneTreeProbHeu4: DUMP all probs\n";
    for(unsigned int i=0; i< tblProbGTSubtreeAt.size(); ++i )
    {
        DumpDoubleVec(tblProbGTSubtreeAt[i]);
    }
#endif
    // sum over for all GT root's prob at different AG branch
    double res = 0.0;
    for(auto x : tblProbGTSubtreeAt[phTree.GetRoot()->GetUserData()] )
    {
        //res += x.second;
        res += x;
    }
//cout << "res: " << res << endl;
//exit(1);
    return log(res);
}



//---------------------------------------------
void AGGeneTreeQuickCoal2 :: Init()
{
    // init net info
    // get topologically ordered nodes in network
    ag.GetListofNodesTopdown( listNodesAG );
    std::reverse( listNodesAG.begin(), listNodesAG.end() );
//cout << "Number of nodes in AG: " << listNodesAG.size() << endl;
    
    // init desendant sets
    listDescNodesAG.clear();
    for( unsigned int i=0; i<listNodesAG.size(); ++i )
    {
        // add self
        listDescNodesAG[listNodesAG[i] ].insert( listNodesAG[i] );
    }
    // traverse from bottom up
    for( unsigned int i=0; i<listNodesAG.size(); ++i )
    {
        // update its parents
        GenealogicalNetworkNode *pn = listNodesAG[i];
        GenealogicalNetworkBranch *pparBr1 = pn->GetAnces1();
        GenealogicalNetworkNode *ppar1 = NULL;
        if( pparBr1 != NULL )
        {
            ppar1 = pparBr1->GetSrcNode();
        }
        if( ppar1 != NULL )
        {
            UnionSetsGen( listDescNodesAG[ppar1], listDescNodesAG[pn] );
        }
        GenealogicalNetworkBranch *pparBr2 = pn->GetAnces2();
        GenealogicalNetworkNode *ppar2 = NULL;
        if( pparBr2 != NULL )
        {
            ppar2 = pparBr2->GetSrcNode();
        }
        if( ppar2 != NULL )
        {
            UnionSetsGen( listDescNodesAG[ppar2], listDescNodesAG[pn] );
        }
    }
    // init ancestral nodes
    listAncesAGNodes.clear();
    for(auto x : listDescNodesAG)
    {
        GenealogicalNetworkNode *pnc = x.first;
        for(auto x2 : x.second )
        {
            listAncesAGNodes[x2].insert(pnc);
        }
    }
    
    // init pop subsets
    mapPopsAtNodes.clear();
    // add leaf pop first
    for( unsigned int i=0; i<listNodesAG.size(); ++i )
    {
        // add self
        GenealogicalNetworkNode *pnag = listNodesAG[i];
        if( pnag->IsLeaf() )
        {
            //
            int pop = pnag->GetTaxonId();
            mapPopsAtNodes[pnag].insert(pop);
            
            setAllPops.insert(pop);
        }
    }
    // now populate all ancestors
    // traverse from bottom up
    for( unsigned int i=0; i<listNodesAG.size(); ++i )
    {
        // update its parents
        GenealogicalNetworkNode *pn = listNodesAG[i];
        GenealogicalNetworkBranch *pparBr1 = pn->GetAnces1();
        GenealogicalNetworkNode *ppar1 = NULL;
        if( pparBr1 != NULL )
        {
            ppar1 = pparBr1->GetSrcNode();
        }
        if( ppar1 != NULL )
        {
            UnionSetsGen( mapPopsAtNodes[ppar1], mapPopsAtNodes[pn] );
        }
        GenealogicalNetworkBranch *pparBr2 = pn->GetAnces2();
        GenealogicalNetworkNode *ppar2 = NULL;
        if( pparBr2 != NULL )
        {
            ppar2 = pparBr2->GetSrcNode();
        }
        if( ppar2 != NULL )
        {
            UnionSetsGen( mapPopsAtNodes[ppar2], mapPopsAtNodes[pn] );
        }
    }
#if 0
    cout << "LIST of pops at each node: \n";
    for(auto x: mapPopsAtNodes)
    {
        cout << "Node: ";
        x.first->Dump();
        cout << "   subset of pops: ";
        DumpIntSet(x.second);
    }
    cout << "Done: init\n";
#endif
}


// calc expected number of branches
void AGGeneTreeQuickCoal2 :: InitExpNumLins(PhylogenyTreeBasic &phTree, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop)
{
    //
    mapAGBrExpNumLinsBottom.clear();
    mapAGBrExpNumLinsTop.clear();
    
    // find out allele counts for this tree
    vector< TreeNode *> listLeafNodes;
    phTree.GetAllLeafNodes(listLeafNodes);
    map<int,int> listPopNumLins;
    for(unsigned int i=0; i<listLeafNodes.size(); ++i)
    {
        int lbl = listLeafNodes[i]->GetIntLabel();
        ++listPopNumLins[lbl];
    }
#if 0
    cout << "** Number of alleles per populations: ";
    for(auto x : listPopNumLins )
    {
        cout << "Pop " << x.first << ":" << x.second << "  ";
    }
    cout << endl;
#endif
    
    // first set the bottom of each leaf branch of AG to be the number of input lineages
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    for(auto x : setAllBranches)
    {
        if( x->IsLeafBranch() )
        {
            GenealogicalNetworkNode *pn = x->GetDestNode();
            int pop = pn->GetTaxonId();
            YW_ASSERT_INFO( listPopNumLins.find(pop) != listPopNumLins.end(), "Fail to find pop" );
            int numLins = listPopNumLins[pop];
            mapAGBrExpNumLinsBottom[x] = numLins;
        }
    }
    
    // now move upwards
    std::vector<GenealogicalNetworkBranch *> listBrs;
    GetAllAGBrsBottomUp( listBrs );
#if 0
cout << "^^^^^^ net: ";
ag.Dump();
ag.DumpMargTrees(false);
cout << "List of branches from bottom up: \n";
for(unsigned int i=0; i<listBrs.size(); ++i)
{
listBrs[i]->Dump();
cout << endl;
}
#endif
    // more stuff
    for(unsigned int i=0; i<listBrs.size(); ++i)
    {
        UpdateExpLinInfoBrTop(listBrs[i], mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
    }
#if 0
    // dump out all top lin info
    cout << "Expected numer of lineages at bottom: \n";
    for(auto x : mapAGBrExpNumLinsBottom)
    {
        cout << x.second << "   for AG branch: ";
        if( x.first != NULL )
        {
            x.first->Dump();
        }
        else
        {
            cout << "  ROOT branch";
        }
        cout << endl;
    }
    cout << "Expected numer of lineages at top: \n";
    for(auto x : mapAGBrExpNumLinsTop)
    {
        cout << x.second << "   for AG branch: ";
        x.first->Dump();
        cout << endl;
    }
#endif
}

// calc expected number of branches
void AGGeneTreeQuickCoal2 :: InitDistNumLins(PhylogenyTreeBasic &phTree, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop)
{
    //
    mapAGBrDistNumLinsBottom.clear();
    mapAGBrDistNumLinsTop.clear();
    
    // find out allele counts for this tree
    vector< TreeNode *> listLeafNodes;
    phTree.GetAllLeafNodes(listLeafNodes);
    map<int,int> listPopNumLins;
    for(unsigned int i=0; i<listLeafNodes.size(); ++i)
    {
        int lbl = listLeafNodes[i]->GetIntLabel();
        ++listPopNumLins[lbl];
    }
#if 0
    cout << "** Number of alleles per populations: ";
    for(auto x : listPopNumLins )
    {
        cout << "Pop " << x.first << ":" << x.second << "  ";
    }
    cout << endl;
#endif
    
    // first set the bottom of each leaf branch of AG to be the number of input lineages
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    for(auto x : setAllBranches)
    {
        if( x->IsLeafBranch() )
        {
            GenealogicalNetworkNode *pn = x->GetDestNode();
            int pop = pn->GetTaxonId();
            YW_ASSERT_INFO( listPopNumLins.find(pop) != listPopNumLins.end(), "Fail to find pop" );
            int numLins = listPopNumLins[pop];
            mapAGBrDistNumLinsBottom[x].resize( numLins+1 );
            mapAGBrDistNumLinsBottom[x][numLins] = 1.0;
            // the others: set to 0
            for(int i=0; i<numLins; ++i)
            {
                mapAGBrDistNumLinsBottom[x][i] = 0.0;
            }
        }
    }
    
    // now move upwards
    std::vector<GenealogicalNetworkBranch *> listBrs;
    GetAllAGBrsBottomUp( listBrs );
#if 0
cout << "^^^^^^ net: ";
ag.Dump();
ag.DumpMargTrees(false);
cout << "List of branches from bottom up: \n";
for(unsigned int i=0; i<listBrs.size(); ++i)
{
listBrs[i]->Dump();
cout << endl;
}
#endif

    
    // more stuff
    for(unsigned int i=0; i<listBrs.size(); ++i)
    {
        UpdateDistLinInfoBrTop(listBrs[i], mapAGBrDistNumLinsBottom, mapAGBrDistNumLinsTop);
    }
#if 0
    // dump out all top lin info
    cout << "Distribution of numer of lineages at bottom: \n";
    for(auto x : mapAGBrDistNumLinsBottom)
    {
        DumpDoubleVec(x.second);
        cout << "   for AG branch: ";
        if( x.first != NULL )
        {
            x.first->Dump();
        }
        else
        {
            cout << "  ROOT branch";
        }
        cout << endl;
    }
    cout << "Distribution of numer of lineages at top: \n";
    for(auto x : mapAGBrDistNumLinsTop)
    {
        DumpDoubleVec(x.second);
        cout << "   for AG branch: ";
        x.first->Dump();
        cout << endl;
    }
#endif
}


// init GT lin not coalseced
void AGGeneTreeQuickCoal2 :: InitGTLinNotCoals(std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop)
{
    mapGTLinNotCoalFromDescToAncInAG.clear();
    
    //
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    // add root branch
    setAllBranches.insert(NULL);
    for(auto x : setAllBranches)
    {
//cout << "---- lineage x: ";
//if( x != NULL )
//x->Dump();
//else cout << "  ROOT branch  ";
        // find all lineages that is ancestral to x
        for(auto y : setAllBranches )
        {
//cout << "+++ lineage y: ";
//if( y != NULL )
//y->Dump();
//else cout << "  ROOT branch  ";
            if( x == y || y == NULL || ( y != NULL && x != NULL && y->IsAncestralTo(x) ) )
            {
                CalcProbGTLinNotCoalBetween(x, y, mapGTLinNotCoalFromDescToAncInAG, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop);
//cout << "Done one...\n";
            }
        }
    }
    
#if 0
    // dump out prob
    cout << "\n\n************PROB no coalescent: \n";
    for(auto x : mapGTLinNotCoalFromDescToAncInAG)
    {
        cout << "----Descendant lineage: ";
        if( x.first != NULL )
            x.first->Dump();
        else cout << "  ROOT branch  ";
        cout << endl;
        for(auto y : x.second)
        {
            cout << "prob: " << y.second << ", for ancestral lineage: ";
            if( y.first != NULL )
                y.first->Dump();
            else cout << "  ROOT branch ";
            cout << endl;
        }
    }
#endif
//exit(1);
}

// init GT lin not coalseced
void AGGeneTreeQuickCoal2 :: InitGTLinNotCoals2(std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop)
{
    mapGTLinNotCoalFromDescToAncInAG.clear();
    
    //
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    // add root branch
    setAllBranches.insert(NULL);
    for(auto x : setAllBranches)
    {
//cout << "---- lineage x: ";
//if( x != NULL )
//x->Dump();
//else cout << "  ROOT branch  ";
        // find all lineages that is ancestral to x
        for(auto y : setAllBranches )
        {
//cout << "+++ lineage y: ";
//if( y != NULL )
//y->Dump();
//else cout << "  ROOT branch  ";
            if( x == y || y == NULL || ( y != NULL && x != NULL && y->IsAncestralTo(x) ) )
            {
                CalcProbGTLinNotCoalBetween2(x, y, mapGTLinNotCoalFromDescToAncInAG, mapAGBrDistNumLinsBottom, mapAGBrDistNumLinsTop);
//cout << "Done one...\n";
            }
        }
    }
    
#if 0
    // dump out prob
    cout << "\n\n************PROB no coalescent: \n";
    for(auto x : mapGTLinNotCoalFromDescToAncInAG)
    {
        cout << "----Descendant lineage: ";
        if( x.first != NULL )
            x.first->Dump();
        else cout << "  ROOT branch  ";
        cout << endl;
        for(auto y : x.second)
        {
            cout << "prob: " << y.second << ", for ancestral lineage: ";
            if( y.first != NULL )
                y.first->Dump();
            else cout << "  ROOT branch ";
            cout << endl;
        }
    }
#endif
//exit(1);
}



void AGGeneTreeQuickCoal2 :: UpdateExpLinInfoBrTop(GenealogicalNetworkBranch *pBr, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop)
{
//cout << "Size of mapAGBrExpNumLinsBottom: " << mapAGBrExpNumLinsBottom.size() << endl;
    if( mapAGBrExpNumLinsBottom.find(pBr) == mapAGBrExpNumLinsBottom.end() )
    {
        cout << "pBr: ";
        if( pBr != NULL )
        {
            pBr->Dump();
        }
        else
        {
            cout << "  ROOT\n";
        }
        cout << "mapAGBrExpNumLinsBottom: \n";
        for(auto x: mapAGBrExpNumLinsBottom)
        {
            cout << "branch: ";
            if( x.first != NULL )
                x.first->Dump();
            else cout << "  ROOT ";
            cout << endl;
        }
    }
    
    YW_ASSERT_INFO( mapAGBrExpNumLinsBottom.find(pBr) != mapAGBrExpNumLinsBottom.end(), "Fail to find the branch bottom info" );
    
    double expNumBrBot = mapAGBrExpNumLinsBottom[pBr];
    
    int numLinsUse = std::round(expNumBrBot);
    if( numLinsUse < 1 )
    {
        //numLinsUse = 1;
        //mapAGBrExpNumLinsTop[pBr] = 0.0;
        //return;
    }
    
    // consider all possible coal along this branch
    double brLen = HAP_MAX_INT;
    if( pBr != NULL )
    {
        brLen = pBr->GetLength();
    }
    double expTopLinNums = 0.0;
    for(int t=1; t<=numLinsUse; ++t)
    {
        double probCoal = UtilsCalcCoalProbBranchAG(numLinsUse, t, brLen);
        expTopLinNums += t*probCoal;
    }
    mapAGBrExpNumLinsTop[pBr] = expTopLinNums;
    
    // also update the parental branch's bottom
    GenealogicalNetworkNode *pNodeSrc = pBr->GetSrcNode();
    GenealogicalNetworkBranch *pBrPar = pNodeSrc->GetParentBrSingle();
    if( pNodeSrc->IsMixNode() == false )
    {
        // add up to the NULL root branch
        mapAGBrExpNumLinsBottom[pBrPar] += expTopLinNums;
    }
    else
    {
        // distribute to two
        double ratio = pNodeSrc->GetMixRatio();
        GenealogicalNetworkBranch *pAnc1 = pNodeSrc->GetAnces1();
        GenealogicalNetworkBranch *pAnc2 = pNodeSrc->GetAnces2();
        YW_ASSERT_INFO( pAnc1 != NULL && pAnc2 != NULL, "Ancestors: cannot be null for admixture node" );
        mapAGBrExpNumLinsBottom[pAnc1] += ratio*expTopLinNums;
        mapAGBrExpNumLinsBottom[pAnc2] += (1.0-ratio)*expTopLinNums;
#if 0
cout << "---- mix node: ";
pNodeSrc->Dump();
cout << " anc1: ";
pAnc1->Dump();
cout << " anc2: ";
pAnc2->Dump();
cout << endl;
#endif
    }
}

void AGGeneTreeQuickCoal2 :: UpdateDistLinInfoBrTop(GenealogicalNetworkBranch *pBr, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop)
{
//cout << "Size of mapAGBrDistNumLinsBottom: " << mapAGBrDistNumLinsBottom.size() << endl;
    if( mapAGBrDistNumLinsBottom.find(pBr) == mapAGBrDistNumLinsBottom.end() )
    {
        cout << "pBr: ";
        if( pBr != NULL )
        {
            pBr->Dump();
        }
        else
        {
            cout << "  ROOT\n";
        }
        cout << "Network: ";
        ag.Dump();
        cout << "mapAGBrDistNumLinsBottom: \n";
        for(auto x: mapAGBrDistNumLinsBottom)
        {
            cout << "branch: ";
            if( x.first != NULL )
                x.first->Dump();
            else cout << "  ROOT ";
            cout << endl;
        }
    }
    
    YW_ASSERT_INFO( mapAGBrDistNumLinsBottom.find(pBr) != mapAGBrDistNumLinsBottom.end(), "Fail to find the branch bottom info2" );
    
    // top of branch: same size as bottom
    mapAGBrDistNumLinsTop[pBr].resize( mapAGBrDistNumLinsBottom[pBr].size() );
    for(unsigned int i=0; i<mapAGBrDistNumLinsTop[pBr].size(); ++i)
    {
        mapAGBrDistNumLinsTop[pBr][i] = 0.0;
    }
    
    for(unsigned int ns = 1; ns<mapAGBrDistNumLinsBottom[pBr].size(); ++ns)
    {
        // consider all possible coal along this branch
        double brLen = HAP_MAX_INT;
        if( pBr != NULL )
        {
            brLen = pBr->GetLength();
        }
        for(int t=1; t<=ns; ++t)
        {
            double probCoal = UtilsCalcCoalProbBranchAG(ns, t, brLen);
            mapAGBrDistNumLinsTop[pBr][t] += probCoal*mapAGBrDistNumLinsBottom[pBr][ns];
        }
    }
    
    // also update the parental branch's bottom
    GenealogicalNetworkNode *pNodeSrc = pBr->GetSrcNode();
    GenealogicalNetworkBranch *pBrPar = pNodeSrc->GetParentBrSingle();
    if( pNodeSrc->IsMixNode() == false )
    {
        // add up to the NULL root branch
        //mapAGBrExpNumLinsBottom[pBrPar] += expTopLinNums;
        
        // update bottom if both desc branch's top has been updated
        // get the other desc
        GenealogicalNetworkNode *pNodeOther = pNodeSrc->GetOtherDescNode(pBr->GetDestNode());
        YW_ASSERT_INFO( pNodeOther != NULL, "fail to find 666" );
        //GenealogicalNetworkBranch *pBrOther = pNodeOther->GetParentBrSingle();
        GenealogicalNetworkBranch *pBrOther = ag.GetBranch(pNodeSrc, pNodeOther);
        YW_ASSERT_INFO( pBrOther != NULL, "fail to find 667" );
        if( mapAGBrDistNumLinsTop[pBrOther].size() > 0 )
        {
            // ready to combine
            mapAGBrDistNumLinsBottom[pBrPar].resize( mapAGBrDistNumLinsTop[pBr].size() + mapAGBrDistNumLinsTop[pBrOther].size() -1  );
            for( unsigned int ii=0; ii<mapAGBrDistNumLinsBottom[pBrPar].size(); ++ii )
            {
                mapAGBrDistNumLinsBottom[pBrPar][ii] = 0.0;
            }
            for( unsigned int ii=0; ii<mapAGBrDistNumLinsTop[pBr].size(); ++ii )
            {
                for( unsigned int jj=0; jj<mapAGBrDistNumLinsTop[pBrOther].size(); ++jj )
                {
                    mapAGBrDistNumLinsBottom[pBrPar][ii+jj] += mapAGBrDistNumLinsTop[pBr][ii] * mapAGBrDistNumLinsTop[pBrOther][jj];
                }
            }
        }
    }
    else
    {
        // distribute to two
        double ratio = pNodeSrc->GetMixRatio();
        GenealogicalNetworkBranch *pAnc1 = pNodeSrc->GetAnces1();
        GenealogicalNetworkBranch *pAnc2 = pNodeSrc->GetAnces2();
        YW_ASSERT_INFO( pAnc1 != NULL && pAnc2 != NULL, "Ancestors: cannot be null for admixture node" );
        
        mapAGBrDistNumLinsBottom[pAnc1].resize( mapAGBrDistNumLinsTop[pBr].size() );
        mapAGBrDistNumLinsBottom[pAnc2].resize( mapAGBrDistNumLinsTop[pBr].size() );
        for(unsigned int ii=0; ii<mapAGBrDistNumLinsTop[pBr].size(); ++ii)
        {
            mapAGBrDistNumLinsBottom[pAnc1][ii] = 0.0;
            mapAGBrDistNumLinsBottom[pAnc2][ii] = 0.0;
        }
        for(unsigned int ii=0; ii<mapAGBrDistNumLinsTop[pBr].size(); ++ii)
        {
            // consider all possible way of splitting: jj to pAnc1
            for(unsigned int jj=0; jj<=ii; ++jj)
            {
                int numLinOther = ii-jj;
                double probSplit = CalcNumNChooseKDouble(ii, jj) * pow( ratio, jj ) * pow(1.0-ratio, numLinOther );
                mapAGBrDistNumLinsBottom[pAnc1][jj] += probSplit*mapAGBrDistNumLinsTop[pBr][ii];
                mapAGBrDistNumLinsBottom[pAnc2][numLinOther] += probSplit*mapAGBrDistNumLinsTop[pBr][ii];
            }
        }
        
#if 0
cout << "---- mix node: ";
pNodeSrc->Dump();
cout << " anc1: ";
pAnc1->Dump();
cout << " anc2: ";
pAnc2->Dump();
cout << endl;
#endif
    }
}


bool AGGeneTreeQuickCoal2 :: IsAGNodeAncestralTo(GenealogicalNetworkNode *pnDesc, GenealogicalNetworkNode *pnAnc) const
{
    auto it = listDescNodesAG.find(pnAnc);
    return  it->second.find(pnDesc) != it->second.end();
}
bool AGGeneTreeQuickCoal2 ::  IsAGNodeAncestralTo(const std::vector<GenealogicalNetworkNode *> &listDesc, GenealogicalNetworkNode *pnAnc) const
{
    // is pnDesc descendant to one of the
    for(auto x : listDesc)
    {
        if( IsAGNodeAncestralTo(x, pnAnc) )
        {
            return true;
        }
    }
    return false;
}

void AGGeneTreeQuickCoal2 :: GetAllAGBrsBottomUp( std::vector<GenealogicalNetworkBranch *> &listBrs ) const
{
    listBrs.clear();
    // first add leaf branches
    set<GenealogicalNetworkBranch *> setAllBranches;
    ag.GetAllBranches( setAllBranches );
    for(auto x : setAllBranches)
    {
        if( x->IsLeafBranch() )
        {
            listBrs.push_back(x);
        }
    }
    int numLeafBrs = listBrs.size();
    // now add others
    for(auto x : setAllBranches)
    {
        if( x->IsLeafBranch() == false )
        {
            listBrs.push_back(x);
        }
    }
    // now sort them from bottom to top
    for(int i=numLeafBrs; i<(int)listBrs.size(); ++i)
    {
        for(int j=i+1; j<(int)listBrs.size(); ++j)
        {
            // swap if br[j] is ancestral to br[i]
            if( IsAGNodeAncestralTo(listBrs[j]->GetSrcNode(), listBrs[i]->GetSrcNode() )  )
            {
                GenealogicalNetworkBranch *pb = listBrs[i];
                listBrs[i] = listBrs[j];
                listBrs[j] = pb;
            }
        }
    }
#if 0
    //
    cout << "List of AG branches in bottom-up order: \n";
    for(unsigned int i=0; i<listBrs.size(); ++i)
    {
        listBrs[i]->Dump();
        cout << endl;
    }
#endif
}

double AGGeneTreeQuickCoal2 :: CalcProbGTLinNotCoalBetween( GenealogicalNetworkBranch *pBrDesc, GenealogicalNetworkBranch *pBrAnc, std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop )
{
#if 0
cout << "CalcProbGTLinNotCoalBetween: pBrDesc: ";
if( pBrDesc != NULL )
pBrDesc->Dump();
else cout << "  ROOT branch  ";
cout << "   pBrAnc: ";
if( pBrAnc != NULL )
pBrAnc->Dump();
else cout << "  ROOT branch ";
cout << endl;
#endif
    // calc prob of a GT lineage that traverses from pBrDesc to the (bottom of) pBrAnc and remains un-coalesced; we use the expected num of lineages for coalescet prob
    // assume pBrDesc is a descedant of pBrAnc
    if( pBrDesc == pBrAnc )
    {
//cout << "ret 1.0\n";
        mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = 1.0;
        return 1.0;
    }
//cout << "here0\n";
    // first, it cannot coalesce within this branch
    auto xit = mapAGBrExpNumLinsBottom.find(pBrDesc);
    YW_ASSERT_INFO( xit != mapAGBrExpNumLinsBottom.end(), "Fail to find 444");
    int numLinsBot = std::round(xit->second);
    if( numLinsBot <= 1 )
    {
        //numLinsBot = 1;
        // minimum num of lins = 2
        numLinsBot = 2;
    }
//cout << "here01\n";
    if( numLinsBot == 1 )
    {
//cout << "ret 1.0 too\n";
        // a single lineage cannot coalescent
        mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = 1.0;
        return 1.0;
    }
//cout << "here...\n";
    double probThis = 0.0;
    double brLen = HAP_MAX_INT;
    if( pBrDesc != NULL )
    {
        brLen = pBrDesc->GetLength();
    }
    for(int nl = 2; nl <= numLinsBot; ++nl)
    {
        double prob1 = UtilsCalcCoalProbBranchAG(numLinsBot, nl, brLen);
        // there are numLinsBot - nl coalescents, this lineage cannot involve in any of these
        probThis += prob1 * CalcProbNotCoalescedFor( numLinsBot, nl );
    }
    // if this is the root branch, done
    double res = 0.0;
    if( pBrDesc == NULL )
    {
        res = probThis;
    }
    else
    {
        // otherwise move up
        vector<pair<GenealogicalNetworkBranch *, double> > listPairAnc;
        GenealogicalNetworkNode *pNodeSrc = pBrDesc->GetSrcNode();
        GenealogicalNetworkBranch *pBrPar = pNodeSrc->GetParentBrSingle();
        if( pNodeSrc->IsMixNode() == false )
        {
            // add up to the NULL root branch
            YW_ASSERT_INFO( pBrAnc == NULL || pBrPar == NULL || pBrAnc == pBrPar || pBrPar->GetSrcNode() == pBrAnc->GetDestNode() || IsAGNodeAncestralTo(pBrPar->GetSrcNode(), pBrAnc->GetDestNode()), "Wrong111");
            listPairAnc.push_back( std::make_pair( pBrPar, 1.0) );
        }
        else
        {
            // distribute to two
            double ratio = pNodeSrc->GetMixRatio();
            GenealogicalNetworkBranch *pAnc1 = pNodeSrc->GetAnces1();
            GenealogicalNetworkBranch *pAnc2 = pNodeSrc->GetAnces2();
            YW_ASSERT_INFO( pAnc1 != NULL && pAnc2 != NULL, "Ancestors: cannot be null for admixture node" );
            if( pBrAnc == NULL || pAnc1 == pBrAnc || IsAGNodeAncestralTo(pAnc1->GetSrcNode(), pBrAnc->GetDestNode()) )
            {
                listPairAnc.push_back( std::make_pair( pAnc1, ratio) );;
            }
            if( pBrAnc == NULL || pAnc2 == pBrAnc || IsAGNodeAncestralTo(pAnc2->GetSrcNode(), pBrAnc->GetDestNode()) )
            {
                listPairAnc.push_back( std::make_pair( pAnc2, 1.0-ratio) );
            }
        }
        YW_ASSERT_INFO( listPairAnc.size() > 0, "Fail to find555" );
        //
        double probRec = 0.0;
        for(unsigned int i=0; i<listPairAnc.size(); ++i)
        {
            // see if can do a quick look up
            auto yit = mapGTLinNotCoalFromDescToAncInAG[listPairAnc[i].first].find(pBrAnc);
            double probTemp = 0.0;
            if( yit != mapGTLinNotCoalFromDescToAncInAG[listPairAnc[i].first].end() )
            {
                probTemp = yit->second;
            }
            else
            {
                probTemp = CalcProbGTLinNotCoalBetween( listPairAnc[i].first, pBrAnc, mapGTLinNotCoalFromDescToAncInAG, mapAGBrExpNumLinsBottom, mapAGBrExpNumLinsTop );
                mapGTLinNotCoalFromDescToAncInAG[listPairAnc[i].first][pBrAnc] = probTemp;
            }
            probRec += probTemp * listPairAnc[i].second;
        }
        res =  probThis * probRec;
    }
//cout << "res: " << res << endl;
    mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = res;
    return res;
}

double AGGeneTreeQuickCoal2 :: CalcProbGTLinNotCoalBetween2( GenealogicalNetworkBranch *pBrDesc, GenealogicalNetworkBranch *pBrAnc, std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop )
{
#if 0
    cout << "CalcProbGTLinNotCoalBetween: pBrDesc: ";
    if( pBrDesc != NULL )
        pBrDesc->Dump();
    else cout << "  ROOT branch  ";
    cout << "   pBrAnc: ";
    if( pBrAnc != NULL )
        pBrAnc->Dump();
    else cout << "  ROOT branch ";
    cout << endl;
#endif
    // calc prob of a GT lineage that traverses from pBrDesc to the (bottom of) pBrAnc and remains un-coalesced; we use the expected num of lineages for coalescet prob
    // assume pBrDesc is a descedant of pBrAnc
    if( pBrDesc == pBrAnc )
    {
        //cout << "ret 1.0\n";
        mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = 1.0;
        return 1.0;
    }
    //cout << "here0\n";
    // first, it cannot coalesce within this branch
    auto xit = mapAGBrDistNumLinsBottom.find(pBrDesc);
    YW_ASSERT_INFO( xit != mapAGBrDistNumLinsBottom.end(), "Fail to find 444");
    
    if( xit->second.size() <= 2 )
    {
    //cout << "ret 1.0 too\n";
        // a single lineage cannot coalescent
        mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = 1.0;
        return 1.0;
    }
    
    double brLen = HAP_MAX_INT;
    if( pBrDesc != NULL )
    {
        brLen = pBrDesc->GetLength();
    }
    double probThis = 0.0;
    for(int numLinsBot = 1; numLinsBot<(int)xit->second.size(); ++numLinsBot  )
    {
        //int numLinsBot = std::round(xit->second);
        //if( numLinsBot <= 1 )
        //{
        //numLinsBot = 1;
        // minimum num of lins = 2
        //    numLinsBot = 2;
        //}
        //cout << "here01\n";
        //if( numLinsBot == 1 )
        //{
        //cout << "ret 1.0 too\n";
        //    // a single lineage cannot coalescent
        //    mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = 1.0;
        //     return 1.0;
        //}
        //cout << "here...\n";
        for(int nl = 2; nl <= numLinsBot; ++nl)
        {
            double prob1 = UtilsCalcCoalProbBranchAG(numLinsBot, nl, brLen);
            // there are numLinsBot - nl coalescents, this lineage cannot involve in any of these
            probThis += prob1 * CalcProbNotCoalescedFor( numLinsBot, nl ) * xit->second[numLinsBot];
        }
    }
    // if this is the root branch, done
    double res = 0.0;
    if( pBrDesc == NULL )
    {
        res = probThis;
    }
    else
    {
        // otherwise move up
        vector<pair<GenealogicalNetworkBranch *, double> > listPairAnc;
        GenealogicalNetworkNode *pNodeSrc = pBrDesc->GetSrcNode();
        GenealogicalNetworkBranch *pBrPar = pNodeSrc->GetParentBrSingle();
        if( pNodeSrc->IsMixNode() == false )
        {
            // add up to the NULL root branch
            YW_ASSERT_INFO( pBrAnc == NULL || pBrPar == NULL || pBrAnc == pBrPar || pBrPar->GetSrcNode() == pBrAnc->GetDestNode() || IsAGNodeAncestralTo(pBrPar->GetSrcNode(), pBrAnc->GetDestNode()), "Wrong111");
            listPairAnc.push_back( std::make_pair( pBrPar, 1.0) );
        }
        else
        {
            // distribute to two
            double ratio = pNodeSrc->GetMixRatio();
            GenealogicalNetworkBranch *pAnc1 = pNodeSrc->GetAnces1();
            GenealogicalNetworkBranch *pAnc2 = pNodeSrc->GetAnces2();
            YW_ASSERT_INFO( pAnc1 != NULL && pAnc2 != NULL, "Ancestors: cannot be null for admixture node" );
            if( pBrAnc == NULL || pAnc1 == pBrAnc || IsAGNodeAncestralTo(pAnc1->GetSrcNode(), pBrAnc->GetDestNode()) )
            {
                listPairAnc.push_back( std::make_pair( pAnc1, ratio) );;
            }
            if( pBrAnc == NULL || pAnc2 == pBrAnc || IsAGNodeAncestralTo(pAnc2->GetSrcNode(), pBrAnc->GetDestNode()) )
            {
                listPairAnc.push_back( std::make_pair( pAnc2, 1.0-ratio) );
            }
        }
        YW_ASSERT_INFO( listPairAnc.size() > 0, "Fail to find555" );
        //
        double probRec = 0.0;
        for(unsigned int i=0; i<listPairAnc.size(); ++i)
        {
            // see if can do a quick look up
            auto yit = mapGTLinNotCoalFromDescToAncInAG[listPairAnc[i].first].find(pBrAnc);
            double probTemp = 0.0;
            if( yit != mapGTLinNotCoalFromDescToAncInAG[listPairAnc[i].first].end() )
            {
                probTemp = yit->second;
            }
            else
            {
                probTemp = CalcProbGTLinNotCoalBetween2( listPairAnc[i].first, pBrAnc, mapGTLinNotCoalFromDescToAncInAG, mapAGBrDistNumLinsBottom, mapAGBrDistNumLinsTop );
                mapGTLinNotCoalFromDescToAncInAG[listPairAnc[i].first][pBrAnc] = probTemp;
            }
            probRec += probTemp * listPairAnc[i].second;
        }
        res =  probThis * probRec;
    }
//cout << "res: " << res << endl;
    mapGTLinNotCoalFromDescToAncInAG[pBrDesc][pBrAnc] = res;
    return res;
}


double AGGeneTreeQuickCoal2 :: CalcProbNotCoalescedFor(int m, int n) const
{
    // when m lineages coalesce into n lineages, the prob of a single lineage is not involved
    // when there are k lineages coalescing into k-1 lineages, a specific lineage is not involved occurs with prob = C(k-1,2)/C(k,2). Then we multiple for each m->m-1->...->n
    YW_ASSERT_INFO( m >=n && n >= 1, "m must be no smaller than n");
    if( n == 1 )
    {
        return 0.0;
    }
    if( m == n )
    {
        return 1;
    }
    if( m == n+1 )
    {
        return ((double) (m-2) )/m  ;
    }
    return  ((double) n*(n-1) )/(m*(m-1));
}

double AGGeneTreeQuickCoal2 :: CalcCoalProbAtBranch(GenealogicalNetworkBranch *pBr, GenealogicalNetworkBranch *pBrChild1, GenealogicalNetworkBranch *pBrChild2, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop) const
{
    int numPrior = 0;
    if( pBrChild1 == pBr )
    {
        ++numPrior;
    }
    if( pBrChild2 == pBr )
    {
        ++numPrior;
    }
    //
    //double coeff = 1.0/(numPrior + 1);
    auto itx = mapAGBrExpNumLinsBottom.find(pBr);
    YW_ASSERT_INFO(itx != mapAGBrExpNumLinsBottom.end(), "Num lins bottom: not found");
    int numLins = std::round(itx->second);
    if( numLins < numPrior + 2)
    {
        numLins = numPrior + 2;
    }
    
    // YW: only require at least two lineages
    //if( numLins < 2 )
    //{
    //    numLins = 2;
    //}
    //if( numLins < numPrior/2 + 2)
    //{
    //    numLins = 2+numPrior/2;
    //}
    
    double coeff = 1.0/(numLins - 1);
    double brLen = HAP_MAX_INT;
    if( pBr != NULL )
    {
        brLen = pBr->GetLength();
    }
    // here we assume, each round, there are two lineages that we aim at coalescing
    double probTot = 0.0;
    for( int nlt = 1; nlt <= numLins-1; ++nlt )
    {
        probTot += UtilsCalcCoalProbBranchAG(numLins, nlt, brLen) * CalcCoalCoeff( numLins, numLins-nlt  );
    }
    return coeff * probTot;
}

double AGGeneTreeQuickCoal2 :: CalcCoalProbAtBranch2(GenealogicalNetworkBranch *pBr, GenealogicalNetworkBranch *pBrChild1, GenealogicalNetworkBranch *pBrChild2, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop) const
{
#if 0
cout << "CalcCoalProbAtBranch2: pBr: ";
if( pBr != NULL )
pBr->Dump();
else cout << "  ROOT ";
cout << endl;
#endif
    
    int numPrior = 0;
    if( pBrChild1 == pBr )
    {
        ++numPrior;
    }
    if( pBrChild2 == pBr )
    {
        ++numPrior;
    }
    //
    //double coeff = 1.0/(numPrior + 1);
    auto itx = mapAGBrDistNumLinsBottom.find(pBr);
    YW_ASSERT_INFO(itx != mapAGBrDistNumLinsBottom.end(), "Num lins bottom: not found");
    int numLins = (int)itx->second.size()-1;
    YW_ASSERT_INFO( numLins >= 0, "Cannot be negative" );
    
    // YW: only require at least two lineages
    //if( numLins < 2 )
    //{
    //    numLins = 2;
    //}
    //if( numLins < numPrior/2 + 2)
    //{
    //    numLins = 2+numPrior/2;
    //}
    
    if( numLins < numPrior )
    {
        // if there is only one lineage, cannot coalesce
        return 0.0;
    }
    
    double coeff = 1.0/(numLins - 1);
    double brLen = HAP_MAX_INT;
    if( pBr != NULL )
    {
        brLen = pBr->GetLength();
    }
    // here we assume, each round, there are two lineages that we aim at coalescing
    double probTot = 0.0;
    for(int ns=2; ns<=numLins; ++ns)
    {
        for( int nlt = 1; nlt <= ns; ++nlt )
        {
            double f1 = UtilsCalcCoalProbBranchAG(ns, nlt, brLen) ;
            double f2 = CalcCoalCoeff( ns, ns-nlt  ) ;
            probTot += f1* f2* mapAGBrDistNumLinsBottom[pBr][ns];
//cout << "ns: " << ns << ", brLen:" << brLen  << ", nlt: " << nlt << ", f1:" << f1 << ", f2:" << f2 << ", num-prob: " << mapAGBrDistNumLinsBottom[pBr][ns] << endl;
        }
    }
//cout << "--- coeff: " << coeff << ", probTot: " << probTot << endl;
    return coeff * probTot;
}

double AGGeneTreeQuickCoal2 :: CalcCoalCoeff(int m, int k) const
{
    // m total lineages which is to coalesce k rounds: 2 of them are what we want to coalesce
    // return prob of success of these 2 targetted lins finally coalescing after k rounds
    // (YW: 04/14/24: not necessiarly at the last round)
    YW_ASSERT_INFO( m >= k+1 && k >= 0, "WRONG777" );
    if(k == 0 )
    {
        return 0.0;
    }
    if( m == 2 )
    {
        return 1.0;
    }
    if( m == 3 )
    {
        return 1.0/3;
    }
    //return 1.0/(0.5*m*(m-1)) *( 2.0*(m-2) +  0.5*(m-2)*(m-3)* CalcCoalCoeff(m-1, k-1) );
    return 1.0/(0.5*m*(m-1)) *(1 +  0.5*(m-2)*(m-3)* CalcCoalCoeff(m-1, k-1) );
}

