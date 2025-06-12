//
//  GenealogicalNetworkClades.cpp
//  
//
//  Created by Yufeng Wu on 1/3/20.
//
//

#include "GenealogicalNetworkClades.hpp"
#include "Utils3.h"
#include "Utils4.h"
#include "PhylogenyTreeBasic.h"
#include <cmath>

//***********************************************************************************
// Genealogical network population label


GenealogicalNetworkHapPop :: GenealogicalNetworkHapPop(const vector<pair<string,vector<int> > > &listHapsInPopsIn)
{
    listHapsInPops = listHapsInPopsIn;
    
//#if 0
    // do a conversion: reduce by 1 for each hapotypes
//#if 0
    for(int i=0; i<(int)listHapsInPops.size(); ++i)
    {
        for(int j=0; j<(int)listHapsInPops[i].second.size(); ++j)
        {
            --listHapsInPops[i].second[j];
        }
    }
//#endif
cout << "GenealogicalNetworkHapPop: listHapsInPops:";
for(int i=0; i<(int)listHapsInPops.size(); ++i)
{
cout << listHapsInPops[i].first << ": ";
DumpIntVec(listHapsInPops[i].second);
}
    
    //
    map<int,string> mapHapToPopUse;
    for(int i=0; i<(int)listHapsInPops.size(); ++i)
    {
        for(int j=0; j<(int)listHapsInPops[i].second.size(); ++j)
        {
            mapHapToPopUse[ listHapsInPops[i].second[j] ] = listHapsInPops[i].first;
        }
    }
//#endif
    
    //map<int,string> mapHapToPopUse;
    //for(int i=0; i<(int)listHapsInPops.size(); ++i)
    //{
    //    mapHapToPopUse[i] = listHapsInPops[i].first;
    //}
    
    
    int numHaps = mapHapToPopUse.size();
    this->mapHapToPop.resize(numHaps);
    for(int i=0; i<numHaps; ++i)
    {
        YW_ASSERT_INFO( mapHapToPopUse.find(i) != mapHapToPopUse.end(), "Haplotype is missing" );
        this->mapHapToPop[i] = mapHapToPopUse[ i ];
    }
}

string GenealogicalNetworkHapPop :: ConvHapToPop(int hap) const
{
    YW_ASSERT_INFO( hap<(int)mapHapToPop.size(), "Out of bound" );
    return mapHapToPop[hap];
}

void GenealogicalNetworkHapPop :: Dump() const
{
    cout << "Populaton informaton:\n";
    for(int i=0; i<(int)listHapsInPops.size(); ++i)
    {
        cout << "Population " << listHapsInPops[i].first << ": ";
        DumpIntVec(listHapsInPops[i].second);
    }
}

//***********************************************************************************
// Genealogical network pop-label clade with threshold

bool GenealogicalNetworkConvCladeRounded :: operator<(const GenealogicalNetworkConvCladeRounded &rhs) const
{
    return this->setPopLabels < rhs.setPopLabels;
}

void GenealogicalNetworkConvCladeRounded :: AddPopLabel(const string &strPop)
{
    setPopLabels.insert(strPop);
}

void GenealogicalNetworkConvCladeRounded :: AddPopLabel(const std::string &strPop, int numCopies)
{
    for(int i=0; i<numCopies; ++i)
    {
        setPopLabels.insert(strPop);
    }
}

double GenealogicalNetworkConvCladeRounded :: CalcDist(const GenealogicalNetworkConvCladeRounded &rhs) const
{
    // distance is based on Jacard distance
    multiset<string> sint;
    JoinMultiset( this->setPopLabels, rhs.setPopLabels, sint );
    return ((double)sint.size())/( this->setPopLabels.size()+rhs.setPopLabels.size()-sint.size()  );
}

void GenealogicalNetworkConvCladeRounded :: Dump() const
{
    for(multiset<string> :: iterator it = setPopLabels.begin(); it != setPopLabels.end(); ++it)
    {
        cout << *it << ",";
    }
}

//***********************************************************************************
// Genealogical network pop-label clade


bool GenealogicalNetworkConvClade :: operator<(const GenealogicalNetworkConvClade &rhs) const
{
    return this->mapPopHaps < rhs.mapPopHaps;
}
void GenealogicalNetworkConvClade :: AddHapForPop(const string &strPop)
{
    AddHapForPop(strPop, 1);
}
void GenealogicalNetworkConvClade :: AddHapForPop(const string &strPop, int ncopy)
{
    if(mapPopHaps.find(strPop) == mapPopHaps.end() )
    {
        mapPopHaps[strPop] = ncopy;
    }
    else
    {
        mapPopHaps[strPop] += ncopy;
    }
}

void GenealogicalNetworkConvClade :: RoundedTo(int numTargetLins, GenealogicalNetworkConvCladeRounded &cladeRounded) const
{
    int numLins = GetNumLins();
    double fac = numTargetLins/((double) numLins);
    
    // create proportionally
    map<string,int> mapPopHapsNew;
    map<string, double> mapPopHapsNewFrac;
    int numLinsTot = 0;
    for( map<string,int> :: const_iterator it = mapPopHaps.begin(); it != mapPopHaps.end(); ++it )
    {
        int numLinsStep = (int)( it->second * fac );
        mapPopHapsNew[it->first] = numLinsStep;
        numLinsTot += numLinsStep;
        
        if( numLinsStep != it->second * fac )
        {
            double rem = it->second *fac - numLinsStep;
            YW_ASSERT_INFO(rem >= 0.0, "Cannot be negative");
            mapPopHapsNewFrac[it->first] = rem;
        }
    }
    YW_ASSERT_INFO( numLinsTot <= numTargetLins, "Overflow" );
    // assign remaining lins if there is any
    while( numLinsTot < numTargetLins )
    {
        // find the lin that has largest remainder
        string strRem;
        double remMax = 0.0;
        for( map<string,double> :: iterator it = mapPopHapsNewFrac.begin(); it != mapPopHapsNewFrac.end(); ++it )
        {
            if( remMax < it->second )
            {
                remMax = it->second;
                strRem = it->first;
            }
        }
        YW_ASSERT_INFO(remMax > 0.0, "Fail to find");
        mapPopHapsNewFrac.erase(strRem);
        ++mapPopHapsNew[strRem];
        ++numLinsTot;
    }
    for( map<string,int> :: const_iterator it = mapPopHapsNew.begin(); it != mapPopHapsNew.end(); ++it )
    {
        if( it->second == 0 )
        {
            continue;
        }
        // calculate copy
        cladeRounded.AddPopLabel(it->first, it->second);
    }
    
    
    
#if 0
    double bucketSzUse = numLins*roundBucketSz;
    // check each pop
    for( map<string,int> :: const_iterator it = mapPopHaps.begin(); it != mapPopHaps.end(); ++it )
    {
        // make sure we should keep it
        if( it->second >= thresDiscard*numLins )
        {
            // calculate copy
            double numCopies = it->second/bucketSzUse;
            int numCopiesInt = (int)numCopies ;
            if( numCopies - numCopiesInt > 0.00001  )
            {
                ++numCopiesInt;
            }
            cladeRounded.AddPopLabel(it->first, numCopiesInt);
        }
    }
#endif
}

double GenealogicalNetworkConvClade :: CalcDist(const GenealogicalNetworkConvClade &rhs) const
{
    // use jacard distance
    //GenealogicalNetworkConvClade cladeInt;
    //FindIntersect(rhs, cladeInt);
    //int numInt = cladeInt.GetNumLins();
    //int num1 = GetNumLins();
    //int num2 = rhs.GetNumLins();
    //return 1.0-((double)numInt)/(num1+num2-numInt);
    map<string,double> mapFreq1, mapFreq2;
    CalcFreqSpec(mapFreq1);
    rhs.CalcFreqSpec(mapFreq2);
    return CalcDistFreqSpec( mapFreq1, mapFreq2 );
}

double GenealogicalNetworkConvClade :: CalcDist( const std::set< const GenealogicalNetworkConvClade *> &setCladePtrs1, const std::set< const GenealogicalNetworkConvClade *> &setCladePtrs2 )
{
    double res = 0.0;
    for(set<const GenealogicalNetworkConvClade *> :: const_iterator it1 = setCladePtrs1.begin(); it1 != setCladePtrs1.end(); ++it1)
    {
        for(set<const GenealogicalNetworkConvClade *> :: const_iterator it2 = setCladePtrs2.begin(); it2 != setCladePtrs2.end(); ++it2)
        {
            double distStep = (*it1)->CalcDist( *(*it2) );
            res += distStep;
        }
    }
    return res/( setCladePtrs1.size()*setCladePtrs2.size() );
}

void GenealogicalNetworkConvClade :: Dump() const
{
    for(map<string,int> :: const_iterator it = mapPopHaps.begin(); it != mapPopHaps.end(); ++it)
    {
        cout << it->first << ":" << it->second << "  ";
    }
}

int GenealogicalNetworkConvClade :: GetNumLins() const
{
    int res = 0;
    for(map<string,int> :: const_iterator it = mapPopHaps.begin(); it != mapPopHaps.end(); ++it)
    {
        res += it->second;
    }
    return res;
}

void GenealogicalNetworkConvClade :: FindIntersect(const GenealogicalNetworkConvClade &rhs, GenealogicalNetworkConvClade &cladeInt ) const
{
    for(map<string,int> :: const_iterator it1 = this->mapPopHaps.begin(); it1 != this->mapPopHaps.end(); ++it1)
    {
        map<string,int> :: const_iterator it2 =  rhs.mapPopHaps.find( it1->first ) ;
        if( it2 != rhs.mapPopHaps.end() )
        {
            int numInt = std::min( it1->second, it2->second );
            cladeInt.AddHapForPop(it1->first, numInt);
        }
    }
}

void GenealogicalNetworkConvClade :: CalcFreqSpec( std::map<std::string, double> &mapFreq ) const
{
    //
    double numLins = GetNumLins();
    mapFreq.clear();
    for(map<string,int> :: const_iterator it1 = this->mapPopHaps.begin(); it1 != this->mapPopHaps.end(); ++it1)
    {
        mapFreq[it1->first] = ((double)it1->second)/numLins;
    }
}

double GenealogicalNetworkConvClade :: CalcDistFreqSpec( const std::map<std::string, double> &mapFreq1, const std::map<std::string, double> &mapFreq2 ) const
{
    //
    double distTot = 0.0;
    
    // get the common things first
    int numComponents = 0;
    for(map<string, double> :: const_iterator it1 = mapFreq1.begin(); it1 != mapFreq1.end(); ++it1)
    {
        double valOther = 0.0;
        map<string,double> :: const_iterator it2 = mapFreq2.find( it1->first );
        if( it2 != mapFreq2.end() )
        {
            valOther = it2->second;
        }
        double distDiffStep = it1->second - valOther;
        distTot += distDiffStep* distDiffStep;
        ++numComponents;
    }
    // deal with only in the second
    for(map<string, double> :: const_iterator it2 = mapFreq2.begin(); it2 != mapFreq2.end(); ++it2)
    {
        map<string,double> :: const_iterator it1 = mapFreq1.find( it2->first );
        if( it1 != mapFreq1.end() )
        {
            continue;
        }
        double distDiffStep = it2->second * it2->second;
        distTot += distDiffStep* distDiffStep;
        ++numComponents;
    }
    
    return std::sqrt( distTot/numComponents);
}

void GenealogicalNetworkConvClade :: GetPops(std::set<std::string> &setPops) const
{
    setPops.clear();
    for(map<string,int> :: const_iterator it1 = this->mapPopHaps.begin(); it1 != this->mapPopHaps.end(); ++it1)
    {
        setPops.insert(it1->first);
    }
}

void GenealogicalNetworkConvClade :: FilterMinorAlleles(double thresMAF)
{
    //
    map<string,double> mapFreq;
    CalcFreqSpec(mapFreq);
    double maxFreq = 0.0;
    for(map<string,double> :: iterator it = mapFreq.begin(); it != mapFreq.end(); ++it)
    {
        if( it->second > maxFreq)
        {
            maxFreq = it->second;
        }
    }
    
    for(map<string,double> :: iterator it = mapFreq.begin(); it != mapFreq.end(); ++it)
    {
        if( it->second < thresMAF && it->second < maxFreq )
        {
            mapPopHaps.erase( it->first );
        }
    }
}

//***********************************************************************************
// Genealogical network raw clade

GenealogicalNetworkRawClade :: GenealogicalNetworkRawClade(const std::set<int> &setHapsIn) : setHaps(setHapsIn)
{
}

GenealogicalNetworkRawClade :: GenealogicalNetworkRawClade(const GenealogicalNetworkRawClade &rhs) : setHaps(rhs.setHaps)
{
}

bool GenealogicalNetworkRawClade :: operator<(const GenealogicalNetworkRawClade &rhs) const
{
    return this->setHaps < rhs.setHaps;
}

void GenealogicalNetworkRawClade :: ConvHapCladeTo( const GenealogicalNetworkHapPop &popInfo, GenealogicalNetworkConvClade &cladeConv ) const
{
    for(set<int> :: iterator it = setHaps.begin(); it != setHaps.end(); ++it)
    {
        string strPop = popInfo.ConvHapToPop(*it);
        cladeConv.AddHapForPop(strPop);
    }
cout << "ConvHapCladeTo: raw clade: ";
this->Dump();
cout << "  convereted clade: ";
cladeConv.Dump();
}

void GenealogicalNetworkRawClade :: Dump() const
{
    DumpIntSet(setHaps);
}

//***********************************************************************************
// Genealogical network clustering

GenealogicalNetworkCladeCluster :: GenealogicalNetworkCladeCluster(const GenealogicalNetworkHapPop &hapPopConvIn) : hapPopConv(hapPopConvIn)
{
}

void GenealogicalNetworkCladeCluster :: Clustering( const std::map< GenealogicalNetworkRawClade, int> &listGTClades  )
{
    // create converted clade first
    map<GenealogicalNetworkConvClade, int> listGTCladesConv;
    ConvRawClades( listGTClades, listGTCladesConv );
    
    // clean up
    const double THRES_MAF = 0.1;
    map<GenealogicalNetworkConvClade, int> listGTCladesConvClean;
    for( map<GenealogicalNetworkConvClade, int>::iterator it = listGTCladesConv.begin(); it != listGTCladesConv.end(); ++it )
    {
        GenealogicalNetworkConvClade clade(it->first);
        clade.FilterMinorAlleles(THRES_MAF);
cout << "Before clean: clade=";
it->first.Dump();
cout << "  after clean: clade=";
clade.Dump();
cout << endl;
        listGTCladesConvClean[clade] = it->second;
    }
    listGTCladesConv = listGTCladesConvClean;
    
    
    // Idea: use rounding to narrow down the list of candidates to merge
    map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > setCladePtrsToClus;
    for( map<GenealogicalNetworkConvClade, int> :: const_iterator it = listGTCladesConv.begin(); it != listGTCladesConv.end(); ++it )
    {
        setCladePtrsToClus[&(it->first)].insert( &(it->first) );
    }
    
    const int MAX_PAIRWISE_SZ = 10;
    const double DEF_NOISE_RED = 0.05;
    const double MAX_DIST_ALLOWED = 0.2;
    double roundBucketSzCurr = 0.25;
    //int szPrevClusCandidates = HAP_MAX_INT;
    while( setCladePtrsToClus.size() > MAX_PAIRWISE_SZ )
    {
cout << "** In Clustering:, number of clusters so far: " << setCladePtrsToClus.size() << endl;
        
        std::set< std::set<const GenealogicalNetworkConvClade *> > listClusCladPtrs;
        CoarseClus( roundBucketSzCurr, DEF_NOISE_RED, setCladePtrsToClus, listClusCladPtrs );
cout << "After coarse clustering: clusters: num of corase clades: " << listClusCladPtrs.size()  << endl;
for(set<set<const GenealogicalNetworkConvClade *> > :: iterator it = listClusCladPtrs.begin(); it != listClusCladPtrs.end(); ++it)
{
cout << "ONE CLUSTER: \n";
for(set<const GenealogicalNetworkConvClade *> :: iterator it2 = it->begin(); it2 != it->end(); ++it2)
{
(*it2)->Dump();
}
}
cout << endl;
        
        map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > setCladePtrsToClusNew;
        
        // if the clade is small enough, now perform direct clusing; otherwise, put it back and continue next round
        for( set< std::set<const GenealogicalNetworkConvClade *> > :: iterator it = listClusCladPtrs.begin(); it != listClusCladPtrs.end(); ++it )
        {
            if( it->size() <= MAX_PAIRWISE_SZ )
            {
                map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > setCladePtrsToClusStep;
                for( set<const GenealogicalNetworkConvClade *> :: iterator it3 = it->begin(); it3 != it->end(); ++it3 )
                {
                    const GenealogicalNetworkConvClade *ptr = *it3;
                    YW_ASSERT_INFO( setCladePtrsToClus.find(ptr) != setCladePtrsToClus.end(), "Fail to find1" );
                    setCladePtrsToClusStep[ptr] = setCladePtrsToClus[ptr];
                }
                
cout << "Before fineClus: number of clusters: " << setCladePtrsToClusStep.size() << endl;
                FineClus( MAX_DIST_ALLOWED, setCladePtrsToClusStep );
cout << "AFTER fineClus: number of clusters: " << setCladePtrsToClusStep.size() << endl;
                
                // put the samller stuff to continue clustering
                for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it4 = setCladePtrsToClusStep.begin(); it4 != setCladePtrsToClusStep.end(); ++it4 )
                {
cout << "adding: ";
it4->first->Dump();
                    for(set<const GenealogicalNetworkConvClade *> :: iterator it5 = it4->second.begin(); it5 != it4->second.end(); ++it5)
                    {
                        UnionSetsGen(setCladePtrsToClusNew[it4->first],  setCladePtrsToClus[ *it5 ]  );
                    }
                    
cout << "Size of updated clade list: " << setCladePtrsToClusNew[it4->first].size() << endl;
                }
            }
            else
            {
cout << "Just put it back.\n";
                // put it back
                for( set<const GenealogicalNetworkConvClade *> :: iterator it2 = it->begin(); it2 != it->end(); ++it2)
                {
                    //setCladePtrsToClusNew.insert(*it2);
                    const GenealogicalNetworkConvClade *ptr = *it2;
                    YW_ASSERT_INFO( setCladePtrsToClus.find(ptr) != setCladePtrsToClus.end(), "Fail to find1" );
                    setCladePtrsToClusNew[ptr] = setCladePtrsToClus[ptr];
                }
            }
        }
        
cout << "New clade ptr to clus size: " << setCladePtrsToClusNew.size() << endl;
        
        // if not shrinking, cluster directly and stop
        if( (int)setCladePtrsToClusNew.size() == setCladePtrsToClus.size()  )
        {
            FineClus(MAX_DIST_ALLOWED, setCladePtrsToClusNew);
            // put the samller stuff to continue clustering
            map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > setCladePtrsToClusRes;
            for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it4 = setCladePtrsToClusNew.begin(); it4 != setCladePtrsToClusNew.end(); ++it4 )
            {
                cout << "adding: ";
                it4->first->Dump();
                for(set<const GenealogicalNetworkConvClade *> :: iterator it5 = it4->second.begin(); it5 != it4->second.end(); ++it5)
                {
                    UnionSetsGen(setCladePtrsToClusRes[it4->first],  setCladePtrsToClus[ *it5 ]  );
                }
                
            }
            setCladePtrsToClus = setCladePtrsToClusRes;
cout << "Final size of updated clade list: " << setCladePtrsToClus.size() << endl;
            
            break;
        }
        
        roundBucketSzCurr *= 0.5;
        setCladePtrsToClus = setCladePtrsToClusNew;
        //szPrevClusCandidates = (int)setCladePtrsToClus.size();
    }
    
    // store all remaining clades
    for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it )
    {
        set<string> setPops;
        //it->first->GetPops(setPops);
        FindConsensus(it->second, setPops);
        mapInfClades[setPops] = it->second.size();
    }
    
cout << "FINAL LIST OF CLUSTERED CLADES: \n";
for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it )
{
cout << "ONE CLADE: ";
for(set<const GenealogicalNetworkConvClade *> :: iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
{
(*it2)->Dump();
cout << ", ";
}
cout << endl;
}
}

void GenealogicalNetworkCladeCluster :: FindConsensus( const std::set<const GenealogicalNetworkConvClade *> &setClades, std::set<std::string> &setPops ) const
{
    setPops.clear();
    map<string, int> mapPopFreq;
    for( set<const GenealogicalNetworkConvClade *> :: const_iterator it = setClades.begin(); it != setClades.end(); ++it )
    {
        set<string> setPopsStep;
        (*it)->GetPops(setPopsStep);
        for(set<string> :: iterator it2 = setPopsStep.begin(); it2 != setPopsStep.end(); ++it2)
        {
            ++mapPopFreq[*it2];
        }
    }
    const double THRES_POP_INC = 0.3;
    for(map<string,int> :: iterator it = mapPopFreq.begin(); it != mapPopFreq.end(); ++it)
    {
        if( it->second >= THRES_POP_INC * setClades.size() )
        {
            setPops.insert(it->first);
        }
    }
}

void GenealogicalNetworkCladeCluster :: Dump() const
{
    cout << "Inferred population clades: \n";
    for(map<set<string>, int> :: const_iterator it = mapInfClades.begin(); it != mapInfClades.end(); ++it)
    {
        cout << "[" << it->second << "]: ";
        for(set<string> :: const_iterator it2 = it->first.begin(); it2 != it->first.end(); ++it2)
        {
            cout << *it2 << " ";
        }
        cout << ",   ";
    }
    cout << endl;
}

void GenealogicalNetworkCladeCluster :: ConvRawClades( const std::map< GenealogicalNetworkRawClade, int> &listGTClades, std::map< GenealogicalNetworkConvClade, int> &listGTCladesConv  ) const
{
    listGTCladesConv.clear();
    for( map<GenealogicalNetworkRawClade,int> :: const_iterator it = listGTClades.begin(); it != listGTClades.end(); ++it )
    {
        GenealogicalNetworkConvClade cladeConv;
cout << "ConvRawClades: frequencey:" << it->second << endl;
        it->first.ConvHapCladeTo(this->hapPopConv, cladeConv);
        if( listGTCladesConv.find( cladeConv ) == listGTCladesConv.end() )
        {
            listGTCladesConv[cladeConv] = 0;
        }
        listGTCladesConv[cladeConv] += it->second;
    }
    
cout << "ConvRawClades: raw clades: \n";
for(map<GenealogicalNetworkConvClade, int> :: iterator it = listGTCladesConv.begin(); it != listGTCladesConv.end(); ++it)
{
cout << "[" << it->second << "]: ";
it->first.Dump();
cout << endl;
}
}

void GenealogicalNetworkCladeCluster :: CoarseClus(double roundBucketSz, double thresDiscard, std::map<const GenealogicalNetworkConvClade *, std::set<const GenealogicalNetworkConvClade *> > setCladePtrsToClus, std::set< std::set<const GenealogicalNetworkConvClade *> > &listClusCladPtrs ) const
{
    int numLinsClus = (int)(1.0/roundBucketSz);
    
    map<GenealogicalNetworkConvCladeRounded, set<const GenealogicalNetworkConvClade *> > mapRoundedToOrig;
    for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it )
    {
        GenealogicalNetworkConvCladeRounded roundedClade;
        it->first->RoundedTo( numLinsClus, roundedClade );
cout << "Rounding clade: ";
it->first->Dump();
cout << " to ";
roundedClade.Dump();
        
        mapRoundedToOrig[roundedClade].insert( it->first );
    }
    listClusCladPtrs.clear();
    for( map<GenealogicalNetworkConvCladeRounded, set<const GenealogicalNetworkConvClade *> > :: iterator it = mapRoundedToOrig.begin(); it != mapRoundedToOrig.end(); ++it )
    {
cout << "Rounded clade frequency: " << it->second.size() << ": ";
it->first.Dump();
cout << endl;
        listClusCladPtrs.insert( it->second );
        // remove from orig list
        for( set<const GenealogicalNetworkConvClade *> :: iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2 )
        {
            setCladePtrsToClus.erase( *it2 );
        }
    }
cout << "Done with CoarseClus: cluster sizes: ";
for(std::set< std::set<const GenealogicalNetworkConvClade *> > :: iterator it = listClusCladPtrs.begin(); it != listClusCladPtrs.end(); ++it)
{
cout << it->size() << "  ";
}
cout << endl;
}

void GenealogicalNetworkCladeCluster :: FineClus( double maxDistAllowed, std::map<const GenealogicalNetworkConvClade *, std::set<const GenealogicalNetworkConvClade *> > &setCladePtrsToClus )
{
cout << "FineClus: maxDist:" << maxDistAllowed << ", list of clades: ";
for(map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: const_iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it)
{
(*it).first->Dump();
cout << "  ";
}
cout << endl;
    
    // Approach: use nearest neighbor to cluster
    map<pair<const GenealogicalNetworkConvClade *, const GenealogicalNetworkConvClade *>, double> mapPairDist;
    for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: const_iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it )
    {
        map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: const_iterator it2 = it;
        ++it2;
        for(; it2 != setCladePtrsToClus.end(); ++it2 )
        {
            double dist = GenealogicalNetworkConvClade :: CalcDist( it->second, it2->second );
            //double dist = (*it)->CalcDist(*(*it2));
            pair<const GenealogicalNetworkConvClade *, const GenealogicalNetworkConvClade *> pp(it->first, it2->first), pp2(it2->first, it->first);
            mapPairDist[pp] = dist;
            mapPairDist[pp2] = dist;
cout << "Dist: " << dist << ", for clades: ";
it->first->Dump();
cout << " ";
it2->first->Dump();
cout << endl;
        }
    }
cout << "Done with distance calclation.\n";
    
    // do clustering
    map< const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > mapCladeClus;
    // initially each in their own clus
    for( map<const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: const_iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it )
    {
        mapCladeClus[ (it->first) ].insert(it->first);
    }
    
    // now loop to cluster
    while(mapCladeClus.size() > 1)
    {
        //
        double distMin = HAP_MAX_INT*1.0;
        const GenealogicalNetworkConvClade *pMerge1 = NULL;
        const GenealogicalNetworkConvClade *pMerge2 = NULL;
        for( map< const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it1 = mapCladeClus.begin(); it1 != mapCladeClus.end(); ++it1 )
        {
            map< const GenealogicalNetworkConvClade *, set<const GenealogicalNetworkConvClade *> > :: iterator it2 = it1;
            ++it2;
            for(; it2 != mapCladeClus.end(); ++it2)
            {
                double distStep = CalcAveDistBtw( mapPairDist, it1->second, it2->second );
                if( distMin > distStep)
                {
                    distMin = distStep;
                    pMerge1 = (it1->first);
                    pMerge2 = (it2->first);
                }
            }
        }
        YW_ASSERT_INFO( pMerge1 != NULL && pMerge2 != NULL, "Fail to find" );
cout << "distMin: " << distMin << ", Merging candidate: ";
pMerge1->Dump();
cout << "  ";
pMerge2->Dump();
        
        if( distMin > maxDistAllowed)
        {
cout << "Don't merge: too different\n";
            break;
        }
        
        // merge
        for( set<const GenealogicalNetworkConvClade *> :: const_iterator it2 = mapCladeClus[pMerge2].begin(); it2 != mapCladeClus[pMerge2].end(); ++it2 )
        {
            mapCladeClus[pMerge1].insert(*it2);
        }
        mapCladeClus.erase(pMerge2);
    }
    // save the updated list
    setCladePtrsToClus = mapCladeClus;
    
cout << "Done with FineClus: number of clusters: " << setCladePtrsToClus.size() << endl;
for(  std::map<const GenealogicalNetworkConvClade *, std::set<const GenealogicalNetworkConvClade *> > :: iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it )
{
cout << "--repreenting clade: ";
it->first->Dump();
cout << "   clustered clades: ";
for(set<const GenealogicalNetworkConvClade *> :: iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
{
(*it2)->Dump();
cout << "  ";
}
cout << endl;
}
}

void GenealogicalNetworkCladeCluster :: FormConsensus( const std::set<const GenealogicalNetworkConvClade *> &setCladePtrsToClus, std::set<string> &cladePop ) const
{
    // approach: use those appear in at least 50% of the clades
    cladePop.clear();
    map<string,int> mapFreqPop;
    for(set<const GenealogicalNetworkConvClade *> :: const_iterator it = setCladePtrsToClus.begin(); it != setCladePtrsToClus.end(); ++it)
    {
        map<string,int> mapPops;
        (*it)->GetPops(mapPops);
        for(map<string,int> :: iterator it2 = mapPops.begin(); it2 != mapPops.end(); ++it2)
        {
            mapFreqPop[it2->first] += it2->second;
        }
    }
    for(map<string,int> :: iterator it = mapFreqPop.begin(); it != mapFreqPop.end(); ++it)
    {
        if( it->second >= setCladePtrsToClus.size()/2 )
        {
            cladePop.insert( it->first );
        }
    }
}

double GenealogicalNetworkCladeCluster :: CalcAveDistBtw( std::map<std::pair<const GenealogicalNetworkConvClade *, const GenealogicalNetworkConvClade *>, double> &mapPairDist, const std::set<const GenealogicalNetworkConvClade *> &setCladePtrsToClus1, const std::set<const GenealogicalNetworkConvClade *> &setCladePtrsToClus2 ) const
{
    double distTot = 0.0;
    for( set<const GenealogicalNetworkConvClade *> :: const_iterator it1 = setCladePtrsToClus1.begin(); it1 != setCladePtrsToClus1.end(); ++it1 )
    {
        for( set<const GenealogicalNetworkConvClade *> :: const_iterator it2 = setCladePtrsToClus2.begin(); it2 != setCladePtrsToClus2.end(); ++it2 )
        {
            pair<const GenealogicalNetworkConvClade *, const GenealogicalNetworkConvClade *> pp(*it1, *it2);
            distTot += mapPairDist[pp];
        }
    }
    
    double dist = distTot/(setCladePtrsToClus1.size()*setCladePtrsToClus2.size());
    
cout << "CalcAveDistBtw: " << dist << ", clades in clsuter 1: ";
for(set<const GenealogicalNetworkConvClade *> :: iterator it = setCladePtrsToClus1.begin(); it != setCladePtrsToClus1.end(); ++it)
{
(*it)->Dump();
cout << "  ";
}
cout << "\n   clades in cluster 2: ";
for(set<const GenealogicalNetworkConvClade *> :: iterator it = setCladePtrsToClus2.begin(); it != setCladePtrsToClus2.end(); ++it)
{
(*it)->Dump();
cout << "  ";
}
cout << endl;
    
    return dist;
}

//***********************************************************************************
void TestGNClades(const std::vector<std::pair<std::string,std::vector<int> > > &listPopRows, std::vector<PhylogenyTreeBasic *> &listGeneTreePtrs)
{
    GenealogicalNetworkHapPop popInfo(listPopRows);
    GenealogicalNetworkCladeCluster clusInf(popInfo);
    
    // collect all clades
    map<set<int>, int> mapCladesFreq;
    for(int i=0; i<(int)listGeneTreePtrs.size(); ++i)
    {
        int numLvs = listGeneTreePtrs[i]->GetNumLeaves();
        set<set<int> > setClades;
        listGeneTreePtrs[i]->GetAllClades(setClades);
        for(set<set<int> > :: iterator it2 = setClades.begin(); it2 != setClades.end(); ++it2)
        {

            if( it2->size() > 1 && it2->size() < numLvs-1 )
            {
                ++mapCladesFreq[ *it2 ];
            }
        }
    }
    //
    map< GenealogicalNetworkRawClade, int> listGTClades;
    for(map<set<int>, int> :: iterator it = mapCladesFreq.begin(); it != mapCladesFreq.end(); ++it )
    {
        GenealogicalNetworkRawClade clade( it->first );
        listGTClades[ clade ] = it->second;
    }
    clusInf.Clustering( listGTClades );
    clusInf.Dump();
exit(1);
}
