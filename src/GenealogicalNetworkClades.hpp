//
//  GenealogicalNetworkClades.hpp
//  
//
//  Created by Yufeng Wu on 1/3/20.
//
//

#ifndef GenealogicalNetworkClades_hpp
#define GenealogicalNetworkClades_hpp

#include <set>
#include <vector>
#include <map>
#include <string>

//***********************************************************************************
class PhylogenyTreeBasic;

void TestGNClades(const std::vector<std::pair<std::string,std::vector<int> > > &listPopRows, std::vector<PhylogenyTreeBasic *> &listGeneTreePtrs);


//***********************************************************************************
// Genealogical network population label

class GenealogicalNetworkHapPop
{
public:
    GenealogicalNetworkHapPop(const std::vector<std::pair<std::string,std::vector<int> > > &listHapsInPopsIn);
    std::string ConvHapToPop(int hap) const;
    void Dump() const;
    
private:
    std::vector<std::pair<std::string,std::vector<int> > > listHapsInPops;
    std::vector<std::string> mapHapToPop;
};

//***********************************************************************************
// Genealogical network pop-label clade with threshold
class GenealogicalNetworkConvCladeRounded
{
public:
    bool operator<(const GenealogicalNetworkConvCladeRounded &rhs) const;
    void AddPopLabel(const std::string &strPop);
    void AddPopLabel(const std::string &strPop, int numCopies);
    double CalcDist(const GenealogicalNetworkConvCladeRounded &rhs) const;
    void Dump() const;
    
private:
    std::multiset<std::string> setPopLabels;
};


//***********************************************************************************
// Genealogical network pop-label clade

class GenealogicalNetworkConvClade
{
public:
    GenealogicalNetworkConvClade() {}
    bool operator<(const GenealogicalNetworkConvClade &rhs) const;
    void AddHapForPop(const std::string &strPop);
    void AddHapForPop(const std::string &strPop, int ncopy);
    void RoundedTo(int numTargetLins, GenealogicalNetworkConvCladeRounded &cladeRounded) const;
    double CalcDist(const GenealogicalNetworkConvClade &rhs) const;
    static double CalcDist( const std::set< const GenealogicalNetworkConvClade *> &setCladePtrs1, const std::set< const GenealogicalNetworkConvClade *> &setCladePtrs2 );
    int GetNumLins() const;
    void GetPops(std::map<std::string, int> &setPops) const {setPops = mapPopHaps;}
    void GetPops(std::set<std::string> &setPops) const;
    void Dump() const;
    void FilterMinorAlleles(double thresMAF);
    
private:
    void FindIntersect(const GenealogicalNetworkConvClade &rhs, GenealogicalNetworkConvClade &cladeInt ) const;
    void CalcFreqSpec( std::map<std::string, double> &mapFreq ) const;
    double CalcDistFreqSpec( const std::map<std::string, double> &mapFreq1, const std::map<std::string, double> &mapFreq2 ) const;
    
    std::map<std::string, int> mapPopHaps;
};


//***********************************************************************************
// Genealogical network raw clade

class GenealogicalNetworkRawClade
{
public:
    GenealogicalNetworkRawClade(const std::set<int> &setHapsIn);
    GenealogicalNetworkRawClade(const GenealogicalNetworkRawClade &rhs);
    bool operator<(const GenealogicalNetworkRawClade &rhs) const;
    void ConvHapCladeTo( const GenealogicalNetworkHapPop &popInfo, GenealogicalNetworkConvClade &cladeConv ) const;
    void Dump() const;
    
private:
    std::set<int> setHaps;
};

//***********************************************************************************
// Genealogical network clustering

class GenealogicalNetworkCladeCluster
{
public:
    GenealogicalNetworkCladeCluster(const GenealogicalNetworkHapPop &hapPopConvIn);
    void Clustering( const std::map< GenealogicalNetworkRawClade, int> &listGTClades  );
    void Dump() const;
    
private:
    void ConvRawClades( const std::map< GenealogicalNetworkRawClade, int> &listGTClades, std::map< GenealogicalNetworkConvClade, int> &listGTCladesConv  ) const;
    void CoarseClus(double roundBucketSz, double thresDiscard, std::map<const GenealogicalNetworkConvClade *, std::set<const GenealogicalNetworkConvClade *> > setCladePtrsToClus, std::set< std::set<const GenealogicalNetworkConvClade *> > &listClusCladPtrs )  const;
    void FineClus( double maxDistAllowed, std::map<const GenealogicalNetworkConvClade *, std::set<const GenealogicalNetworkConvClade *> > &setCladePtrsToClus );
    void FormConsensus( const std::set<const GenealogicalNetworkConvClade *> &setCladePtrsToClus, std::set<std::string> &cladePop ) const;
    double CalcAveDistBtw( std::map<std::pair<const GenealogicalNetworkConvClade *, const GenealogicalNetworkConvClade *>, double> &mapPairDist, const std::set<const GenealogicalNetworkConvClade *> &setCladePtrsToClus1, const std::set<const GenealogicalNetworkConvClade *> &setCladePtrsToClus2 ) const;
    void FindConsensus( const std::set<const GenealogicalNetworkConvClade *> &setClades, std::set<std::string> &setPops ) const;
    
    
    const GenealogicalNetworkHapPop &hapPopConv;
    std::map<std::set<std::string>, int> mapInfClades;
};


#endif /* GenealogicalNetworkClades_hpp */
