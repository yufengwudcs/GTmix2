//
//  AGGeneTreeQuickCoal.hpp
//  
//
//  Created by Yufeng Wu on 2/16/24.
//

#ifndef AGGeneTreeQuickCoal_hpp
#define AGGeneTreeQuickCoal_hpp

#include <vector>
#include <set>
#include <map>

class GenealogicalNetwork;
class GenealogicalNetworkNode;
class GenealogicalNetworkBranch;
class PhylogenyTreeBasic;
class TreeNode;
class TaxaMapper;


//***********************************************************************************
// Placing coalescents in the given tree onto network using parsimony

class AGGeneTreeQuickCoal2
{
public:
    AGGeneTreeQuickCoal2( GenealogicalNetwork &agIn );
    double CalcGeneTreeProbHeu(PhylogenyTreeBasic &phTree);
    double CalcGeneTreeProbHeu2(PhylogenyTreeBasic &phTree);
    double CalcGeneTreeProbHeu3(PhylogenyTreeBasic &phTree);
    double CalcGeneTreeProbHeu4(PhylogenyTreeBasic &phTree);
    
private:
    void Init();
    bool IsAGNodeAncestralTo(GenealogicalNetworkNode *pnDesc, GenealogicalNetworkNode *pnAnc) const;
    bool IsAGNodeAncestralTo(const std::vector<GenealogicalNetworkNode *> &listDesc, GenealogicalNetworkNode *pnAnc) const;
    void GetAllAGBrsBottomUp( std::vector<GenealogicalNetworkBranch *> &listBrs ) const;
    void UpdateExpLinInfoBrTop(GenealogicalNetworkBranch *pBr, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop);
    void UpdateDistLinInfoBrTop(GenealogicalNetworkBranch *pBr, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrExpNumLinsTop);
    double CalcProbGTLinNotCoalBetween( GenealogicalNetworkBranch *pBrDesc, GenealogicalNetworkBranch *pBrAnc, std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop );
    double CalcProbGTLinNotCoalBetween2( GenealogicalNetworkBranch *pBrDesc, GenealogicalNetworkBranch *pBrAnc, std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop );
    double CalcProbNotCoalescedFor(int m, int n) const;
    double CalcCoalProbAtBranch(GenealogicalNetworkBranch *pBr, GenealogicalNetworkBranch *pBrChild1, GenealogicalNetworkBranch *pBrChild2, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop) const;
    double CalcCoalProbAtBranch2(GenealogicalNetworkBranch *pBr, GenealogicalNetworkBranch *pBrChild1, GenealogicalNetworkBranch *pBrChild2, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop) const;
    double CalcCoalCoeff(int m, int k) const;
    
    // calc expected number of branches
    void InitExpNumLins(PhylogenyTreeBasic &phTree, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop);
    // calc distribution of number of branches
    void InitDistNumLins(PhylogenyTreeBasic &phTree, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop);
    // init GT lin not coalseced
    void InitGTLinNotCoals(std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsBottom, std::map< GenealogicalNetworkBranch *, double> &mapAGBrExpNumLinsTop);
    // init GT lin not coalseced
    void InitGTLinNotCoals2(std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > &mapGTLinNotCoalFromDescToAncInAG, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsBottom, std::map< GenealogicalNetworkBranch *, std::vector<double> > &mapAGBrDistNumLinsTop);
    
    //
    GenealogicalNetwork &ag;
    std::set<int> setAllPops;
    
    // info about ag
    std::vector<GenealogicalNetworkNode *> listNodesAG;
    std::map<GenealogicalNetworkNode *, std::set<GenealogicalNetworkNode *> > listDescNodesAG;      // for each node in AG, the set of descendant nodes
    std::map<GenealogicalNetworkNode *, std::set<GenealogicalNetworkNode *> > listAncesAGNodes;     // set of ancestral nodes for each node in AG
    // keep track population subsets that are reachable at each network node
    std::map<GenealogicalNetworkNode *, std::set<int> > mapPopsAtNodes;
    
    
    // for each AG branch, what is the expected number of branches at bottom and top of this branch
    //std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsBottom;
    //std::map< GenealogicalNetworkBranch *, double> mapAGBrExpNumLinsTop;
    
    // prob of a GT lineage not involved in coalescents from one point of AG to an ancestral position of AG
    //std::map<GenealogicalNetworkBranch *, std::map<GenealogicalNetworkBranch *,double> > mapGTLinNotCoalFromDescToAncInAG;
};

#endif /* AGGeneTreeQuickCoal_hpp */
