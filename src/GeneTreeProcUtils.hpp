//
//  GeneTreeProcUtils.hpp
//  
//
//  Created by Yufeng Wu on 4/24/19.
//
//

#ifndef GeneTreeProcUtils_hpp
#define GeneTreeProcUtils_hpp

#include <fstream>
#include <iostream>
#include "Utils3.h"

class PhylogenyTreeBasic;
class TaxaMapper;

//*****************************************************************************

bool ReadInGeneTrees( ifstream &inFile, vector<PhylogenyTreeBasic *> &listGeneTreePtrs, TaxaMapper *pTmapperTaxaIds );
void AdjustLabelsForPopInfo( vector<PhylogenyTreeBasic *> &listGeneTreePtrs, const vector<pair<string, vector<int> > > &listPopInfo, TaxaMapper *pTMapper );
void AdjustLabelsTaxaMapper( vector<PhylogenyTreeBasic *> &listGeneTreePtrs, TaxaMapper *pTMapperTree );
void ConsPopHapToPopId( const vector<pair<string, vector<int> > > &listPopInfo, map<int,int> &mapPopRowToPopId  );
void SetNumTreesToPick(int n);
void PickTrees(vector<PhylogenyTreeBasic *> &listGeneTreePtrsToProc);

#endif /* GeneTreeProcUtils_hpp */
