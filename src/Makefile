CC = g++
LIBPATH =

#CFLAGS = -O3 -Wall -g -D HAVE_PTHREADS
CFLAGS = -O3 -std=c++11 -g -D HAVE_PTHREADS

all: gtmix


Utils.o : Utils.cpp  Utils.h
	$(CC)  $(CFLAGS) -c Utils.cpp

Utils2.o : Utils2.cpp  Utils2.h Utils.h
	$(CC)  $(CFLAGS) -c Utils2.cpp

Utils3.o : Utils3.cpp  Utils3.h Utils.h
	$(CC)  $(CFLAGS) -c Utils3.cpp

Utils4.o : Utils4.cpp  Utils4.h Utils.h
	$(CC)  $(CFLAGS) -c Utils4.cpp

UtilsNumerical.o : UtilsNumerical.cpp  UtilsNumerical.h
	$(CC)  $(CFLAGS) -c UtilsNumerical.cpp

#Utils-basic.o : Utils-basic.cpp  Utils-basic.h
#	$(CC)  $(CFLAGS) -c Utils-basic.cpp

TreeBuilder.o : TreeBuilder.cpp  TreeBuilder.h
	$(CC)  $(CFLAGS) -c TreeBuilder.cpp

#GraphUtils.o : GraphUtils.cpp  GraphUtils.h
#	$(CC)  $(CFLAGS) -c GraphUtils.cpp

UnWeightedGraph.o : UnWeightedGraph.cpp  UnWeightedGraph.h
	$(CC)  $(CFLAGS) -c UnWeightedGraph.cpp

#SNPUtils.o : SNPUtils.cpp SNPUtils.h
#	$(CC)  $(CFLAGS) -c SNPUtils.cpp

BioSequenceMatrix.o : BioSequenceMatrix.cpp BioSequenceMatrix.h
	$(CC)  $(CFLAGS) -c BioSequenceMatrix.cpp

#BinaryMatrix.o : BinaryMatrix.cpp BinaryMatrix.h
#	$(CC)  $(CFLAGS) -c BinaryMatrix.cpp

#MultiBinaryMatrices.o : MultiBinaryMatrices.cpp MultiBinaryMatrices.h
#	$(CC)  $(CFLAGS) -c MultiBinaryMatrices.cpp

#GTHapAnalyzer.o : GTHapAnalyzer.cpp GTHapAnalyzer.hpp
#	$(CC)  $(CFLAGS) -c GTHapAnalyzer.cpp

MarginalTree.o : MarginalTree.cpp MarginalTree.h
	$(CC)  $(CFLAGS) -c MarginalTree.cpp

PhylogenyTreeBasic.o : PhylogenyTreeBasic.cpp PhylogenyTreeBasic.h
	$(CC)  $(CFLAGS) -c PhylogenyTreeBasic.cpp

GeneTreeProcUtils.o : GeneTreeProcUtils.cpp GeneTreeProcUtils.hpp
	$(CC)  $(CFLAGS) -c GeneTreeProcUtils.cpp

DeepCoalescence.o : DeepCoalescence.cpp DeepCoalescence.h
	$(CC)  $(CFLAGS) -c DeepCoalescence.cpp

GeneSpeciesTreeProb.o : GeneSpeciesTreeProb.cpp GeneSpeciesTreeProb.h
	$(CC)  $(CFLAGS) -c GeneSpeciesTreeProb.cpp

#ApproxGeneTreeProb.o : ApproxGeneTreeProb.cpp ApproxGeneTreeProb.h
#	$(CC)  $(CFLAGS) -c ApproxGeneTreeProb.cpp

GenealogicalNetwork.o : GenealogicalNetwork.cpp GenealogicalNetwork.h
	$(CC)  $(CFLAGS) -c GenealogicalNetwork.cpp

GenealogicalNetworkHeuSearch.o : GenealogicalNetworkHeuSearch.cpp GenealogicalNetworkHeuSearch.h
	$(CC)  $(CFLAGS) -c GenealogicalNetworkHeuSearch.cpp

ApproxGeneTreeProb2.o : ApproxGeneTreeProb2.cpp ApproxGeneTreeProb2.h
	$(CC)  $(CFLAGS) -c ApproxGeneTreeProb2.cpp

#GenealogicalNetworkProb.o : GenealogicalNetworkProb.cpp GenealogicalNetworkProb.h
#	$(CC)  $(CFLAGS) -c GenealogicalNetworkProb.cpp

GenealogicalNetworkProbStore.o : GenealogicalNetworkProbStore.cpp GenealogicalNetworkProbStore.h
	$(CC)  $(CFLAGS) -c GenealogicalNetworkProbStore.cpp

#ApproxGeneTreeProb.o : ApproxGeneTreeProb.cpp ApproxGeneTreeProb.hpp
#	$(CC)  $(CFLAGS) -c ApproxGeneTreeProb.cpp

AGFastCoalBuilder.o : AGFastCoalBuilder.cpp AGFastCoalBuilder.hpp
	$(CC)  $(CFLAGS) -c AGFastCoalBuilder.cpp

AGGeneTreeQuickCoal.o : AGGeneTreeQuickCoal.cpp AGGeneTreeQuickCoal.hpp
	$(CC)  $(CFLAGS) -c AGGeneTreeQuickCoal.cpp

GenealogicalNetworkClades.o : GenealogicalNetworkClades.cpp GenealogicalNetworkClades.hpp
	$(CC)  $(CFLAGS) -c GenealogicalNetworkClades.cpp

AGFastLenOpt.o : AGFastLenOpt.cpp AGFastLenOpt.hpp
	$(CC)  $(CFLAGS) -c AGFastLenOpt.cpp

AGGeneTreeProb.o : AGGeneTreeProb.cpp AGGeneTreeProb.hpp
	$(CC)  $(CFLAGS) -c AGGeneTreeProb.cpp

AGmain.o : AGmain.cpp
	$(CC)  $(CFLAGS) -c AGmain.cpp

# Now final link

gtmix: Utils.o Utils2.o  Utils3.o Utils4.o UtilsNumerical.o TreeBuilder.o UnWeightedGraph.o PhylogenyTreeBasic.o  MarginalTree.o GeneTreeProcUtils.o \
   GeneSpeciesTreeProb.o GenealogicalNetwork.o GenealogicalNetworkHeuSearch.o GenealogicalNetworkProbStore.o  DeepCoalescence.o  \
   GenealogicalNetworkClades.o ApproxGeneTreeProb2.o AGFastCoalBuilder.o AGGeneTreeQuickCoal.o AGFastLenOpt.o AGGeneTreeProb.o AGmain.o
	$(CC) -o gtmix2 Utils.o Utils2.o  Utils3.o Utils4.o UtilsNumerical.o TreeBuilder.o UnWeightedGraph.o  PhylogenyTreeBasic.o MarginalTree.o GeneTreeProcUtils.o \
  GeneSpeciesTreeProb.o  GenealogicalNetwork.o GenealogicalNetworkHeuSearch.o GenealogicalNetworkProbStore.o  DeepCoalescence.o \
   GenealogicalNetworkClades.o ApproxGeneTreeProb2.o AGFastCoalBuilder.o AGGeneTreeQuickCoal.o AGFastLenOpt.o AGGeneTreeProb.o AGmain.o \
        -lm -lpthread -std=c++11
