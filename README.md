# GTmix2
##Inference of admixture graph (AG) based on coalescent theory

GTmix2 is the more efficient version of the original GTmix program (https://github.com/yufengwudcs/GTmix).

## How to build GTmix2?

GTmix2 is simple to build. Once you download the zip file from this GitHub repository. All you need is to (i) open a terminal to go to the unzipped folder, (ii) change directory to the "src" folder and  (iii) type "make". That is,

cd src
make

The executable gtmix2 is located within the src folder.

## How to run GTmix2?

Start by test whether GTmix2 runs on your computer as follows. From the src folder, type:

./gtmix2 -n 1 -P listPops4p1a.txt test-trees.nwk

Explanation: 
"-n 1":                you want to reconstruct an AG with one admixture. 
"-P listPops4p1a.txt:  specifiy the population definiton file, which specifies which haplotype (numerically numbered by default) belongs to which population. 
"test-trees.nwk"       the file containing the genealogical trees.

If it runs, you should see:

*** GTmix2 ver. 1.0.0.0, June 12, 2025 ***

Fix number of admixture nodes to: 1
Number of gene trees to use: 11
Scaffold tree: ((A:1.0000,C:1.0000):0.6364,(B:1.0000,D:1.0000):0.3636)
Maximum number of gene alleles: 4
Initial network probability: -25.7869
Number of neighboring nets to evaluate: 50
****** Highest log-probabiliyt of optimized network (searching over network space): -23.5208
Time needed to find the optimal network: 0
Optimal network: output to file: optimal-network.gml
List of marginal trees in the optimal network:
[0.5] ((B:0.0020,D:0.0010):0.4000,(A:1.0000,C:1.0000):0.7500)
[0.5] (D:0.4010,(C:1.0000,(A:0.5000,B:0.0020):0.5000):0.7500)
Admixture population: B
Cache probability computation: [0] cache hits among total 0 queries, num of storing operations: 0
AGProcessedNetsDepot: total number of novel networks processed: 0, number of skipped networks: 0
ApproxGTPCache2: Num of listdbVals entries: 14
ApproxGTPCache2: Num of listCachedCoalProb entries: 0
Elapsed time = 0 seconds.

Pay attention to the following:
List of marginal trees in the optimal network:
[0.5] ((B:0.0020,D:0.0010):0.4000,(A:1.0000,C:1.0000):0.7500)The


## How does GTmix2 work?
I haven't written up the methodologies of GTmix2 yet. To get some ideas, I recommend to check out the original GTmix repository (https://github.com/yufengwudcs/GTmix). 

## How to cite this work?
For now, cite the original GTmix paper:

Inference of Population Admixture Network from Local Gene Genealogies: a Coalescent-based Maximum Likelihood Approach, Yufeng Wu, Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i326–i334, https://doi.org/10.1093/bioinformatics/btaa465.
