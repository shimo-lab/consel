# consel
CONSEL: for assessing the confidence of phylogenetic tree selection

Hidetoshi Shimodaira

http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/

WHAT IS CONSEL?

CONSEL is a program package consists of small programs written in C language. It calculates the probability value (i.e., p-value) to assess the confidence in the selection problem. Although CONSEL is applicable to any selection problem, it is mainly designed for the phylogenetic tree selection. CONSEL does not estimate the phylogenetic tree by itself, but CONSEL does read the output of the other phylogenetic packages, such as Molphy, PAML, PAUP*, TREE-PUZZLE, and PhyML. CONSEL calculates the p-value using several testing procedures; the bootstrap probability, the Kishino-Hasegawa test, the Shimodaira-Hasegawa test, and the weighted Shimodaira-Hasegawa test. In addition to these conventional tests, CONSEL calculates the p-value based on the approximately unbiased test using the multi-scale bootstrap technique. This newly developed method gives less biased results than the conventional methods.

## Parallelization Options

Specify the number of threads used by makermt.c by using the flag OMP_NUM_THREADS=# before the binary file. 
(example: OMP_NUM_THREADS=8 makermt --paml mam15) for 8 threads to be used.