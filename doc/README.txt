CONSEL: a program for assessing the confidence of phylogenetic tree selection
=============================================================================


CONSEL v0.20

Time-stamp: "2011-05-12 16:30:57 shimo"

Hidetoshi Shimodaira


 CHANGE OF AFFILIATION
=======================

Hidetoshi Shimodaira has moved to

Department of Mathematical and Computing Sciences 
Tokyo Institute of Technology 
2-12-1 Ookayama, Meguroku, Tokyo 152-8552, Japan
shimo@is.titech.ac.jp
http://www.is.titech.ac.jp/~shimo/

The former institution is

The Institute of Statistical Mathematics
4-6-7 Minami-Azabu, Minatoku, Tokyo 106-8569, Japan
shimo@ism.ac.jp
http://www.ism.ac.jp/~shimo/


 WHAT IS CONSEL?
=================

CONSEL is a program package consists of small programs written in C
language. It calculates the probability value (i.e., p-value) to
assess the confidence in the selection problem. Although CONSEL is
applicable to any selection problem, it is mainly designed for the
phylogenetic tree selection. CONSEL does not estimate the phylogenetic
tree by itself, but CONSEL does read the output of the other
phylogenetic packages, such as Molphy, PAML, PAUP*, TREE-PUZZLE,
PHYML. CONSEL calculates the p-value using several testing procedures;
the bootstrap probability, the Kishino-Hasegawa test, the
Shimodaira-Hasegawa test, and the weighted Shimodaira-Hasegawa
test. In addition to these conventional tests, CONSEL calculates the
p-value based on the approximately unbiased test using the multi-scale
bootstrap technique. This newly developed method gives less biased
results than the conventional methods.


 PORTABILITY
=============

CONSEL has been developed on a debian Linux box. It should work on
the other UNIX machines with gcc installed. It has been tested on
* Linux / debian (2.2.19, 2.4.18)
* SGI IRIX (6.5.6)
* MS-DOS (Windows 98, 2000, XP)


 MS-DOS USER
=============

For MS-DOS (Windows 98, 2000, XP) users, the executable binary of CONSEL
is available. Please download "cnslb01.exe" for the self-extracting
archive. Read the document "DOSUSER.txt" in it.


 UNIX USER
=============

For UNIX users, the source file of CONSEL is available. Please
download "cnsls01.tgz" for the tar-gzipped archive. Read the
document "UNIXUSER.txt" in it.


 DOCUMENTATION
===============

A tentative user's guide of CONSEL is available at the same directory
as this "README.txt" file. The PDF file is "program.pdf" and its TeX
file is "program.tex".


 CREDITS
=========

The program CONSEL is written by Hidetoshi Shimodaira. The raw-data
reading functions are written by Hiroshi Mine. Some other functions
are borrowed from other packages and modified for CONSEL;

Molphy (Jun Adachi): luinverse.

Meschach library (David E. Stewart): v_sort, skipjunk, and several
ideas on the design of input/output routines.

MT19937 random number generator (Takuji Nishimura and Makoto
Matsumoto) is now used.

I would like to appreciate the following people for making CONSEL.

Masami Hasegawa: for letting me know the biological issues.
Hiroshi Mine: for technical advice on debian linux and Windows.
Joe Felsenstein: brought me to the field of phylogeny.
Brad Efron: for suggestions to develop the method of the AU test.
and all the other people who helped me for the development of CONSEL.


 REFERENCE
===========

(1) Hidetoshi Shimodaira (2000) "Another calculation of the p-value
    for the problem of regions using the scaled bootstrap
    resamplings," Technical Report No. 2000-35, Stanford University.

(2) Hidetoshi Shimodaira (2002) "An approximately unbiased test of
    phylogenetic tree selection," Systematic Biology, 51, 492-508.

(3) Hidetoshi Shimodaira and Masami Hasegawa (2001) "CONSEL: for
    assessing the confidence of phylogenetic tree selection,"
    Bioinformatics, 17, 1246-1247.

(4) Hidetoshi Shimodaira (2004) "Approximately unbiased tests of
    regions using multistep-multiscale bootstrap resampling," 
    Annals of Statistics, 32, 2616-2641.

(5) Hidetoshi Shimodaira (2008) "Testing Regions with Nonsmooth
    Boundaries via Multiscale Bootstrap," Journal of Statistical
    Planning and Inference, 138, 1227-1241.
    http://dx.doi.org/10.1016/j.jspi.2007.04.001

(6) Hidetoshi Shimodaira (2006-2008) "scaleboot: Approximately Unbiased
    P-values via Multiscale Bootstrap" 
    http://cran.r-project.org/web/packages/scaleboot/index.html

The theory is described in the references (1), (2) and (4).  Please
make a citation to the reference (3) for the use of the program
CONSEL.  Recently a new theory which generalizes those for CONSEL is
developed. This new theory is described in (5) and implemented in the
R package of (6).


Hidetoshi Shimodaira
Department of Mathematical and Computing Sciences 
Tokyo Institute of Technology 
http://www.is.titech.ac.jp/~shimo/

----------

This file is part of CONSEL.

CONSEL is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

CONSEL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CONSEL; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

