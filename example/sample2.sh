#! /bin/sh

PATH=..:../src:../bin:$PATH

treeass --outgroup 6 mam15 > mam15-tree.log
consel -a mam15 mam15 mam15e
catpv mam15e > mam15e-pv.txt
cat mam15e-pv.txt
