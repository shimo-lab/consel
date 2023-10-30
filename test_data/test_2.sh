#!/bin/bash

mkdir -p results
TESTDIR=~/GroupProject/consel/test_data
BINLOC=~/GroupProject/consel/src
PATH=..:../src:../bin:$PATH

seqmt --puzzle $TESTDIR/RAxML_perSiteLLs_clpP_ge.trees.sitelh $TESTDIR/clpP_ge_CONSEL.mt
gprof $BINLOC/seqmt gmon.out > $TESTDIR/results/prof2_seqmt.txt
makermt $TESTDIR/clpP_ge_CONSEL.mt
gprof $BINLOC/makermt gmon.out > $TESTDIR/results/prof2_makermt.txt
consel $TESTDIR/clpP_ge_CONSEL.rmt
gprof $BINLOC/consel gmon.out > $TESTDIR/results/prof2_consel.txt
catpv $TESTDIR/clpP_ge_CONSEL.pv > $TESTDIR/clpP_ge_CONSEL.consel
