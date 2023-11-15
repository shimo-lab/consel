#!/bin/bash

mkdir -p results
TESTDIR=~/consel/test_data
BINLOC=~/consel/src
PATH=..:../src:../bin:$PATH

seqmt --puzzle $TESTDIR/RAxML_perSiteLLs_ycf1_ge.trees.sitelh $TESTDIR/ycf1_ge_CONSEL.mt
gprof $BINLOC/seqmt gmon.out > $TESTDIR/results/prof3_seqmt.txt
makermt $TESTDIR/ycf1_ge_CONSEL.mt
gprof $BINLOC/makermt gmon.out > $TESTDIR/results/prof3_makermt.txt
consel $TESTDIR/ycf1_ge_CONSEL.rmt
gprof $BINLOC/consel gmon.out > $TESTDIR/results/prof3_consel.txt
catpv $TESTDIR/ycf1_ge_CONSEL.pv > $TESTDIR/ycf1_ge_CONSEL.consel
