#!/bin/bash

mkdir -p results
TESTDIR=~/consel/test_data
BINLOC=~/consel/src
PATH=..:../src:../bin:$PATH

seqmt --puzzle $TESTDIR/RAxML_perSiteLLs_atpF_ge.trees.sitelh $TESTDIR/atpF_ge_CONSEL.mt
gprof $BINLOC/seqmt gmon.out > $TESTDIR/results/prof1_seqmt.txt
makermt $TESTDIR/atpF_ge_CONSEL.mt
gprof $BINLOC/makermt gmon.out > $TESTDIR/results/prof1_makermt.txt
consel $TESTDIR/atpF_ge_CONSEL.rmt
gprof $BINLOC/consel gmon.out > $TESTDIR/results/prof1_consel.txt
catpv $TESTDIR/atpF_ge_CONSEL.pv > $TESTDIR/atpF_ge_CONSEL.consel
