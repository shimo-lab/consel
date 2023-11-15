#!/bin/bash

mkdir -p results
TESTDIR=~/consel/test_data
BINLOC=~/consel/src
PATH=..:../src:../bin:$PATH

seqmt --puzzle $TESTDIR/RAxML_perSiteLLs_ycf2_ge.trees.sitelh $TESTDIR/ycf2_ge_CONSEL.mt
gprof $BINLOC/seqmt gmon.out > $TESTDIR/results/prof4_seqmt.txt
makermt $TESTDIR/ycf2_ge_CONSEL.mt
gprof $BINLOC/makermt gmon.out > $TESTDIR/results/prof4_makermt.txt
consel $TESTDIR/ycf2_ge_CONSEL.rmt
gprof $BINLOC/consel gmon.out > $TESTDIR/results/prof4_consel.txt
catpv $TESTDIR/ycf2_ge_CONSEL.pv > $TESTDIR/ycf2_ge_CONSEL.consel