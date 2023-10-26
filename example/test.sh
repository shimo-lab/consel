#! /bin/sh

mkdir -p tmp

PATH=..:../src:../bin:$PATH
ERR=error.log
echo $0 >> $ERR
date >> $ERR

TMP=tmp
SRC=result

FILE=test1
randrep -m mu1 $TMP/$FILE
consel $TMP/$FILE
catpv $TMP/$FILE > $TMP/$FILE.txt
diff $TMP/$FILE.txt $SRC/$FILE.txt >> $ERR

gprof ~/Path/To/Executeable/consel gmon.out > prof1.txt

FILE=test2
makermt -p k10s --paml mam15 $TMP/$FILE
consel $TMP/$FILE
catpv $TMP/$FILE > $TMP/$FILE.txt
diff $TMP/$FILE.txt $SRC/$FILE.txt >> $ERR

gprof ~/Path/To/Executeable/consel gmon.out > prof2.txt

PREV=test2
FILE=test3
treeass mam15 $TMP/$FILE > $TMP/$FILE.log
consel -a $TMP/$FILE $TMP/$PREV $TMP/$FILE
catpv $TMP/$FILE > $TMP/$FILE.txt
diff $TMP/$FILE.txt $SRC/$FILE.txt >> $ERR

gprof ~/Path/To/Executeable/consel gmon.out > prof3.txt

date >> $ERR
