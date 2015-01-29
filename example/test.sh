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

FILE=test2
makermt -p k10s --paml mam15 $TMP/$FILE
consel $TMP/$FILE
catpv $TMP/$FILE > $TMP/$FILE.txt
diff $TMP/$FILE.txt $SRC/$FILE.txt >> $ERR

PREV=test2
FILE=test3
treeass mam15 $TMP/$FILE > $TMP/$FILE.log
consel -a $TMP/$FILE $TMP/$PREV $TMP/$FILE
catpv $TMP/$FILE > $TMP/$FILE.txt
diff $TMP/$FILE.txt $SRC/$FILE.txt >> $ERR

date >> $ERR
