#! /bin/sh

PATH=..:../src:../bin:$PATH

makermt --paml mam15
consel mam15
catpv mam15 > mam15-pv.txt
cat mam15-pv.txt

