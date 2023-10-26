#!/bin/bash

cd src
make clean
make
cd ..
cd example
./test.sh
cd ..
