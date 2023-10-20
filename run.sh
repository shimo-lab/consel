#!/bin/bash

cd src
make
cd ..
cd example
./test.sh
cd ..
