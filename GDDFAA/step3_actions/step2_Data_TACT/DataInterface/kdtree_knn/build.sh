#!/bin/bash
rm -rdf ./build > /dev/null 2>&1
mkdir build 
cd build
CXX=g++ cmake ..
make
cd ..
