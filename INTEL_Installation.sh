#!/bin/bash
echo '##############################################'
echo '#                                            #'
echo '#               Intel Compiler               #'
echo '# Compiling AlignmentByNumbers initiation... #'
echo '#                                            #'


cwd=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/Code/"
cd $cwd


icpc -Wall -O3 -std=c++11 -o AlignmentByNumbers *.cpp *.c *.cc
mv AlignmentByNumbers ../Exec

echo '#                                            #'
echo '# Compiling AlignmentByNumbers completed...  #'
echo '#                                            #'
echo '##############################################'

