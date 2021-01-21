#!bin/sh

source ~/.bashrc

./benchmark.ex -f input.dat -startID 1 -endID 235  2>&1 |tee -a prob1-235.log
