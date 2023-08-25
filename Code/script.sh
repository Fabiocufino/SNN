#!/bin/bash

for ((c=1; c<=24; c++))
do
    nohup root -l -b < ./starting_parameters/start_$c.cmd >&./out/log_$c.txt 2>&1 & 
done
