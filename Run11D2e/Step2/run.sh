#!/bin/bash
for i in 3;do
        for j in 1;do
            root -b -q -x 'run.C('$i','$j','$RANDOM','$RANDOM')'
        done
done