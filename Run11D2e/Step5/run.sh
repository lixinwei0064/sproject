#!/bin/bash
for i in {0,1};do
    for j in {0..3};do
        root -b -q -x 'combine.C('$i','$j')'
    done
done