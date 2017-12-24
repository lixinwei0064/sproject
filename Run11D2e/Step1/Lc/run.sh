#!/bin/bash

for i in {0..3}; do   
    root -b -q -x 'run.C('$i','$RANDOM')'
done