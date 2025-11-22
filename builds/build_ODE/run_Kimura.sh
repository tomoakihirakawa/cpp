#!/bin/bash

for i in {0..14}; do
    ./MDOF_beam_vibration 15 $i
done
