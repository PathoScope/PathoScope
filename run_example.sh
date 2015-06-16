#!/bin/bash

echo "python pathoscope/pathoscope2.py  ID -alignFile example/MAP_3852_align.sam -expTag 3852 -outDir results"
python pathoscope/pathoscope2.py  ID -alignFile example/MAP_3852_align.sam -expTag 3852 -outDir results

echo "python pathoscope/pathoscope2.py  REP -samfile results/updated_MAP_3852_align.sam -outDir results"
python pathoscope/pathoscope2.py  REP -samfile results/updated_MAP_3852_align.sam -outDir results

echo "done"
