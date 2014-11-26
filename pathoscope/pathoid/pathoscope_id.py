#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Gets alignment file (currently support sam or BLAST-m8 format (.bl8)) and runs EM algorithm.
# Outputs the pathogens rank in the sample as a report that can be opened in Excel.
# Optionally outputs an updated alignment file (sam/bl8)

#usage information: pathoscope.py -h

#	Pathoscope - Predicts strains of genomes in Nextgen seq alignment file (sam/bl8)
#	Copyright (C) 2013  Johnson Lab - Boston University
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, sys
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir) 

import argparse
from pathoscope.pathoid import PathoID
from time import time

# ===========================================================
# main ()
parser = argparse.ArgumentParser(description="Pathoscope")
parser.add_argument('-o', action='store_true', default=False, dest='out_matrix',
					help='Output alignment matrix')
parser.add_argument('-noUpdatedAlignFile', action='store_true', default=False, dest='noalign',
					help='Do not generate an updated alignment file')
parser.add_argument('-verbose', action='store_true', default=False, dest='verbose',
					help='Prints verbose text while running')
parser.add_argument('-s', action='store', default=0.01, type=float,
					dest='score_cutoff', help='Score Cutoff')
parser.add_argument('-emEpsilon', action='store', default=1e-7, type=float,
					dest='emEpsilon', help='EM Algorithm Epsilon cutoff')
parser.add_argument('-maxIter', action='store', default=50, type=int,
					dest='maxIter', help='EM Algorithm maximum iterations')
parser.add_argument('-e', action='store', default='testset', dest='exp_tag',
					help='Experiment tag')
parser.add_argument('-outdir', action='store', default='.', dest='outdir',
					help='Output Directory')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-t', action='store', default='sam', dest='ali_format',
					help='Alignment Format: sam/bl8/gnu-sam')
parser.add_argument('-f', action='store', dest='ali_file', required=True,
					help='Alignment file path')

inputArgs = parser.parse_args()
start = time();
PathoID.pathoscope_reassign(inputArgs.out_matrix, inputArgs.verbose, 
	inputArgs.score_cutoff, inputArgs.exp_tag, inputArgs.ali_format, inputArgs.ali_file,
	inputArgs.outdir, inputArgs.emEpsilon, inputArgs.maxIter, not(inputArgs.noalign))
elapsed = time() - start;
if inputArgs.verbose:
	print "EM Elapsed Time: %d" % (elapsed)
