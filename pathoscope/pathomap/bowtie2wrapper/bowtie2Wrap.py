#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Bowtie2 wrapper performs the alignment using Bowtie2 aligner 
# with appropriate parameters and splits files bigger than 4.3GB.

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
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0,pathoscopedir) 
from pathoscope.utils import samUtils

# ===========================================================
class Bowtie2Options:
	DEFAULT_OPTION = "--very-sensitive-local -k 100 --score-min L,20,1.0" # Another option --score-min L,205,0.0 
	btHome = None # Uses bowtie2 in the system path by default
	verbose = True
	btBin = "bowtie2"
	btIndexer = "bowtie2-build"
	pairedReadFlag = False
	readFile = ""
	readFilePair1 = ""
	readFilePair2 = ""
	outAlignFile = ""
	outDir = "."
	indexDir = "."
	numThreads = 8
	additionalOptions = DEFAULT_OPTION
	refFile = ""
	btIndexPrefix = ""
	def __init__(self): {}

# Runs bowtie2 alignment with the given bowtie2 options
def run_bowtie2(bowtie2Options):
	align_exist=1
	if os.path.exists(bowtie2Options.outAlignFile):
		if bowtie2Options.verbose:
			print "Bowtie2 alignment file already exist: " + bowtie2Options.outAlignFile
		return align_exist
	
	btBinPath = bowtie2Options.btBin
	if bowtie2Options.btHome is not None:
		btBinPath = bowtie2Options.btHome + os.sep + bowtie2Options.btBin
	outAlignFilePath = bowtie2Options.outDir + os.sep + bowtie2Options.outAlignFile
	if bowtie2Options.pairedReadFlag:
		cmd = "%s -x %s -1 %s -2 %s -p %s %s -S %s" % (btBinPath, 
			bowtie2Options.btIndexPrefix, 
			bowtie2Options.readFilePair1, bowtie2Options.readFilePair2,
			bowtie2Options.numThreads, bowtie2Options.additionalOptions,
			outAlignFilePath)
	else:
		cmd = "%s -x %s -U %s -p %s %s -S %s" %	(btBinPath, 
			bowtie2Options.btIndexPrefix, 
			bowtie2Options.readFile,
			bowtie2Options.numThreads, bowtie2Options.additionalOptions,
			outAlignFilePath)
	if bowtie2Options.verbose:
		print cmd
	return os.system(cmd)

# Create bowtie2 Index for the given reference file
def create_bowtie2_index(bowtie2Options):
	index_exist=1
	btIndexPrefixPath = bowtie2Options.indexDir + os.sep + bowtie2Options.btIndexPrefix
	btIndexPath = btIndexPrefixPath + ".1.bt2"
	if os.path.exists(btIndexPath):
		if bowtie2Options.verbose:
			print "Bowtie2 index already exist for: " + btIndexPath
		return index_exist
	btIndexerPath = bowtie2Options.btIndexer
	if bowtie2Options.btHome is not None:
		btIndexerPath = bowtie2Options.btHome + os.sep + bowtie2Options.btIndexer
	cmd = "%s %s %s" % (btIndexerPath, bowtie2Options.refFile, btIndexPrefixPath)
	if bowtie2Options.verbose:
		print cmd
	return os.system(cmd)


# Returns an alignment file with the reads that align in targetAlignFile and
# that do not align in filterAlignFile
def filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile):
	h_readId = find_readsAlignScore(filterAlignFiles)
	with open(targetAlignFile,'r') as in1:
		with open(outAlignFile,'w') as out1:
			for ln in in1:
				if (ln[0] == '@' or ln[0] == '#'):
					out1.write(ln)
					continue
				l = ln.split('\t')
				readId=l[0]
				aScore = samUtils.findSamAlignHiScore(l)
				if (aScore is not None):
					score = h_readId.get(readId, None)
					if (score == None) or (score < aScore):
						#This read is (not/having low scores) in the filterAlignFiles 
						out1.write(ln)
	return outAlignFile
	
# Returns a hash table of all read ids to alignment score 
# from all the given alignment files
# Used by filter_sam() function
def find_readsAlignScore(filterAlignFiles):
	h_readId = {}
	for filterAlignFile in filterAlignFiles:
		with open(filterAlignFile,'r') as in1:
			for ln in in1:
				if (ln[0] == '@' or ln[0] == '#'):
					continue
				l = ln.split('\t')
				readId=l[0]
				aScore = samUtils.findSamAlignHiScore(l)
				if aScore is not None:
					score = h_readId.get(readId,None)
					if (score == None) or (score < aScore):
						h_readId[readId] = aScore
	return h_readId

def extractRead(appendAlignFile, readFile):
	h_readId = {}
	with open(appendAlignFile,'r') as in1:
		with open(readFile,'w') as out1:
			for ln in in1:
				if (ln[0] == '@' or ln[0] == '#'):
					continue
				l = ln.split('\t')
				aScore = samUtils.findSamAlignScore(l)
				if (aScore is None):
					continue
				readId=l[0]
				score = h_readId.get(readId,None)
				if (score is None):
					h_readId[readId] = aScore
					sequence = l[9]
					quality = l[10]
					out1.write('@%s\n%s\n+\n%s\n' % (readId, sequence, quality))

