#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# PathoMap performs the alignment through wrappers for each type of aligners.

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

import os, sys, math, shutil
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir) 

from pathoscope.pathomap.bowtie2wrapper import bowtie2Wrap
from pathoscope.utils import seqParse

# ===========================================================
class PathoMapOptions:
	MAX_REF_FILE_SIZE = 4.3e9
	verbose = False
	outDir = "."
	indexDir = "."
	numThreads = 8
	outAlignFile = "outalign.sam"
	inReadFile = ""
	inReadFilePair1 = ""
	inReadFilePair2 = ""
	targetRefFiles = []
	filterRefFiles = []
	targetIndexPrefixes = []
	filterIndexPrefixes = []
	targetAlignFiles = []
	filterAlignFiles = []
	targetAlignParameters = None
	filterAlignParameters = None
	btHome = None
	exp_tag = ""

# Main entry function to PathoMap that does all the processing
def processPathoMap(pathoMapOptions):
	procPathoMapOptions = copyPathoMapOptions(pathoMapOptions)
	# Splitting reference files if bigger than MAX_REF_FILE_SIZE
	ptargetRefFiles = []
	for filePath in pathoMapOptions.targetRefFiles:
		if pathoMapOptions.verbose:
			print "Checking whether the file: " + filePath + " needs to be split"
		files = splitCheck(filePath, pathoMapOptions.MAX_REF_FILE_SIZE);
		for f in files:
			ptargetRefFiles.append(f)
	procPathoMapOptions.targetRefFiles = ptargetRefFiles
	pfilterRefFiles = []
	for filePath in pathoMapOptions.filterRefFiles:
		if pathoMapOptions.verbose:
			print "Checking whether the file: " + filePath + " needs to be split"
		files = splitCheck(filePath, pathoMapOptions.MAX_REF_FILE_SIZE);
		for f in files:
			pfilterRefFiles.append(f)
	procPathoMapOptions.filterRefFiles = pfilterRefFiles
	
	# Creating Index if it does not exist
	bowtie2Options = bowtie2Wrap.Bowtie2Options()
	bowtie2Options.verbose = procPathoMapOptions.verbose
	bowtie2Options.btHome = procPathoMapOptions.btHome
	bowtie2Options.indexDir = procPathoMapOptions.indexDir
	for filePath in ptargetRefFiles:
		bowtie2Options.refFile = filePath
		(_, tail) = os.path.split(filePath)
		(base, _) = os.path.splitext(tail)
		bowtie2Options.btIndexPrefix = base
		if pathoMapOptions.verbose:
			print "Creating bowtie2 index for: " + filePath
		bowtie2Wrap.create_bowtie2_index(bowtie2Options)
		procPathoMapOptions.targetIndexPrefixes.append(base)
	for filePath in pfilterRefFiles:
		bowtie2Options.refFile = filePath
		(_, tail) = os.path.split(filePath)
		(base, _) = os.path.splitext(tail)
		bowtie2Options.btIndexPrefix = base
		if pathoMapOptions.verbose:
			print "Creating bowtie2 index for: " + filePath
		bowtie2Wrap.create_bowtie2_index(bowtie2Options)
		procPathoMapOptions.filterIndexPrefixes.append(base)
	
	# Creating the Alignment file
	bowtie2Options = bowtie2Wrap.Bowtie2Options()
	bowtie2Options.verbose = procPathoMapOptions.verbose
	bowtie2Options.btHome = procPathoMapOptions.btHome
	bowtie2Options.numThreads = procPathoMapOptions.numThreads
	bowtie2Options.outDir = procPathoMapOptions.outDir
	bowtie2Options.indexDir = procPathoMapOptions.indexDir
	bowtie2Options.readFile = procPathoMapOptions.inReadFile
	bowtie2Options.readFilePair1 = procPathoMapOptions.inReadFilePair1
	bowtie2Options.readFilePair2 = procPathoMapOptions.inReadFilePair2
	if (len(procPathoMapOptions.inReadFilePair1)>0 and
		len(procPathoMapOptions.inReadFilePair2)>0):
		bowtie2Options.pairedReadFlag = True
	if procPathoMapOptions.targetAlignParameters is not None:
		bowtie2Options.additionalOptions = procPathoMapOptions.targetAlignParameters
	for indexPrefix in procPathoMapOptions.targetIndexPrefixes:
		bowtie2Options.btIndexPrefix = procPathoMapOptions.indexDir + os.sep + indexPrefix
		bowtie2Options.outAlignFile = procPathoMapOptions.exp_tag + indexPrefix + ".sam"
		if pathoMapOptions.verbose:
			print "Creating bowtie2 alignment: " + bowtie2Options.outAlignFile
		bowtie2Wrap.run_bowtie2(bowtie2Options)
		procPathoMapOptions.targetAlignFiles.append(procPathoMapOptions.outDir + os.sep + 
			bowtie2Options.outAlignFile)

	# Appending the Alignment files and Filtering
	if len(procPathoMapOptions.targetAlignFiles) > 1:
		appendAlignFile = procPathoMapOptions.outDir + os.sep + procPathoMapOptions.exp_tag + "appendAlign.sam"
		if pathoMapOptions.verbose:
			print "Appending alignment files to: " + appendAlignFile
		append_sam_file(appendAlignFile, procPathoMapOptions.targetAlignFiles)
	else:
		appendAlignFile = procPathoMapOptions.targetAlignFiles[0]

	if len(procPathoMapOptions.filterIndexPrefixes) > 0: 
		bowtie2Options.readFile = procPathoMapOptions.outDir + os.sep + procPathoMapOptions.exp_tag + "appendAlign.fq"
		bowtie2Options.readFilePair1 = ""
		bowtie2Options.readFilePair2 = ""
		bowtie2Options.pairedReadFlag = False
		if procPathoMapOptions.filterAlignParameters is not None:
			bowtie2Options.additionalOptions = procPathoMapOptions.filterAlignParameters
		bowtie2Wrap.extractRead(appendAlignFile, bowtie2Options.readFile)
		for indexPrefix in procPathoMapOptions.filterIndexPrefixes:
			bowtie2Options.btIndexPrefix = procPathoMapOptions.indexDir + os.sep + indexPrefix
			bowtie2Options.outAlignFile = procPathoMapOptions.exp_tag + indexPrefix + ".sam"
			if pathoMapOptions.verbose:
				print "Creating bowtie2 alignment: " + bowtie2Options.outAlignFile
			bowtie2Wrap.run_bowtie2(bowtie2Options)
			procPathoMapOptions.filterAlignFiles.append(procPathoMapOptions.outDir + os.sep + 
				bowtie2Options.outAlignFile)
	# Filtering the Alignment file
	outAlignFile = procPathoMapOptions.outDir + os.sep + procPathoMapOptions.outAlignFile
	if pathoMapOptions.verbose:
		print "Filtering and creating the alignment: " + outAlignFile
	if len(procPathoMapOptions.filterAlignFiles) > 0:
		filter_alignment(appendAlignFile, procPathoMapOptions.filterAlignFiles, outAlignFile)
	elif ((len(procPathoMapOptions.targetAlignFiles) > 1) or \
		(len(procPathoMapOptions.targetIndexPrefixes) > 0)):
		os.rename(appendAlignFile, outAlignFile)
	else: # Input appendAlignFile provided by user, hence make a copy for outAlignFile
		shutil.copy(appendAlignFile, outAlignFile)
	

# Make a copy of core pathoMapOptions
def copyPathoMapOptions(pathoMapOptions):
	procPathoMapOptions = PathoMapOptions()
	procPathoMapOptions.verbose = pathoMapOptions.verbose
	procPathoMapOptions.outDir = pathoMapOptions.outDir
	procPathoMapOptions.indexDir = pathoMapOptions.indexDir
	procPathoMapOptions.outAlignFile = pathoMapOptions.outAlignFile
	procPathoMapOptions.inReadFile = pathoMapOptions.inReadFile
	procPathoMapOptions.inReadFilePair1 = pathoMapOptions.inReadFilePair1
	procPathoMapOptions.inReadFilePair2 = pathoMapOptions.inReadFilePair2
	procPathoMapOptions.targetRefFiles = pathoMapOptions.targetRefFiles
	procPathoMapOptions.filterRefFiles = pathoMapOptions.filterRefFiles
	procPathoMapOptions.targetIndexPrefixes = pathoMapOptions.targetIndexPrefixes
	procPathoMapOptions.filterIndexPrefixes = pathoMapOptions.filterIndexPrefixes
	procPathoMapOptions.targetAlignFiles = pathoMapOptions.targetAlignFiles
	procPathoMapOptions.filterAlignFiles = pathoMapOptions.filterAlignFiles
	procPathoMapOptions.targetAlignParameters = pathoMapOptions.targetAlignParameters
	procPathoMapOptions.filterAlignParameters = pathoMapOptions.filterAlignParameters
	procPathoMapOptions.btHome = pathoMapOptions.btHome
	procPathoMapOptions.exp_tag = pathoMapOptions.exp_tag
	return procPathoMapOptions

# If the given file size is greater than maxSize, then it splits
# and returns a list of split file paths where each file is less than maxSize
def splitCheck(filePath, maxSize):
	files = []
	fileSize = os.stat(filePath).st_size
	nSplit = 1
	if (fileSize > maxSize):
		nSplit = int(math.ceil(1.0*fileSize/float(maxSize)))
	if nSplit==1:
		files.append(filePath)
		return files
	(base, ext) = os.path.splitext(filePath)
	#check if we have already done this splitting
	for i in range(nSplit):
		fiPath=base+'_'+str(i)+ext
		splitReq=False
		if not os.path.exists(fiPath):
			splitReq=True
			break
	fps = []
	for i in range(nSplit):
		fiPath=base+'_'+str(i)+ext
		files.append(fiPath)
		if splitReq:
			fps.append(open(fiPath,'w'))
	if splitReq:
		with open(filePath,'r') as fp:
			j=0
			if ext=='.fq':
				for r in seqParse.parse(fp,'fastq'):
					fps[j%nSplit].write('>%s %s\n%s\n%s\n' % (r.id, r.description, r.seq, r.qual))
					j+=1
			else:
				for r in seqParse.parse(fp,'fasta'):
					fps[j%nSplit].write('>%s %s\n%s\n' % (r.id, r.description, r.seq))
					j+=1
		for i in range(nSplit):
			fps[i].close()
	return files


def filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile):
	return bowtie2Wrap.filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile)


# Appends all the appendFiles to outfile
def append_file(outfile, appendFiles):
	with open(outfile,'w') as out1:
		for file1 in appendFiles:
			if (file1 is not None):
				with open(file1,'r') as in2:
					for ln in in2:
						out1.write(ln)

# Appends all the sam appendFiles to outfile
def append_sam_file(outfile, appendFiles):
	with open(outfile,'w') as out1:
		# First, writing the header by merging headers from all files
		for file1 in appendFiles:
			if (file1 is not None):
				with open(file1,'r') as in2:
					for ln in in2:
						if ln[0] == '@':
							out1.write(ln)
		# Writing the body by merging body from all files
		for file1 in appendFiles:
			if (file1 is not None):
				with open(file1,'r') as in2:
					for ln in in2:
						if ln[0] != '@':
							out1.write(ln)

