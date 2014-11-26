#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Functions to read sam alignment file and obtain some information

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

import re, math, os


def findAlignmentReadPercentage_lite(samfile,taxonomyLevel1):
	sfp = open(samfile,'r') #open the original sam file
	hu_ref_read = {}
	h_readId = {}
	readCnt = 0
	h_refRead = {}
	for sln in sfp:
		if (sln[0] == '@' or sln[0] == '#'):
			continue
		s = sln.split('\t')
		if s[2] == '*':
			continue

		readId = s[0]
		if taxonomyLevel1==1:
			mObj=re.search(r'ti:(\S+)\|gi',s[2])
			refId = mObj.group(1)
		else:
			refId = s[2]
		
		#update h_refRead
		hkey = refId+':'+s[0]
		if hu_ref_read.get(hkey,-1) == -1:
			hu_ref_read[hkey]=1
			if h_refRead.get(refId,['X'])[0] == 'X':
				h_refRead[refId] = []
			h_refRead[refId].append(s[0])

		if h_readId.get(s[0],-1) == -1:
			h_readId[s[0]] = 1
			readCnt+=1
			
	sfp.close()
	return h_refRead,readCnt
	
# ===========================================================
def findAlignmentReadPercentage(samFile):
	U = {}
	NU = {}
	h_readId = {}
	h_readSequence = {}
	h_refId = {}
	h_refRead = {}
	h_refScore = {}
	h_initRefScore = {}
	h_gisPerTi = {}
	h_tiRef = {}
	genomes = []
	reads =[]
	gCnt = 0
	rCnt = 0
	maxScore = None
	minScore = None

	with open(samFile,'r') as in1:
		for ln in in1:
			if (ln[0] == '@' or ln[0] == '#'):
				continue
			l = ln.split('\t')
			readId=l[0]
			refId=l[2]
			mObj=re.search(r'ti\|(\d+)\|org\|([^|]+)\|gi\|(\d+)\|',refId)
			if mObj:
				refId = "ti|"+mObj.group(1)+"|org|"+mObj.group(2)
				ti = mObj.group(1)
				gi = mObj.group(3)
				gis = h_gisPerTi.get(ti,[])
				gis.append(gi)
				h_gisPerTi[ti] = gis
				refIdName = h_tiRef.get(ti,[])
				if len(refIdName)==0:
					refIdName.append(refId)
					name = mObj.group(2)
					refIdName.append(name)
					h_tiRef[ti] = refIdName
			else:
				mObj=re.search(r'ti\|(\d+)\|gi\|(\d+)\|',refId)
				if mObj and mObj.group(1)!="-1":
					refId = "ti|"+mObj.group(1)
					ti = mObj.group(1)
					gi = mObj.group(2)
					gis = h_gisPerTi.get(ti,[])
					gis.append(gi)
					h_gisPerTi[ti] = gis
					refIdName = h_tiRef.get(ti,[])
					if len(refIdName)==0:
						refIdName.append(refId)
						name = refId
						refIdName.append(name)
						h_tiRef[ti] = refIdName
			if refId == '*':
				continue
			
			aScore = findSamAlignScore(l)
			if aScore is None:
				continue
			mapq = float(l[4])
			mapq2 = mapq/(-10.0)
			pScore = 1.0 - pow(10,mapq2)
			#pScore = int(round(pScore*100)) # Converting to integer to conserve memory space
			#if pScore < 1:
			if pScore < 0.01:
				continue
			if ((maxScore == None) or (aScore > maxScore)):
				maxScore = aScore
			if ((minScore == None) or (aScore < minScore)):
				minScore = aScore
			
			gIdx = h_refId.get(refId,-1)
			if gIdx == -1:
				gIdx = gCnt
				h_refId[refId] = gIdx
				h_refRead[refId] = []
				h_refScore[refId] = 0
				h_initRefScore[refId] = 0
				genomes.append(refId)
				gCnt += 1

			if not readId in h_refRead.get(refId,[-1]): #no redun reads allowed!
				h_refRead[refId].append(readId)

			rIdx = h_readId.get(readId,-1)
			if rIdx == -1:
				#hold on this new read
				#first, wrap previous read profile and see if any previous read has a same profile with that!
				rIdx = rCnt
				h_readId[readId] = rIdx
				reads.append(readId)
				h_readSequence[readId] = l[9]
				rCnt += 1
				U[rIdx] = [[gIdx], [aScore], [pScore], aScore]
			else:
				if (rIdx in U):
					if gIdx in U[rIdx][0]:
						continue
					NU[rIdx] = U[rIdx]
					del U[rIdx]
				if gIdx in NU[rIdx][0]:
					continue
				NU[rIdx][0].append(gIdx)
				NU[rIdx][1].append(aScore)
				NU[rIdx][2].append(pScore)
				if aScore > NU[rIdx][3]:
					NU[rIdx][3] = aScore
	#			length = len(NU[rIdx][1])
	#			NU[rIdx][2] = [1.0/length]*length

	del h_refId, h_readId
	(U, NU) = rescale_samscore(U, NU, maxScore, minScore)
	Uweights = 0.0
	for rIdx in U:
		gIdx = U[rIdx][0][0]
		U[rIdx] = [gIdx, U[rIdx][1][0]] #keep gIdx and alignment score for weights
		Uweights += U[rIdx][1]
		h_refScore[genomes[gIdx]] += 1.0*U[rIdx][1]
		h_initRefScore[genomes[gIdx]] += 1.0*U[rIdx][1]
	NUweights = 0.0
	for rIdx in NU:
		pScoreSum = sum(NU[rIdx][2])
		NU[rIdx][2] = [1.0*k/pScoreSum for k in NU[rIdx][2]] #Normalizing pScore
		initpScoreSum = sum(NU[rIdx][1])
		NU[rIdx][1] = [1.0*k/initpScoreSum for k in NU[rIdx][1]] #Normalizing pScore
		NUweights += NU[rIdx][3]
		for j in range(len(NU[rIdx][0])):
			gIdx = NU[rIdx][0][j]
			h_refScore[genomes[gIdx]] += 1.0*NU[rIdx][2][j]*NU[rIdx][3]
			h_initRefScore[genomes[gIdx]] += 1.0*NU[rIdx][1][j]*NU[rIdx][3]
	weight = Uweights + NUweights
	for refId in h_refScore:
		h_refScore[refId] /= weight
		h_initRefScore[refId] /= weight
	pi = [1.0*h_refScore[refId] for refId in genomes] #Normalizing pScore
	initPi = [1.0*h_initRefScore[refId] for refId in genomes] #Normalizing pScore
	
	return h_refRead, h_refScore, reads, h_readSequence, h_gisPerTi, h_tiRef, U, NU, genomes, pi, initPi

#===========================
'''
Objective: to convert sam to bam file. 
'''
def sam2bam(samFile, samtoolsHome):
	samtoolsBin = 'samtools'
	(head, tail) = os.path.split(samFile)
	(base, _) = os.path.splitext(tail)
	if len(head)>0:
		bamfile = head + os.sep + base+'.bam'
		sortBamPrefix = head + os.sep + 'sorted_'+base
	else:
		bamfile = base+'.bam'
		sortBamPrefix = 'sorted_'+base
	if samtoolsHome is not None:
		samtoolsPath = samtoolsHome + os.sep + samtoolsBin
	else:
		samtoolsPath = samtoolsBin
	cmd='%s view -bS %s > %s' % (samtoolsPath, samFile, bamfile)
	print cmd
	os.system(cmd)
	#samtools sorting .bam file
	cmd='%s sort %s %s' % (samtoolsPath, bamfile, sortBamPrefix)
	print cmd
	os.system(cmd)
	sortBamFile = sortBamPrefix+'.bam'
	#samtools index .bam file
	cmd='%s index %s' % (samtoolsPath, sortBamFile)
	print cmd
	os.system(cmd)
	return (sortBamFile)

#===========================================================
def samtools_consensus(bamFile, samtoolsHome):
	samtoolsBin = 'samtools'
	bcftoolsBin = 'bcftools'
	vcfutilsBin = 'vcfutils.pl'
	if samtoolsHome is not None:
		samtoolsPath = samtoolsHome + os.sep + samtoolsBin
		bcftoolsPath = samtoolsHome + os.sep + bcftoolsBin
		vcfutilsPath = samtoolsHome + os.sep + vcfutilsBin
	else:
		samtoolsPath = samtoolsBin
		bcftoolsPath = bcftoolsBin
		vcfutilsPath = vcfutilsBin
	(head, tail) = os.path.split(bamFile)
	(base, _) = os.path.splitext(tail)
	if len(head)>0:
		refConsFq = head + os.sep + base+'_cns.fq'
	else:
		refConsFq = base+'_cns.fq'
	#[example]samtools mpileup -uD -m 3 -F 0.0002 -f NC_001803.fa 287A_cDNA_33_at_noBacHg-topReadId-em_bwasw_hit.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq
	cmd='%s mpileup -A -u -D -Q 1 %s | %s view -cg - | %s vcf2fq > %s' % (samtoolsPath, bamFile, bcftoolsPath, vcfutilsPath, refConsFq)
	print cmd
	os.system(cmd)
	return refConsFq

# ===========================================================
# Finds the alignment score + read length for the given entry from sam file
def findSamAlignScore(l):
	useMapq = True
	readL = 1.*len(l[9])
	for i in range(11, len(l)):
		if useMapq and l[i].startswith('AS:i:'):
			aScore=int(l[i][5:])
			useMapq = False
		elif l[i].startswith('YS:i:'):
			aScore+=int(l[i][5:])
			readL = 2*readL # For paired end we simply multiply read length by 2
			break
	if useMapq:
		aScore = None
	else:
		aScore = aScore+readL
	return aScore

# ===========================================================
# Finds the alignment score + read length for the given entry from sam file
# Returns the higher of the two alignment scores for paired end entries
def findSamAlignHiScore(l):
	useMapq = True
	readL = 1.*len(l[9])
	for i in range(11, len(l)):
		if useMapq and l[i].startswith('AS:i:'):
			aScore=int(l[i][5:])
			useMapq = False
		elif l[i].startswith('YS:i:'):
			aScore1=int(l[i][5:])
			if (aScore1 > aScore):
				aScore = aScore1
			break
	if useMapq:
		aScore = None
	else:
		aScore = aScore+readL
	return aScore

# ===========================================================
# Rescaling the sam alignment score and taking exponent
def rescale_samscore(U, NU, maxScore, minScore):
	if (minScore < 0):
		scalingFactor = 100.0 / (maxScore - minScore)
	else:
		scalingFactor = 100.0 / (maxScore)
	for rIdx in U:
		if (minScore < 0):
			U[rIdx][1][0] = U[rIdx][1][0] - minScore
		U[rIdx][1][0] = math.exp(U[rIdx][1][0] * scalingFactor)
		U[rIdx][3] = U[rIdx][1][0]
	for rIdx in NU:
		NU[rIdx][3] = 0.0
		for i in range(0, len(NU[rIdx][1])):
			if (minScore < 0):
				NU[rIdx][1][i] = NU[rIdx][1][i] - minScore
			NU[rIdx][1][i] = math.exp(NU[rIdx][1][i] * scalingFactor)
			if NU[rIdx][1][i] > NU[rIdx][3]:
				NU[rIdx][3] = NU[rIdx][1][i]
	return (U, NU)
