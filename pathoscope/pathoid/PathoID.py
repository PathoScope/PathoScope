#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Functions to read alignment file (sam/gnu-sam or bl8), run EM algorithm
# and output report file that can be opened in Excel and also 
# output updated alignment file (sam/gnu-sam or bl8)  

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

import os, math, csv, re, sys
from pathoscope.utils import pathoUtilsA
from pathoscope.utils import samUtils
from pathoscope.pathoreport import PathoReportA

class PathoIdOptions:
	ali_file = ""
	verbose = False
	score_cutoff = 0.01
	exp_tag = ""
	ali_format = "sam"
	outdir = ""
	emEpsilon = 0.01 
	maxIter = 50
	piPrior = 0
	thetaPrior = 0
	out_matrix_flag = True
	noalign = False
	def __init__(self, ali_file):
		self.ali_file = ali_file

# ===========================================================
def conv_align2GRmat(aliDfile,pScoreCutoff,aliFormat):
	in1 = open(aliDfile,'r')
	U = {}
	NU = {}
	h_readId = {}
	h_refId = {}
	genomes = []
	read =[]
	gCnt = 0
	rCnt = 0

	maxScore = None
	minScore = None
	for ln in in1:
		if (ln[0] == '@' or ln[0] == '#'):
			continue

		l = ln.split('\t')
		
		readId=l[0]
		if (aliFormat == 0 or aliFormat == 1): # gnu-sam or sam
			if int(l[1])&0x4 == 4: # bitwise FLAG - 0x4 : segment unmapped
				continue
			refId=l[2]
		elif (aliFormat == 2): # bl8
			refId=l[1]
		
		if refId == '*':
			continue
		
		#refId=refId.split("ti:")[-1]
		mObj=re.search(r'ti\|(\d+)\|org\|([^|]+)\|gi',refId)
		if mObj:
			refId = "ti|"+mObj.group(1)+"|org|"+mObj.group(2)
		else:
			mObj=re.search(r'ti\|(\d+)\|gi',refId)
			if mObj and mObj.group(1)!="-1":
				refId = "ti|"+mObj.group(1)

		(pScore, skipFlag) = find_entry_score(ln, l, aliFormat, pScoreCutoff)
		if skipFlag:
			continue
		if ((maxScore == None) or (pScore > maxScore)):
			maxScore = pScore
		if ((minScore == None) or (pScore < minScore)):
			minScore = pScore
		
		gIdx = h_refId.get(refId,-1)
		if gIdx == -1:
			gIdx = gCnt
			h_refId[refId] = gIdx
			genomes.append(refId)
			gCnt += 1

		rIdx = h_readId.get(readId,-1)
		if rIdx == -1:
			#hold on this new read
			#first, wrap previous read profile and see if any previous read has a same profile with that!
			rIdx = rCnt
			h_readId[readId] = rIdx
			read.append(readId)
			rCnt += 1
			U[rIdx] = [[gIdx], [pScore], [float(pScore)], pScore]
		else:
			if (rIdx in U):
				if gIdx in U[rIdx][0]:
					continue
				NU[rIdx] = U[rIdx]
				del U[rIdx]
			if gIdx in NU[rIdx][0]:
				continue
			NU[rIdx][0].append(gIdx)
			NU[rIdx][1].append(pScore)
			if pScore > NU[rIdx][3]:
				NU[rIdx][3] = pScore
#			length = len(NU[rIdx][1])
#			NU[rIdx][2] = [1.0/length]*length

	in1.close()

	if (aliFormat == 1): # sam
		(U, NU) = samUtils.rescale_samscore(U, NU, maxScore, minScore)

	del h_refId, h_readId
	for rIdx in U:
		U[rIdx] = [U[rIdx][0][0], U[rIdx][1][0]] #keep gIdx and score only
	for rIdx in NU:
		pScoreSum = sum(NU[rIdx][1])
		NU[rIdx][2] = [k/pScoreSum for k in NU[rIdx][1]] #Normalizing pScore

	return U, NU, genomes, read

# ===========================================================
# Entry function to PathoID
# Does the reassignment and generates a tsv file report
def pathoscope_reassign(pathoIdOptions):
	out_matrix = pathoIdOptions.out_matrix_flag
	verbose = pathoIdOptions.verbose
	scoreCutoff = pathoIdOptions.score_cutoff
	expTag = pathoIdOptions.exp_tag
	ali_format = pathoIdOptions.ali_format
	ali_file = pathoIdOptions.ali_file
	outdir = pathoIdOptions.outdir
	emEpsilon = pathoIdOptions.emEpsilon
	maxIter = pathoIdOptions.maxIter
	upalign = not(pathoIdOptions.noalign)
	piPrior = pathoIdOptions.piPrior
	thetaPrior = pathoIdOptions.thetaPrior
	
	if float(os.stat(ali_file).st_size)<1.0:
		print 'the alignment file [%s] is empty.' % ali_file
		sys.exit(1)

	if ali_format == 'gnu-sam':
		aliFormat = 0
		if verbose:
			print "parsing gnu-sam file/likelihood score/reads and mapped genomes..."
	elif ali_format == 'sam': #standard sam
		aliFormat = 1
		if verbose:
			print "parsing sam file/likelihood score/reads and mapped genomes..."
	elif ali_format == 'bl8': #blat m8 format
		aliFormat = 2
		if verbose:
			print "parsing bl8 file/likelihood score/reads and mapped genomes..."
	else:
		print "unknown alignment format file..."
		return
	(U, NU, genomes, reads) = conv_align2GRmat(ali_file,scoreCutoff,aliFormat)
	
	nG = len(genomes)
	nR = len(reads)
	if verbose:
		print "EM iteration..."
		print "(Genomes,Reads)=%dx%d" % (nG, nR)
		print "Delta Change:"
	
	if out_matrix:
		if verbose:
			print "writing initial alignment ..."
		out_initial_align_matrix(genomes, reads, U, NU, expTag, ali_file, outdir)	

	(bestHitInitialReads, bestHitInitial, level1Initial, level2Initial) = \
		PathoReportA.computeBestHit(U, NU, genomes, reads)
	
	(initPi, pi, _, NU) = pathoscope_em(U, NU, genomes, maxIter, emEpsilon, verbose,
		piPrior, thetaPrior)
	tmp = zip(initPi,genomes)
	tmp = sorted(tmp,reverse=True) #similar to sort row
	
	if out_matrix:
		initialGuess = outdir + os.sep + expTag + '-initGuess.txt'
		oFp = open(initialGuess,'wb')
		csv_writer = csv.writer(oFp, delimiter='\t')
		csv_writer.writerows(tmp)
		oFp.close()
	
	del tmp
	
	(bestHitFinalReads, bestHitFinal, level1Final, level2Final) = \
		PathoReportA.computeBestHit(U, NU, genomes, reads)

	if out_matrix:
		finalGuess = outdir + os.sep + expTag + '-finGuess.txt'
		oFp = open(finalGuess,'wb')
		tmp = zip(pi,genomes)
		tmp = sorted(tmp,reverse=True)
		csv_writer = csv.writer(oFp, delimiter='\t')
		csv_writer.writerows(tmp)
		oFp.close()

	finalReport = outdir + os.sep + expTag +'-'+ ali_format + '-report.tsv'
	header = ['Genome', 'Final Guess', 'Final Best Hit', 'Final Best Hit Read Numbers', \
		'Final High Confidence Hits', 'Final Low Confidence Hits', 'Initial Guess', \
		'Initial Best Hit', 'Initial Best Hit Read Numbers', \
		'Initial High Confidence Hits', 'Initial Low Confidence Hits']
	(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) = PathoReportA.write_tsv_report(
		finalReport, nR, nG, pi, genomes, initPi, bestHitInitial, bestHitInitialReads, 
		bestHitFinal, bestHitFinalReads, level1Initial, level2Initial, level1Final, 
		level2Final, header)
	
	reAlignfile = ali_file
	if upalign:
		reAlignfile = rewrite_align(U, NU, ali_file, scoreCutoff, aliFormat, outdir)

	return (finalReport, x2, x3, x4, x5, x1, x6, x7, x8, x9, x10, x11, reAlignfile)


# ===========================================================
# This is the main EM algorithm
# ===========================================================
def pathoscope_em(U, NU, genomes, maxIter, emEpsilon, verbose, piPrior, thetaPrior):
	G = len(genomes)

	### Initial values
	pi = [1./G for _ in genomes]
	initPi = pi
	theta = [1./G for _ in genomes]
	
	pisum0=[0 for _ in genomes]
	Uweights = [U[i][1] for i in U] # weights for unique reads...
	maxUweights = 0
	Utotal = 0
	if Uweights:
		maxUweights = max(Uweights)
		Utotal = sum(Uweights)
	for i in U: 
		pisum0[U[i][0]]+=U[i][1]
	
	#pisum0/Utotal would be the weighted proportions of unique reads assigned to each genome (weights are alignment scores)_. 
	
	#need to change the structure for U matrix
	#notes NU weights are unnormalized
	# pull out a weighted likelihood score
	
	### data structure
	### 3 unique reads: 2 reads to genome1 and 1 read to genome4: readnum:[genome,reascore]
	# U = {0: 0, 1: 0, 2: 3}
	# U = {0: [0,1], 1: [0,.5], 2: [3,1]}
	### non-unique reads: 3 total reads  readnum:[[genomes],[qij],[xij]]
	#NU = {0: [[0, 2, 3], [0.4, 0.2, 0.4] , [0.33, 0.33, 0.33],.4, 1: [[0, 1], [0.6, 0.4] , [0.5,0.5]], 2: [[1, 3], [0.5, 0.5] , [0.5,0.5]]}
	### Genome hash
	#genomes = {0:"ecoli", 1:"strep", 2:"anthrax", 3:"plague"}
	#NUweights = [max(NU[i][1]) for i in NU] # weights for non-unique reads...
	NUweights = [NU[i][3] for i in NU] # weights for non-unique reads...
	maxNUweights = 0
	NUtotal = 0
	if NUweights:
		maxNUweights = max(NUweights)
		NUtotal = sum(NUweights)
	priorWeight = max(maxUweights, maxNUweights)
	lenNU=len(NU)
	if lenNU==0:
		lenNU=1

	for i in range(maxIter):  ## EM iterations--change to convergence 
		pi_old = pi
		thetasum=[0 for k in genomes]
		
		# E Step 

		for j in NU: #for each non-uniq read, j
			z = NU[j] 
			ind = z[0] #a set of any genome mapping with j
			pitmp = [pi[k] for k in ind]			### get relevant pis for the read
			thetatmp = [theta[k] for k in ind]  	### get relevant thetas for the read
			xtmp = [1.*pitmp[k]*thetatmp[k]*z[1][k] for k in range(len(ind))]  ### Calculate unormalized xs
			xsum = sum(xtmp)
			if xsum == 0:
				xnorm = [0.0 for k in xtmp]  		### Avoiding dividing by 0 at all times
			else:
				xnorm = [1.*k/xsum for k in xtmp]  	### Normalize new xs
			
			NU[j][2] = xnorm     					## Update x in NU 
			
			for k in range(len(ind)):
				#thetasum[ind[k]] += xnorm[k]   		### Keep running tally for theta
				thetasum[ind[k]] += xnorm[k]*NU[j][3]   ### Keep weighted running tally for theta
					
		# M step	
		pisum = [thetasum[k]+pisum0[k] for k in range(len(thetasum))]   ### calculate tally for pi
		pip = piPrior*priorWeight # pi prior - may be updated later
		#pi = [(1.*k+pip)/(len(U)+len(NU)+pip*len(pisum)) for k in pisum]  		## update pi
		#pi = [1.*k/G for k in pisum]  		## update pi
		totaldiv = Utotal+NUtotal
		if totaldiv==0:
			totaldiv=1
		pi = [(1.*k+pip)/(Utotal+NUtotal+pip*len(pisum)) for k in pisum]  		## update pi
		if (i == 0):
			initPi = pi
		
		thetap = thetaPrior*priorWeight # theta prior - may be updated later
		NUtotaldiv = NUtotal
		if NUtotaldiv==0:
			NUtotaldiv=1
		theta = [(1.*k+thetap)/(NUtotaldiv+thetap*len(thetasum)) for k in thetasum]
		#theta = [(1.*k+thetap)/(lenNU+thetap*len(thetasum)) for k in thetasum]

		cutoff = 0.0
		for k in range(len(pi)):
			cutoff += abs(pi_old[k]-pi[k])
		if verbose:
			print "[%d]%g" % (i,cutoff)
		if (cutoff <= emEpsilon or lenNU==1):
			break

	return initPi, pi, theta, NU

def out_initial_align_matrix(ref, read, U, NU, expTag, ali_file, outdir):
	genomeId = outdir + os.sep + expTag + '-genomeId.txt'
	oFp = open(genomeId,'wb')
	csv_writer = csv.writer(oFp, delimiter='\n')
	csv_writer.writerows([ref])
	oFp.close()

	readId = outdir + os.sep + expTag + '-readId.txt'
	oFp = open(readId,'wb')
	csv_writer = csv.writer(oFp, delimiter='\n')
	csv_writer.writerows([read])
	oFp.close()
	
# ===========================================================
# Generates the updated alignment file with the updated score
def rewrite_align(U, NU, aliDfile, pScoreCutoff, aliFormat, outdir):
	pathoUtilsA.ensure_dir(outdir)
	f = os.path.basename(aliDfile)
	reAlignfile = outdir + os.sep + 'updated_' + f
	
	with open(reAlignfile,'w') as of:
		with open(aliDfile,'r') as in1:
			h_readId = {}
			h_refId = {}
			genomes = []
			read =[]
			gCnt = 0
			rCnt = 0
		
			mxBitSc = 700
			sigma2 = 3
			for ln in in1:
				if (ln[0] == '@' or ln[0] == '#'):
					of.write(ln)
					continue
		
				l = ln.split('\t')
				
				readId=l[0]
				if (aliFormat == 0 or aliFormat == 1): # gnu-sam or sam
					#refId=l[2].split("ti:")[-1]
					refId=l[2]
					if int(l[1])&0x4 == 4: # bitwise FLAG - 0x4 : segment unmapped
						continue
				elif (aliFormat == 2): # bl8
					refId=l[1]
				
				if refId == '*':
					continue

				mObj=re.search(r'ti\|(\d+)\|org\|([^|]+)\|gi',refId)
				if mObj:
					refId = "ti|"+mObj.group(1)+"|org|"+mObj.group(2)
				else:
					mObj=re.search(r'ti\|(\d+)\|gi',refId)
					if mObj and mObj.group(1)!="-1":
						refId = "ti|"+mObj.group(1)
				
				(_, skipFlag) = find_entry_score(ln, l, aliFormat, pScoreCutoff)
				if skipFlag:
					continue
				
				gIdx = h_refId.get(refId,-1)
				if gIdx == -1:
					gIdx = gCnt
					h_refId[refId] = gIdx
					genomes.append(refId)
					gCnt += 1
		
				rIdx = h_readId.get(readId,-1)
				if rIdx == -1:
					#hold on this new read
					#first, wrap previous read profile and see if any previous read has a same profile with that!
					rIdx = rCnt
					h_readId[readId] = rIdx
					read.append(readId)
					rCnt += 1
					if rIdx in U:
						of.write(ln)
						continue
							
				if rIdx in NU:
					if (aliFormat == 0): # gnu-sam
						scoreComponents = l[12].split(':')
						(upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
						scoreComponents[2] = str(upPscore*pscoreSum)
						if (scoreComponents[2] < pScoreCutoff):
							continue
						l[12] = ':'.join(scoreComponents)
						ln = '\t'.join(l)
						of.write(ln)
					elif (aliFormat == 1): # sam
						(upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
						if (upPscore < pScoreCutoff):
							continue
						if (upPscore >= 1.0):
							upPscore = 0.999999
						mapq2 = math.log10(1 - upPscore)
						l[4] = str(int(round(-10.0*mapq2)))
						ln = '\t'.join(l)
						of.write(ln)
					elif (aliFormat == 2): # bl8
						(upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
						score = upPscore*pscoreSum
						if score <= 0.0:
							continue
						bitSc = math.log(score)
						if bitSc > mxBitSc:
							bitSc = mxBitSc
						l[10] = str(bitSc*sigma2)
						ln = '\t'.join(l)
						of.write(ln)

	return reAlignfile

# ===========================================================
# Function to find the updated score after pathoscope reassignment
def find_updated_score(NU, rIdx, gIdx):
	index = NU[rIdx][0].index(gIdx);
	pscoreSum = 0.0
	for pscore in NU[rIdx][1]:
		pscoreSum += pscore
	pscoreSum /= 100
	upPscore = NU[rIdx][2][index]
	return (upPscore, pscoreSum)

# ===========================================================
# Internal function to calculate the score from the alignment file entries
def find_entry_score(ln, l, aliFormat, pScoreCutoff):
	mxBitSc = 700
	sigma2 = 3
	skipFlag = False
	if (aliFormat == 0): # gnu-sam
		pScore = float(l[12].split(':')[2])
		if (pScore < pScoreCutoff):
			skipFlag = True
	elif (aliFormat == 1): # sam
		pScore = samUtils.findSamAlignScore(l)
		if pScore is None:
			mapq = float(l[4])
			mapq2 = mapq/(-10.0)
			pScore = 1.0 - pow(10,mapq2)
			if (pScore < pScoreCutoff):
				skipFlag = True
	elif (aliFormat == 2): # bl8
		eVal = float(l[10])
		if (eVal > pScoreCutoff):
			skipFlag = True
		bitSc = float(l[11])/sigma2
		if bitSc > mxBitSc:
			bitSc = mxBitSc
		pScore = math.exp(bitSc)
	#pScore = int(round(pScore*100)) # Converting to integer to conserve memory space
	#if pScore < 1:
	#	skipFlag = true
	return (pScore, skipFlag)

