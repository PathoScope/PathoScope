#!/usr/bin/python
# Initial Author: Solaiappan Manimaran
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

import os, re, csv, sys
from pathoscope.utils import samUtils
from pathoscope.utils import seqParse
from pathoscope.utils import pathoUtilsA
from pathoscope.pathodb import dbUtils
from pathoscope.pathoreport import xmlReport
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir) 

# ===========================================================
class PathoReportOptions:
	MIN_CONTIG_LEN = 101
	verbose = False
	contigFlag = False
	outDir = "."
	samtoolsHome = None
	samFile = None
	mysqlConf = None
	minContigLen = MIN_CONTIG_LEN
	def __init__(self, samFile):
		self.samFile = samFile


# Main entry function to PathoMap that does all the processing
def processPathoReport(pathoReportOptions):
	h_gisPerTi = {}
	h_annoT = {}
	h_ti_contig = {}
	if pathoReportOptions.samFile is not None:
		(_, tail) = os.path.split(pathoReportOptions.samFile)
		(base, _) = os.path.splitext(tail)
		tsvFile = base+'.tsv'
		outTsv = pathoReportOptions.outDir + os.sep + tsvFile

		xmlFile = base+'.xml'
		outXml = pathoReportOptions.outDir + os.sep + xmlFile
		
		run_param_str =''
		(h_refRead, h_refScore, reads, h_readSequence, h_gisPerTi, h_tiRef, U, NU, genomes, pi, initPi) = \
			samUtils.findAlignmentReadPercentage(pathoReportOptions.samFile)
		
		(bestHitFinalReads, bestHitFinal, level1Final, level2Final) = \
			computeBestHit(U, NU, genomes, reads)
		
		for j in NU:
			NU[j][2] = NU[j][1]
		(bestHitInitialReads, bestHitInitial, level1Initial, level2Initial) = \
			computeBestHit(U, NU, genomes, reads)
		header = ['Genome', 'MapQ Guess', 'MapQ Best Hit', 'MapQ Best Hit Read Numbers', \
			'MapQ High Confidence Hits', 'MapQ Low Confidence Hits', 'Alignment Guess', \
			'Alignment Best Hit', 'Alignment Best Hit Read Numbers', \
			'Alignment High Confidence Hits', 'Alignment Low Confidence Hits']
		nR = len(reads)
		nG = len(h_refRead)
		(_, _, _, _, _, _, _, _, _, _, _) = write_tsv_report(outTsv, nR, nG, 
			pi, genomes, initPi, bestHitInitial, bestHitInitialReads, bestHitFinal, 
			bestHitFinalReads, level1Initial, level2Initial, level1Final, level2Final, header)
		
		if pathoReportOptions.contigFlag:
			bamFile = samUtils.sam2bam(pathoReportOptions.samFile, pathoReportOptions.samtoolsHome)
			#3)for each genome, get covered fragment and also make sure that each fragment has a high quality
			refConsFq = samUtils.samtools_consensus(bamFile, pathoReportOptions.samtoolsHome)
			#3.1) get delimiter and send a query to mysql to retrieve annotation
			(h_annoT, h_ti_contig) =get_genome_annotation_in_mysql(\
				refConsFq, pathoReportOptions.minContigLen, 
				pathoReportOptions.mysqlConf, h_annoT, h_ti_contig)
			xmlReport.writePathoXML(run_param_str, outXml, h_annoT, h_ti_contig, 
				h_refRead, h_refScore, h_gisPerTi, h_tiRef, reads, h_readSequence, 
				pathoReportOptions.samFile, pathoReportOptions.mysqlConf)
		else:
			h_annoT = simple_genome_annotation(h_gisPerTi, pathoReportOptions.mysqlConf, h_annoT)
			xmlReport.writePathoXML(run_param_str, outXml, h_annoT, h_ti_contig, 
				h_refRead, h_refScore, h_gisPerTi, h_tiRef, reads, h_readSequence, 
				pathoReportOptions.samFile, pathoReportOptions.mysqlConf)
	

#===============================================
'''
objective: this is a core def in retrieve_genome_annotation_from_sam. Having template(consensus fastq) ref covered by pileup reads, we want to retrieve gene annotations from mysql
'''
def get_genome_annotation_in_mysql(\
	refConsFq, minContigLen, MySqlConf, h_annoT, h_ti_contig):
	
	START,END = range(2)
	SUBGI,GENE,LOCS_TAG,PROID,STBP,EDBP = range(6)
	NAs = 'X'
	useMysql=True
	con = None
	#(hostname,port,user,passwd,defaultDb)=range(5)
	(_,_,_,passwd,_)=range(5)
	if MySqlConf[passwd]==NAs: #then, we do not use mysql
		useMysql=False
	if useMysql:
		con = dbUtils.init_mysql_innocentive(MySqlConf,0)

	fp = open(refConsFq,'r')
	#debugCnt = 0 #debug
	for r in seqParse.parse(fp,'fastq'): # for each covered genome

		covRange = selectConsensusContigs(r,minContigLen,-1) #disable checking seq complexity of contig
		
		if not covRange:
			continue
		C = len(covRange)

		#extract ti and gi
		refName = r.id
		mObj=re.search(r'ti\|(\d+)\|org\|([^|]+)\|gi\|(\d+)\|',r.id)
		if mObj:
			ti = mObj.group(1)
			gi = mObj.group(3)
		else:
			mObj=re.search(r'ti\|(\d+)\|gi\|(\d+)\|',r.id)
			if mObj and mObj.group(1)!="-1":
				ti = mObj.group(1)
				gi = mObj.group(2)
			else:
				mObj=re.search(r'gi\|(\d+)\|',r.id)
				if mObj:
					gi = mObj.group(1)

		if not h_ti_contig.get(ti,[]):
			h_ti_contig[ti]=[]
			
		for c in range(C):
			#contig = r[covRange[c][0]:covRange[c][1]+1]
			contigSeq = str(r.seq[covRange[c][0]:covRange[c][1]+1])
			#cqual = contig.letter_annotations["phred_quality"]
			#cLen = len(cqual)
			cLen = covRange[c][1]-covRange[c][0]+1
			#cqual_ave = 1.*sum(cqual)/cLen
			
			#h_ti_contig[ti].append([refName,cLen,str(contig.seq)])
			h_ti_contig[ti].append([refName,cLen,contigSeq])
		
		if con:
			mysql_sel_cmd = 'select sub_gi, gene, locus_tag, protein_id, stbp, edbp from giDelimT where gi = %s' % gi
			cur = con.cursor()
			cur.execute(mysql_sel_cmd)
			entr=cur.fetchall()
			if entr:
				#subgi2query=[]
				#subgiAnnot=[]
				#print r.id #debug
				#print covRange #debug
				for j in entr: #select which subgi sits within the covered genomic regions
					aStbp=int(j[STBP]); aEdbp=int(j[EDBP])

					A=aEdbp-aStbp+1
					notCoveredA=A
					minCoveredA2 = notCoveredA - 100
					
					reportA=False
					
					for i in range(C):
						#print '[subgi%s:%d - %d][cov:%d-%d]' % (gi,aStbp,aEdbp,covRange[START][i],covRange[END][i])
						notCoveredA -= pathoUtilsA.segments_intersect(aStbp,aEdbp,covRange[i][START],covRange[i][END])
						if notCoveredA<minCoveredA2:
							reportA=True
							break

					if reportA:
						selCmd = 'select ref_name, product from giAnnoT where gi = %s' % j[SUBGI]
						cur = con.cursor()
						cur.execute(selCmd)
						entr2 = cur.fetchone()
						ref_name=NAs; product=NAs
						if entr2:
							ref_name = entr2[0]; product = entr2[1]
						if h_annoT.get(ti,-1)==-1:
							h_annoT[ti]=[]
						h_annoT[ti].append([j[SUBGI],j[GENE],j[LOCS_TAG],j[PROID],ref_name,product])

	fp.close()
	if con:
		dbUtils.mysql_close(con)
	return h_annoT,h_ti_contig


#===============================================
'''
objective: this is a core def in retrieve_genome_annotation_from_sam. Having template(consensus fastq) ref covered by pileup reads, we want to retrieve gene annotations from mysql
'''
def simple_genome_annotation(h_gisPerTi, mySqlConf, h_annoT):
	
	#SUBGI,GENE,LOCS_TAG,PROID,STBP,EDBP = range(6)
	SUBGI,GENE,LOCS_TAG,PROID = range(4)
	NAs = 'X'
	useMysql=True
	con = None
	#(hostname,port,user,passwd,defaultDb)=range(5)
	(_,_,_,passwd,_)=range(5)
	if mySqlConf[passwd]==NAs: #then, we do not use mysql
		useMysql=False
	if useMysql:
		con = dbUtils.init_mysql_innocentive(mySqlConf,0)
	if con:
		for ti in h_gisPerTi:
			gis = h_gisPerTi[ti]
			for gi in gis:
				mysql_sel_cmd = 'select sub_gi, gene, locus_tag, protein_id, stbp, edbp from giDelimT where gi = %s' % gi
				cur = con.cursor()
				cur.execute(mysql_sel_cmd)
				entr=cur.fetchall()
				if entr:
					for j in entr: #select which subgi sits within the covered genomic regions
						selCmd = 'select ref_name, product from giAnnoT where gi = %s' % j[SUBGI]
						cur = con.cursor()
						cur.execute(selCmd)
						entr2 = cur.fetchone()
						ref_name=NAs; product=NAs
						if entr2:
							ref_name = entr2[0]; product = entr2[1]
						if h_annoT.get(ti,-1)==-1:
							h_annoT[ti]=[]
						h_annoT[ti].append([j[SUBGI],j[GENE],j[LOCS_TAG],j[PROID],ref_name,product])

	if con:
		dbUtils.mysql_close(con)
	return h_annoT

#===========================
def selectConsensusContigs(fqRec,minContigLen,kolCompxCutoff):
	
	covRanges2=[]
	covRanges=[]
	
	rSeq = str(fqRec.seq)
	#print rSeq #debug
	tmp=[match.start()+1 for match in re.finditer('^[^nN]|[nN][^nN]',rSeq)]
	if not tmp:
		return covRanges2
	if tmp[0]==1:
		tmp[0]=0
	covRanges.append(tmp)
	
	tmp=[match.start() for match in re.finditer('[^nN][nN]|[^nN]$',rSeq)]
	if not tmp:
		return covRanges2
	covRanges.append(tmp)
	
	C = len(covRanges[0])
	
	for c in range(C):
		if (covRanges[1][c]-covRanges[0][c]+1) <= minContigLen:
			continue

		if kolCompxCutoff>0:
			subSeq=rSeq[covRanges[0][c]:covRanges[1][c]+1]
			kx = pathoUtilsA.kolmogorov(subSeq)
			if kx > kolCompxCutoff:
				covRanges2.append([covRanges[0][c],covRanges[1][c]])
		else:
			covRanges2.append([covRanges[0][c],covRanges[1][c]])

	del covRanges
	return covRanges2
	
#===========================
'''
Computes the best hit read metrics 
'''
def computeBestHit(U, NU, genomes, read):
	bestHitReads=[0.0 for _ in genomes]
	level1Reads=[0.0 for _ in genomes]
	level2Reads=[0.0 for _ in genomes]
	for i in U: 
		bestHitReads[U[i][0]]+=1
		level1Reads[U[i][0]] += 1
	for j in NU:
		z = NU[j]
		ind = z[0]
		xnorm = z[2]
		bestGenome = max(xnorm)
		numBestGenome = 0
		for i in range(len(xnorm)):
			if (xnorm[i] == bestGenome):
				numBestGenome += 1
		if (numBestGenome == 0):
			numBestGenome = 1
		for i in range(len(xnorm)):
			if (xnorm[i] == bestGenome):
				bestHitReads[ind[i]] += 1.0/numBestGenome
				if (xnorm[i] >= 0.5):
					level1Reads[ind[i]] += 1
				elif (xnorm[i] >= 0.01):
					level2Reads[ind[i]] += 1
		
	nG = len(genomes)
	nR = len(read)
	bestHit = [bestHitReads[k]/nR for k in range(nG)]
	level1 = [level1Reads[k]/nR for k in range(nG)]
	level2 = [level2Reads[k]/nR for k in range(nG)]
	return bestHitReads, bestHit, level1, level2

# ===========================================================
# Function to create the tsv file report
def write_tsv_report(finalReport, nR, nG, pi, genomes, initPi, bestHitInitial, bestHitInitialReads, 
		bestHitFinal, bestHitFinalReads, level1Initial, level2Initial, level1Final, level2Final,
		header):
	with open(finalReport, 'wb') as oFp:
		tmp = zip(pi,genomes, initPi, bestHitInitial, bestHitInitialReads, bestHitFinal, 
			bestHitFinalReads, level1Initial, level2Initial, level1Final, level2Final)
		tmp = sorted(tmp,reverse=True) # Sorting based on Final Guess
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = zip(*tmp)
		for i in range(len(x10)):
			if (x1[i] < 0.01 and x10[i] <= 0 and x11[i] <= 0):
				break
			if i == (len(x10)-1):
				i += 1
		tmp = zip (x2[:i], x1[:i], x6[:i], x7[:i], x10[:i], x11[:i], x3[:i], x4[:i], x5[:i], x8[:i], x9[:i]) # Changing the column order here
		csv_writer = csv.writer(oFp, delimiter='\t')
		header1 = ['Total Number of Aligned Reads:', nR, 'Total Number of Mapped Genomes:', nG]
		csv_writer.writerow(header1)
		csv_writer.writerow(header)
		csv_writer.writerows(tmp)
	return (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)

