#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Some miscellaneous utility functions used by pathoscope

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

import os
import zlib
import subprocess as sp

# ===========================================================
def file_len(fname):
	with open(fname) as f:
		for i, _ in enumerate(f):
			pass
	return (i + 1)
	
# ===========================================================
def ensure_dir(d):
	if not os.path.exists(d):
		os.makedirs(d)

#========================================
'''
objective: compute overlapped seg len
min1                 max1
|-------------------|
      |----------------|
     min2               max2
'''
def segments_intersect(min1, max1, min2, max2):
	return max(0, min(max1, max2) - max(min1, min2))

def kolmogorov(s):
	l = float(len(s))
	if l==0:
		return -1.0
	compr = zlib.compress(s)
	c = float(len(compr))
	kCompx = c/l
	#print kCompx #debug
	return kCompx

def separateDirFn(fullPathFile,delimiter):
	loc=fullPathFile.rfind(delimiter)
	if loc:
		file1=fullPathFile[:loc]
		file2=fullPathFile[loc+1:]
	else:
		file1=fullPathFile
		file2=''
	return (file1,file2)
	
def ncbi_eutil(gi,genbank_id,entry2sel):
	
	NOT_VALID=-1
	TAXONOMY_ID=1
	
	if entry2sel==TAXONOMY_ID:
		field2sel='TSeq_taxid'
	
	#print 'sending query to ncbi eutil %s(gi:%d)' % (genbank_id,gi) #debug
	
	queries=[gi,genbank_id]
	taxonId=NOT_VALID
	
	for q in queries:
		curlCmd='curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=xml\" | grep %s | cut -d \'>\' -f 2 | cut -d \'<\' -f 1 | tr -d \"\\n\"' % (q,field2sel)
		
		#print curlCmd
		proc = sp.Popen(curlCmd, stdout=sp.PIPE, shell=True)
		(msg, _) = proc.communicate() # err is the second object returned
		msg=msg.strip()

		if msg:
			taxonId=int(msg)
			break
	
	return taxonId

#=================================
def ex_wget_download(protozoaD):
	
	protozoaDfn='%s/ncbi_protozoa.fa' % protozoaD
	
	if not os.path.exists(protozoaDfn):
		cmd = 'cd %s; wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/protozoa.*.rna.fna.gz\n' % protozoaD;
		cmd = '%swget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/protozoa.*.genomic.fna.gz\n' % cmd
		cmd = '%sgunzip %s/*.gz\n' % (cmd,protozoaD)
		cmd = '%scat %s/protozoa*.fna > %s/ncbi_protozoa.fa\n' % (cmd,protozoaD,protozoaD)
		cmd = '%srm -rf %s/protozoa*.fna\n' % (cmd,protozoaD)
		
		print 'downloading data...'
		print cmd
		os.system(cmd)
		print 'done'
	
	return protozoaDfn 

#==================================
def wget_download2(url,downFbase,downExt,downloadD,operation,outName):
	
	if not downloadD:
		os.makedirs(downloadD)

	outDfn = downloadD + '/' + outName
	if os.path.exists(outDfn):
		return outDfn

	downFn='%s.%s' % (downFbase,downExt)
	downDfn='%s/%s' % (downloadD,downFn)
	if not os.path.exists(downDfn):
		cmd='cd %s' % downloadD
		cmd='%s;wget %s' % (cmd,url)
		os.system(cmd)
	
	cmd = 'cd %s;' % downloadD
	if downExt=='gz':
		extOp='gunzip'
		
	elif downExt=='tar.gz':
		extOp='tar zxvf'
	
	if operation=='select':
		cmd='%s %s %s' % (cmd, extOp, downFn)
		print cmd
		os.system(cmd)
		outDfn='%s/%s' % (downloadD,outName)
	else:
		cmd = 'ls -la %s' % downloadD
		print cmd
		os.system(cmd)
		outDfn='X'
		
	return outDfn

def file_tag(filename,tag,newFileExt):
	(filebase,fileExt)=separateDirFn(filename,'.')
	if newFileExt=='X':
		newFileExt=fileExt
		
	if tag == 'X':
		taggedFn=filebase+'.'+newFileExt
	else:
		taggedFn=filebase+'_'+tag+'.'+newFileExt
	return taggedFn
	
# ===========================================================
def search_cat_in_online_taxonomy(ti,ncbiLineage,tmpFn):
	#catIdx=['N','A','B','E','EP','EH','EF','V','U','O']
	field2sel='Lineage'
	curlCmd='curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%d&rettype=fasta&retmode=xml\" | grep -iP \'%s\' > %s' % (ti,field2sel,tmpFn)
	os.system(curlCmd)
	
	NOT_VALID='N'
	hit_cat=NOT_VALID
	
	for c in ncbiLineage:
		cmd='grep -iP \'%s\' %s' % (c,tmpFn)
		proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
		(msg, _) = proc.communicate() # err is the second object returned
		msg=msg.strip()
		if msg:
			hit_cat=c
			break
	return hit_cat

