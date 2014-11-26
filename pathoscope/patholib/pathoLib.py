#!/usr/bin/python
# input: ncbi nr or nt file, categories you want to select
#

import os, re, pickle
import subprocess as sp
from time import time
from pathoscope.utils import seqParse
from pathoscope.utils import pathoUtilsA
from pathoscope.pathodb import dbUtils

'''
objective: read app_build_nt_tgt() input argument and check consistency
input: 
-taxon_ids: something like 9606,4732,1234 in string format
-nodes_dmp: ncbi phylogeny tree database
output:
-taxon_ids2: something like [9606,4732,1234] in list format
Note:
'''
def parse_input_app_build_nt_tgt(taxon_ids):
	taxon_idsL=taxon_ids.split(',')
	NAs='X'
	GET_ALL_TAX=-2
	
	if taxon_idsL[0]!=NAs:
		taxon_ids2=[]
		for i in taxon_idsL:
			taxon_ids2.append(int(i))
	else:
		taxon_ids2=[GET_ALL_TAX]
	return taxon_ids2

#=============================
'''
objective: given fasta file having gi or some unique seq id, append taxonomy id in the front of header
keyword: gi, taxonomy, mapping
input:
-ncbiNt:
-taxon_ids:
-gi2taxDump:
-nodes_fn:
-enable_onlineF:

output:
-ncbiNt_ti
-ncbiNt_invalid

author: cjhong@bu.edu
date: 06/16/2013
note:
'''
def append_ti_into_fasta_app(ncbiNt, taxon_ids, exclude_taxon_ids, subTaxF, MySqlConf, 
		enable_descF, enable_onlineF, outPrefix, outDir):
	
	invalSelFlag = False # Not generating invalid fasta file by default
	GET_ALL_TAX=-2
	NAs='X'
	(hostname,port,user,passwd,defaultDb)=range(5)
	useHash=False
	
	workD = outDir + os.sep + 'ncbiDB'
	if not os.path.exists(workD):
		os.makedirs(workD)
	
	valSel = outDir + os.sep + outPrefix+'_ti.fa'
	invalSel = outDir + os.sep + outPrefix+'_ti_inval.fa'
	
	if MySqlConf[passwd]==NAs: #then, we use download gi2tax instead of using mysql query
		useHash=True
		gi2taxDump=getGi2TaxDump_online(workD)
		nodes_fn=getNodesDump_online(workD)
	
	if not useHash:
		con = dbUtils.init_mysql_innocentive(MySqlConf,0) #DON'T REMOVE LAST 0 or CHANGE TO 1
	
	if taxon_ids[0]==GET_ALL_TAX:
		subTaxons=[GET_ALL_TAX]
	else:
		#select all taxons under the given taxonomy lists
		subTaxons=taxon_ids
		if subTaxF:
			#build dic of ti nodes in phylogeny
			if useHash:
				h_tax_node2par = taxon_rel_in_hash(nodes_fn)
			else:
				h_tax_node2par = dbUtils.taxon_rel_in_hash_mysql(con, workD)
			subTaxons=get_allsub_taxons_phylo(h_tax_node2par,taxon_ids, exclude_taxon_ids)

	#append taxonomy id to fasta header
	if useHash:
		append_ti_into_fasta_hash(ncbiNt, gi2taxDump, subTaxons, enable_descF, enable_onlineF,
			valSel, invalSel, invalSelFlag)
	else:
		append_ti_into_fasta_mysql(con, ncbiNt, subTaxons, enable_descF, enable_onlineF,
			valSel, invalSel, invalSelFlag)
		dbUtils.mysql_close(con)
	
	return (valSel,invalSel)

#==============================
'''
objective: a wrapper to generate gi to taxonomy mapping list
keyword: gi, taxonomy, mapping
input:
-gi2taxDump:
output:
-mxGi: maximum number of gi in integer
-gi2Taxid: mapping list

author: cjhong@bu.edu
date: 06/16/2013
note:
'''
def gi2tax_list(gi2taxDump):
	
	#calculate max GI
	mxGi=read_gi_taxid_nucldmp_max(gi2taxDump)

	#read gi_taxid_nucl.dmp in list
	(_,gi2Taxid)=gi2taxid_offline(mxGi,gi2taxDump,True)
	
	return (mxGi,gi2Taxid)
	
#=================================
'''
objective: get a max gi number at the end of gi2taxDump where we assume that gi is in ascending order
input:
-gi2taxDump: ncbi gi2tax_nucl.dump file
output:
-mxGi: get a maximum gi
note: to be disappeared
'''
def read_gi_taxid_nucldmp_max(gi2taxDump):
	loadSuffix='maxGi'
	file2sv=gi2taxDump+'_'+loadSuffix
	print "find max gi..."
	if not os.path.isfile(file2sv):

		cmd='tail -n 1 %s | awk -F\"\\t\" \'{print $1}\'' % gi2taxDump
		proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
		(msg, err) = proc.communicate()
		mxGi=int(msg.strip())
			
		with open(file2sv,'w') as f:
			pickle.dump(mxGi,f)
	else:
		with open(file2sv) as f:
			mxGi=pickle.load(f)
	print "done!"
	return mxGi
	
#=================================
'''
objective: mapping gi to taxon in list using ncbi database (core)
keyword: gi, taxonomy, mapping
input:
-mxGi:
-gi2taxDump: three columns mapping (gi\tti1\tti2\n)
-rerunF:
output:
-file2sv: saved list file to reuse (optional)
-gi2Taxid: mapping list

author: cjhong@bu.edu
date: 06/16/2013
note:
'''
def gi2taxid_offline(mxGi,gi2taxDump,rerunF):
	loadSuffix='gi2tax'
	file2sv=gi2taxDump+'_'+loadSuffix
	print "converting gi to taxon_id..."
	print mxGi
	
	if rerunF or (not os.path.isfile(file2sv)):
		fp=open(gi2taxDump,'r')
		gi2Taxid=[]
		gi2Taxid=[0 for x in xrange(mxGi+1)]
		for line in fp:
			entry=line.split('\t')
			ti=int(entry[1])
			gi2Taxid[int(entry[0])]=ti
		fp.close()
		if not rerunF:
			with open(file2sv,'w') as f:
				pickle.dump(gi2Taxid,f)
	else:
		with open(file2sv,'r') as f:
			gi2Taxid=pickle.load(f)
	print "done!"
	return (file2sv,gi2Taxid)

	

#===========================
'''
objective: read ncbi phylogeny tree node file and hash it
keyword: ncib taxonomy tree hash dictionary
input:
-nodes_fn
output:
-h_tax_node2par
author: cjhong@bu.edu
date: 06/16/2013
note:
'''
def read_ncbi_phylogeny_mysql(con,downloadD):
	taxon_pyv=downloadD+'/nodes_hash.pyv'
	if os.path.isfile(taxon_pyv):
		fp=open(taxon_pyv,'r')
		h_tax_node2par = pickle.load(fp)
	else:
		#h_tax_node2par = taxon_rel_in_hash(nodes_fn)
		h_tax_node2par = dbUtils.taxon_rel_in_hash_mysql(con)
		fp=open(taxon_pyv,'w')
		pickle.dump(h_tax_node2par,fp)
	fp.close()
	return h_tax_node2par

'''
objective: read ncbi phylogeny database, store it into dictionary variable
input:
-nodes_fn: it follows the format like {gi}\t{ti1}\t{ti2}\n
output:
-h_tax_node2par: {gi1:ti1,gi2:ti2,gi3:ti1,gi5:ti2}
note: to be replaced by taxon_rel_in_hash_mysql()
'''
#=================================
def taxon_rel_in_hash(nodes_fn):
	print 'storing taxon tree to hash memory...'
	fp=open(nodes_fn,'r')
	TOP_LEVEL=-1
	[NODE, PARENT]=range(2)
	NOT_VISIT,INCLUDE,EXCLUDE = range(3)
	tax=[-1,0]
	h_tax_node2par={}
	for n in fp:
		mObj=re.search(r'(\d+)\t\|\t(\d+)\t',n)
		tax[NODE]=int(mObj.group(1))
		tax[PARENT]=int(mObj.group(2))
		
		if not h_tax_node2par.get(tax[NODE],[-1,False])[1]: #have you seen before
			h_tax_node2par[tax[NODE]]=[-1,NOT_VISIT] #first initialize

		if tax[NODE]!=tax[PARENT]:
			h_tax_node2par[tax[NODE]][0]=tax[PARENT]
		else: #only one case 1:1
			h_tax_node2par[tax[NODE]][0]=TOP_LEVEL

	fp.close()
	print 'done'
	return h_tax_node2par
	
#===========================
'''
objective: read ncbi phylogeny tree node file and hash it
keyword: ncib taxonomy tree hash dictionary
input:
-nodes_fn
output:
-h_tax_node2par
author: cjhong@bu.edu
date: 06/16/2013
note: to be replaced by read_ncbi_phylogeny_mysql()
'''
def read_ncbi_phylogeny(nodes_fn):
	taxon_py=nodes_fn+'_hash.pyv'
	if os.path.isfile(taxon_py):
		fp=open(taxon_py,'r')
		h_tax_node2par = pickle.load(fp)
	else:
		h_tax_node2par = taxon_rel_in_hash(nodes_fn)
		fp=open(taxon_py,'w')
		pickle.dump(h_tax_node2par,fp)
	fp.close()
	return h_tax_node2par

'''
objective: given fa file containing gi in front, store them into hash
'''
#=================================
def hash_gi_in_fasta(fa):
	fp=open(fa,'r')
	h_rid={}
	for r in seqParse.parse(fp,'fasta'):
		mObj=re.search(r'^gi\|(\d+)\|',r.id)
		gi=int(mObj.group(1))
		h_rid[gi]=True
	fp.close()
	return h_rid
	

'''
objective:
input:
-fa:
-catTag:
-ncbiNt_ti:
-h_tax2cat:
'''
#=================================
def register_fa_category(fa,catTag,ncbiNt_ti,h_tax2cat):
	
	h_rid=hash_gi_in_fasta(fa)
	fp=open(ncbiNt_ti,'r')
	for r in seqParse.parse(fp,'fasta'):
		mObj=re.search(r'^ti\|(\d+)\|gi\|(\d+)\|.*',r.id)
		gi=int(mObj.group(2))
		if h_rid.get(gi,-1)!=-1:
			ti=int(mObj.group(1))
			if ti!=-1:
				#print ti #debug
				h_tax2cat[ti]=catTag
	fp.close()
	return h_tax2cat

#==================================
'''
objective: given ncbi taxonomy node in hash format and taxonomy id in list, retrieve all sub subspecies or taxonomies under the taxonomy.
keyword: taxonomy, subtrees
input:
-h_tax_node2par:
-hiTaxons:

output:
-subtax: a set of sub taxon in list format

author: cjhong@bu.edu
date: 06/16/2013
note:
'''
def get_allsub_taxons_phylo(h_tax_node2par,hiTaxons, exclude_taxon_ids):
	#select all lower level taxon_id under hiTaxons of our interest and save them into hash
	NOT_VISIT,INCLUDE,EXCLUDE = range(3)
	TOP_LEVEL=-1
	TOP_LEVEL0=1
	NA = -1
	
	h_tax_node2par_upd={}
	print 'traversing phylogeny tree to get sub taxonomy ids...'
	
	for node1 in h_tax_node2par:
		if node1 == TOP_LEVEL0:
			continue
		node=node1
		nodes_in_trace=[]
		
		while True:
			nodes_in_trace.append(node)
			if (node in hiTaxons):
				search_result=INCLUDE
				break
			elif (node in exclude_taxon_ids):
				search_result=EXCLUDE
				break

			search_result = h_tax_node2par_upd.get(node,NOT_VISIT)
			if search_result==INCLUDE:#was found prev
				break
			elif search_result==EXCLUDE:#was found prev but not valid
				break
			else: # keep going upward..
				#print 'node=%d' % node
				node=h_tax_node2par.get(node,[TOP_LEVEL,EXCLUDE])[0]
				if node == TOP_LEVEL:
					search_result=EXCLUDE
					break

		for node2 in nodes_in_trace:
			h_tax_node2par_upd[node2]=search_result
	print 'done'
	
	#collect tax_id.fa in each taxon_of_interest to print out
	print 'printing fa file ...'
	subtax=[]
	for node in h_tax_node2par_upd:
		status=h_tax_node2par_upd.get(node)
		if status==INCLUDE:
			subtax.append(node)
	print 'done'
	return subtax

def check_if_nt_has_ti(nt):
	fp = open(nt,'r')
	i = fp.next().strip()
	hasTiF=False
	mObj = re.search(r'>ti\|(\d+)\|',i)
	if mObj:
			hasTiF=True
			print '%s has already taxonomy id.' % nt
	fp.close()
	return hasTiF

def append_ti_into_fasta_mysql(con, nt, Ti2sel, enable_descF, enable_onlineF,
		nt2, noTaxIdFa, invalSelFlag):
	
	NOT_VALID=-1
	GET_ALL_TAX=-2
	TAXON_ID=1
	
	#check if nt has ti tagged already
	tiReadyF=False
	if check_if_nt_has_ti(nt):
		tiReadyF=True

	get_all_taxF=False
	if Ti2sel[0]==GET_ALL_TAX:
		get_all_taxF=True

	
	print 'selecting some reference genome sequences in [%s]' % nt
	
	if (invalSelFlag):
		fp1 = open(noTaxIdFa,'w')
	with open(nt2,'w') as fp2:
		with open(nt,'r') as fp:
			for r in seqParse.parse(fp,'fasta'):
				if tiReadyF:
					mObj=re.search(r'ti\|(\d+)\|',r.id)
					if not mObj:
						continue
					ti=int(mObj.group(1))
					if ti!=NOT_VALID and (get_all_taxF or (ti in Ti2sel)):
						if enable_descF and r.description:
							fp2.write('>%s\n%s\n' % (r.description, r.seq))
						else:
							fp2.write('>%s\n%s\n' % (r.id, r.seq))
				else:
					mObj=re.search(r'gi\|(\d+)\|\S+\|(\S+)',r.id)
					if not mObj:
						continue
					gi=int(mObj.group(1))
					
					with con:
						cur=con.cursor()
						sqlcmd='select taxon from giAnnoT where gi=%d' %gi
						cur.execute(sqlcmd)
						entr = cur.fetchone()
						if entr:
							ti=int(entr[0])
						elif enable_onlineF:
							seqId=int(mObj.group(2))
							ti=pathoUtilsA.ncbi_eutil(gi,seqId,TAXON_ID) #updated ti
						else:
							ti=NOT_VALID
					
					if ti==NOT_VALID:
						if (invalSelFlag):
							fp1.write('>ti|-1|%s\n%s\n' % (r.description,r.seq))
					else:
						if get_all_taxF or (ti in Ti2sel):
							organismName, _ = dbUtils.findOrganismLineage(con, ti)
							organismName = re.sub('\s+', '_', organismName)
							if enable_descF and r.description:
								fp2.write('>ti|%d|org|%s|%s\n%s\n' % (ti, organismName, 
									r.description, r.seq))
							else:
								fp2.write('>ti|%d|org|%s|%s\n%s\n' % (ti, organismName, 
									r.id, r.seq))
	
	print 'check %s' % nt2
	if (invalSelFlag):
		fp1.close()
		print 'check %s' % noTaxIdFa
	print 'done.'
	
#=============================
'''
objective: append taxonomy id in front of each sequence header
keyword: add taxonomy id, ncbi gi nt database
input:
-nt:
-maxGi:
-gi2ti:
-Ti2sel:
-enable_onlineF:

output:
-nt2: updated reference sequence in fasta
-noTaxIdFa: fasta file where gi does not match to any offline database or online database

author: cjhong@bu.edu
date: 06/16/2013
note: to be replaced by append_ti_into_fasta_mysql()
'''
def append_ti_into_fasta_hash(nt, gi2taxFn, Ti2sel, enable_descF, enable_onlineF,
		nt2, noTaxIdFa, invalSelFlag):
	
	NOT_AVAIL=0
	NOT_VALID=-1
	GET_ALL_TAX=-2
	TAXONOMY_ID=1
	
	#check if nt has ti tagged already
	tiReadyF=False
	if check_if_nt_has_ti(nt):
		tiReadyF=True
	
	if not tiReadyF:
		(maxGi,gi2ti)=gi2tax_list(gi2taxFn)
	
	get_all_taxF=False
	if Ti2sel[0]==GET_ALL_TAX:
		get_all_taxF=True
	
	if os.path.exists(nt2):
		return (nt2,noTaxIdFa)
		
	print 'selecting some reference genome sequences in [%s]...' % nt
	
	if (invalSelFlag):
		fp1 = open(noTaxIdFa,'w')
	with open(nt2,'w') as fp2:
		with open(nt,'r') as fp:
			if tiReadyF:
				for r in seqParse.parse(fp,'fasta'):
					#print r.id #debug
					mObj=re.search(r'ti\|(\d+)\|',r.id)
					if not mObj:
						continue
					ti=int(mObj.group(1))
					if get_all_taxF or (ti in Ti2sel):
						if enable_descF and r.description:
							fp2.write('>%s\n%s\n' % (r.description, r.seq))
						else:
							fp2.write('>%s\n%s\n' % (r.id, r.seq))
			else:
				for r in seqParse.parse(fp,'fasta'):
					mObj=re.search(r'gi\|(\d+)\|\S+\|(\S+)',r.id)
					if not mObj:
						continue
					gi=int(mObj.group(1))
					if gi>maxGi or gi2ti[gi]==NOT_AVAIL:
						if enable_onlineF:
							genbank_id=mObj.group(2) #telling exactly, it must be any gene name in a database
							#genbank_id=entries[3] #telling exactly, it must be any gene name in a database
							ti=pathoUtilsA.ncbi_eutil(gi,genbank_id,TAXONOMY_ID) #updated ti
						else:
							ti=NOT_VALID
					else:
						ti=gi2ti[gi]
						
					if gi<maxGi:
						gi2ti[gi]=ti
						
					if ti==NOT_VALID:
						if invalSelFlag:
							fp1.write('>ti|-1|%s\n%s\n' % (r.description, r.seq))
					else:
						if get_all_taxF or (ti in Ti2sel):
							if enable_descF:
								fp2.write('>ti|%d|%s\n%s\n' % (ti, r.description, r.seq))
							else:
								fp2.write('>ti|%d|%s\n%s\n' % (ti, r.id, r.seq))

	print 'check %s' % nt2
	if (invalSelFlag):
		fp1.close()
		print 'check %s' % noTaxIdFa
	print 'done.'

#======================================
'''
objective: compute a distance between two nodes(upper and lower)
'''
def check_subtax_dist(utax,dtax,hPhylo):
	tax=dtax
	TOP_LEVEL=-1
	NAi=-2
	distance=0
	if utax==dtax:
		return distance
	while True:
		tax=hPhylo.get(tax,[NAi,-1])[0]
		if tax==NAi or tax==TOP_LEVEL:
			distance=NAi
			break
		elif tax==utax:
			distance+=1
			break
		else:
			distance+=1
	return distance

#=================================
'''
objective: this app is for innocentive project done in offline and to be maintained to reproduce human genome and target database for blood and tissue samples especially. Note that this is application and it may not be useful for everybody.
keyword: innocentive, blood and tissue, human genome, filtering
input:
-ncbiNt: ncbi target database. A suggested header format should be like >gi|nnnnnn|some_db|ssssss|...
-ncbiGi2TaxConvFn:
-ncbiCatFn:
-nodes_fn:
-outD:
output:
-nt_tgt_fasta: non-host genome target database for host genome
-unclutured_fasta: Unclassified categories which mainly corresponds to environment samples
-host_genome_fasta: all human genome sequences under tax_id:9606
author: cjhong@bu.edu
date: 06/18/13
note: TODO it should be updated corresponding to a new append_ti_into_fasta_app()
'''
def build_innocentive_hg19_tgt_db(ncbiNt,ncbiGi2TaxConvFn,ncbiCatFn,nodes_fn,outD):
	
	if not os.path.exists(outD):
		os.makedirs(outD)
	
	##############################
	#1. read gi2ti.dump
	##############################
	(ncbiNt_ti,ncbiNt_ti_invalid)=append_ti_into_fasta_app(ncbiNt,[-2],ncbiGi2TaxConvFn,nodes_fn,1)

	##############################
	#3. read categories.dump
	##############################
	# [things not belong to the other categories, archea, bac, eukarayote, protozoa, human, fungi, virus(viroids), unclassified, others, plasmid_transpozon]
	catIdx=['N','A','B','E','EP','EH','EF','V','U','O','OPT']
	(h_tax2cat)=getCatFromTaxid2(ncbiCatFn)

	##############################
	#4. handle protozoa fasta from ncbi ftp
	##############################
	protozoaD = outD +'/protozoa'
	if not os.path.exists(protozoaD):
		os.makedirs(protozoaD)
	#4.1 read ncbi protozoa fasta file to get taxon_id
	protozoaDfn = pathoUtilsA.ex_wget_download(protozoaD)

	#4.2 register ncbi collected protozoa and update category

	cat2append='EP'
	h_tax2cat2=register_fa_category(protozoaDfn,cat2append,ncbiNt_ti,h_tax2cat)

	##############################################################
	#5 read taxonomy tree and extract human and protozoa and fungi
	##############################################################
	h_tax_node2par=read_ncbi_phylogeny(nodes_fn)

	#6. tag homo sapiens
	hiTaxons=[9606]
	subTaxons=get_allsub_taxons_phylo(h_tax_node2par,hiTaxons)
	cat2append='EH';
	print '----------------\n'
	print cat2append
	h_tax2cat2=update_categories(subTaxons,h_tax2cat2,cat2append)

	#7 tag all protozoa sub found from web resources
	hiTaxons=[5758, 5741, 5794, 35581, 6029, 5878, 5738, 5653, 37104, 6029, 5794, 554915]
	subTaxons=get_allsub_taxons_phylo(h_tax_node2par,hiTaxons)
	cat2append='EP'
	print '----------------\n'
	print cat2append
	h_tax2cat2=update_categories(subTaxons,h_tax2cat2,cat2append)

	#8 get all fungi taxon_ids
	hiTaxons=[4751]
	subTaxons=get_allsub_taxons_phylo(h_tax_node2par,hiTaxons)

	#8.1 append all fungi taxon_ids to cat (update)
	cat2append='EF'
	print '----------------\n'
	print cat2append
	h_tax2cat2=update_categories(subTaxons,h_tax2cat2,cat2append)


	##############################################################
	#6 extracting Others(plasmid and transpozon)
	##############################################################
	hiTaxons=[2387,36549]
	subTaxons=get_allsub_taxons_phylo(h_tax_node2par,hiTaxons)
	cat2append='OPT'
	print '----------------\n'
	print cat2append
	h_tax2cat2=update_categories(subTaxons,h_tax2cat2,cat2append)

	############################################$
	#9. now, let us group nt2 into each category
	############################################$
	h_tax2cat3=group_print_ti_cat_fa(catIdx,h_tax2cat2,ncbiNt_ti,outD)

	#####################################################
	#10. merge them into two files, (EH.fa and nt_tgt.fa)
	#####################################################
	tgtD = '%s/blood_tissue' % outD
	os.makedirs(tgtD)
	nt_tgt_fasta='%s/nt_tgt.fa' % tgtD
	cmd='cd %s\ncat A.fa B.fa EP.fa EF.fa V.fa OPT.fa N.fa > %s\n' % (outD,nt_tgt_fasta)
	os.system(cmd)

	uncD = '%s/unclutured' % outD
	os.makedirs(uncD)
	cmd='mv %s/U.fa %s/' % (outD,uncD)
	os.system(cmd)
	unclutured_fasta='%s/U.fa' % uncD

	humD = '%s/9606' % outD
	os.makedirs(humD)
	cmd='mv %s/EH.fa %s/' % (outD,humD)
	os.system(cmd)
	human_fasta='%s/EH.fa' % humD

	#optional: save hash variable, tax2cat, to use later
	h_tax2cat3_pyv = ncbiCatFn+'_with_eupathodb_onine.pyv'
	with open(h_tax2cat3_pyv,'w') as f:
		pickle.dump([h_tax2cat3,catIdx],f)

	return (nt_tgt_fasta,unclutured_fasta,human_fasta)
	
#=================================
def update_categories(taxonIds,h_ti2cat,cat2append):
	for ti in taxonIds:
		print ti
		h_ti2cat[ti]=cat2append
	return h_ti2cat

	
#=================================
def get_genebank_file_gi(gi,gb_report):
	cmd='curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=gb&retmode=xml\" -o %s' % (gi,gb_report)
	os.system(cmd)
	gotIt=True
	#check if it is empty
	if os.stat(gb_report).st_size == 0: 
		gotIt=False
	return gotIt



#=================================
def csv_update_anno_gi(fp,gi,ref_name,seq_len,taxon,product,has_sub,pkey):
	fp.write('%d\t%d\t%s\t%d\t%d\t%s\t%d\n' % (pkey,gi,ref_name,seq_len,taxon,product,has_sub))
	return (pkey+1,fp)
	


#=================================
def csv_update_ti(fp,ti,organism,lineage,pkey):
	dbSz = 1
	fp.write('%d\t%d\t%d\t%s\t%s\n' % (pkey,ti,dbSz,organism,lineage))
	return (pkey+1,fp)
	

	
#======================
def mysql_str_sanity(str2):
	NAs='X'
	if not str2:
		print 'it is empty string!'
		return NAs
	str2 = re.sub('\'', '', str2)
	return str2
	

#=================================
def csv_update_delim(fp,gi,sub_gi,gene,locus_tag,protein_id,strand,stbp,edbp,pkey):
	fp.write('%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n' % (pkey,gi,sub_gi,gene,locus_tag,protein_id,strand,stbp,edbp))
	return (pkey+1,fp)
	


def getGi2TaxDump_online(downloadD):
	if not os.path.exists(downloadD):
		os.makedirs(downloadD)
	downFbase='gi_taxid_nucl.dmp'
	downExt='gz'
	gi2taxFtp='ftp://ftp.ncbi.nih.gov/pub/taxonomy/%s.%s' % (downFbase,downExt)
	gi2taxDump=pathoUtilsA.wget_download2(gi2taxFtp,downFbase,downExt,downloadD,'select','gi_taxid_nucl.dmp')
	return gi2taxDump
	
def getNodesDump_online(downloadD):
	if not os.path.exists(downloadD):
		os.makedirs(downloadD)
	downFbase='taxdump'
	downExt='tar.gz'
	phyloNodeFtp='ftp://ftp.ncbi.nih.gov/pub/taxonomy/%s.%s' % (downFbase,downExt)
	nodesDfn=pathoUtilsA.wget_download2(phyloNodeFtp,downFbase,downExt,downloadD,'select','nodes.dmp')
	return nodesDfn
	
## parses a single ncbi entry. Input is a list with the lines of an ncbi entry as elements. Output is a list of things we want from the entry.
def parse_ncbi_entry(entry): 
	gi,ref_name,taxonId,organism,lineage,product,stbp,edbp='0','','0','','','','0','0'  ## set all values of interest equal to an empty string
	sub_gis=[]
	E=len(entry)
	for x in enumerate(entry):
		tmp=x[1].strip()
		if len(tmp)<1:  # if nothing is in the line continue to the next one
			continue
		if tmp[:10] == "DEFINITION":
			if not product:
				product = tmp[12:]
				product = mysql_str_sanity(product)
		if tmp[:7] == "VERSION":
			if not ref_name:
				ref_name = tmp.split()[1]
				ref_name = mysql_str_sanity(ref_name)
				gi = tmp.split(":")[1]
			continue
		if tmp[:10] == "ORGANISM  ":
			if not organism:
				#print x[1] #debug
				organism = x[1].split("ORGANISM  ")[1][:-1]
				organism = mysql_str_sanity(organism)
				lineage = ''
				j=x[0]+1
				while (j<E) and (entry[j].strip()[:9] != "REFERENCE") and (entry[j].strip()[:8] != "FEATURES"):
					if j > (x[0]+1):
						lineage += " "
					lineage += entry[j].strip()
					j+=1
					#print entry[j] #debug
				lineage = mysql_str_sanity(lineage)
			continue
		if tmp[:7] == 'source ':
			if (stbp=='0' or edbp=='0'):
				#print tmp
				mObj=re.search(r'source\s+(\d+)\.\.(\d+)',tmp)
				if mObj:
					stbp = mObj.group(1)
					edbp = mObj.group(2)
			continue
		if tmp[:16] == '/db_xref="taxon:':
			if taxonId=='0':
				taxonId = tmp.split(":")[1][:-1]
			continue
		if tmp[:6] == "CDS   ":
			j=x[0]
			(strand_sub,cds_sub_locs)=parse_gene_location(j,entry)
			gi_sub,gene_sub,locus_tag_sub,product_sub,protein_id_sub = "0","","","",""
			cdstmp=[]
			#print j #debug
			#print '[t]%s' % tmp #debug
			while True and j<E:
				#print '[%d]%s' % (j,entry[j]) #debug
				if re.search(r'(translation|ORIGIN)',entry[j]):
					break
				else:
					cdstmp.append(entry[j].strip())
				j+=1

			for y in cdstmp:
				#print y #debug
				if y[:5] == "/gene":
					#gene_sub=y.split('="')[1][:-1]
					mObj = re.search(r'gene=[\'\"](.+)[\'\"]',y)
					if mObj:
						gene_sub = mObj.group(1)
					continue
				if y[:10] == "/locus_tag":
					mObj = re.search(r'locus_tag=\"(\S+)\"',y)
					if mObj:
						locus_tag_sub = mObj.group(1)
					continue
				if y[:8] == "/product":
					#product_sub=y.split('="')[1][:-1]
					if y[-1]=='"':
						mObj = re.search(r'product=[\'\"](.+)[\'\"]',y)
					else:
						mObj = re.search(r'product=[\'\"](.+)',y)
					if mObj:
						product_sub = mObj.group(1)
						product_sub = mysql_str_sanity(product_sub)
					continue 
				if y[:11] == "/protein_id":
					mObj = re.search(r'protein_id=[\'\"](.+)[\'\"]',y)
					if mObj:
						protein_id_sub = mObj.group(1)
					#protein_id_sub=y.split('="')[1][:-1]
					continue 
				if y[:8] == "/db_xref":
					if gi_sub=='0':
						#print y #debug
						mObj=re.search(r'GI:(\d+)',y)
						if mObj:
							gi_sub=mObj.group(1)
					continue
			if gi_sub!='0':
				for c in cds_sub_locs:
					sub_gis.append([int(gi_sub),strand_sub,int(c[0]),int(c[1]),gene_sub,locus_tag_sub,product_sub,protein_id_sub])
				
	return([int(gi),ref_name,int(taxonId),organism,lineage,product,int(stbp),int(edbp),sub_gis])

#============================$
def parse_gene_location(idx,entry):
	#check complement
	line = entry[idx]
	strand_sub='+'
	if re.search(r'complement',line):
		strand_sub='-'
	#check if it has join delimiter
	idx2 = idx
	cds_sub_locs=[]
	while True:
		line = entry[idx2]
		idx2+=1
		mObj=re.findall(r'(\d+)\.\.(\d+)',line)
		if mObj:
			for m in mObj:
				cds_sub_locs.append(list(m))
		mObj=re.findall(r'(\d+)\.\.>(\d+)',line)
		if mObj:
			for m in mObj:
				cds_sub_locs.append(list(m))
		else:
			break
	
	return (strand_sub,cds_sub_locs)

#==========================
'''
note: TODO
1)to include locus_tag
2)to download(protein) ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*.protein.gpff.gz
3)to download(rna) ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*.rna.gbff.gz
4)to download(tsa) ftp://ftp.ncbi.nih.gov/genbank/tsa/*.gbff.gz
5)add one more column to gi_annoT to indicate if this is (genbank, refSeq, tsa) and (primary|nucleotide, protein, rna), so that we can compute dbSize efficiently
'''
def gb2prepare_load_data_file(MySqlConf,tiNtDfn,downloadD):
	
	HOST_NAME,MYSQL_PORT,USER,PASSWORD,DEFAULT_DB = range(5)
	#TODO[to develop to maintain pathoDB , differential updates]----------
	#to clean up downloadD first
	downloadD_gff = downloadD + '/gbff'
	
	gbDnExt=['x','x','x']
	gbDnExt[0] = 'seq.gz'
	gbDnExt[1] = 'gbff.gz'
	gbDnExt[2] = 'protein.gpff.gz'
	
	
	if not os.path.exists(downloadD_gff):
		os.makedirs(downloadD_gff)
	else:
		cmd = 'rm -rf %s/*.gbff.gz*\n' % (downloadD_gff)
		cmd = '%srm -rf %s/*.gpff.gz*\n' % (cmd,downloadD_gff)
		cmd = '%srm -rf %s/*.seq.gz*\n' % (cmd,downloadD_gff)
		os.system(cmd)
			
	#download genabank flat format files
	downExt='gz'
	downFbase='*.seq'
	#downFbase='gbest?.seq' #debug
	gbFlatFtpD='ftp://ftp.ncbi.nih.gov/genbank'
	gbFlatFtp=gbFlatFtpD+'/'+downFbase+'.'+downExt
	
	dummy=pathoUtilsA.wget_download2(gbFlatFtp,downFbase,downExt,downloadD_gff,'nothing','X')
	
	#download refseq flat format files (genomic and rna)
	downExt='gz'
	downFbase='*.gbff'
	#downFbase='complete.?.*.gbff' #debug
	refSeqFlatFtpD='ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete'
	refSeqFlatFtp=refSeqFlatFtpD+'/'+downFbase+'.'+downExt
	dummy=pathoUtilsA.wget_download2(refSeqFlatFtp,downFbase,downExt,downloadD_gff,'nothing','X')
	
	
	#download refseq flat format files (protein)
	downExt='gz'
	downFbase='*.gpff'
	#downFbase='complete.?.protein.gpff' #debug
	refSeqFlatFtpD='ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete'
	refSeqFlatFtp=refSeqFlatFtpD+'/'+downFbase+'.'+downExt
	dummy=pathoUtilsA.wget_download2(refSeqFlatFtp,downFbase,downExt,downloadD_gff,'nothing','X')
			
	#download tsa flat format files
	downExt='gz'
	downFbase='*.gbff'
	#downFbase='tsa.GAA?.1.gbff' #debug
	refSeqFlatFtpD='ftp://ftp.ncbi.nih.gov/genbank/tsa'
	refSeqFlatFtp=refSeqFlatFtpD+'/'+downFbase+'.'+downExt
	dummy=pathoUtilsA.wget_download2(refSeqFlatFtp,downFbase,downExt,downloadD_gff,'nothing','X')
	#except the following format files
	cmd = 'rm -rf %s/*.mstr.gbff.gz' % downloadD_gff
	os.system(cmd)
	

	#the following two lines should be conistent with the one defined in parse_ncbi_entry()
	GI,REF_NAME,TAXON_ID,ORGANISM,LINEAGE,PRODUCT,STBP,EDBP,SUB_GI=range(9)
	GI_SUB,STRAND_SUB,STBP_SUB,EDBP_SUB,GENE_SUB,LOCUS_TAG_SUB,PRODUCT_SUB,PROTEIN_ID_SUB=range(8)
	
	NAi=0
	NAs='X'

	##################################################################$
	#processing genbank flat file and transfer all annotation to mysql
	##################################################################$
	ANNO_T,DELIM_T,TAX_T=range(3)
	
	gi_annoT_fn=downloadD_gff+'/giAnnoT2load.csv'
	delimT_fn=downloadD_gff+'/delimT2load.csv'
	taxT_fn=downloadD_gff+'/tax2load.csv'
	
	h_taxLookup = {}
	#read tax2load to dictionary (debug) ------------>
	if False:
		fp=open(taxT_fn,'r')
		for i in fp:
			words = i.split('\t')
			h_taxLookup[words[1]]=1
		fp.close()
	#<--------------------
	
	
	fps=[-1,-1,-1]
	fps[ANNO_T]=open(gi_annoT_fn,'w')
	fps[DELIM_T]=open(delimT_fn,'w')
	fps[TAX_T]=open(taxT_fn,'w')

	print 'transferring gene bank report to mysql...'
	gbFlatTmp = '%s/gb2process.tmp' % downloadD_gff
	pkey_anno = 1
	pkey_delim = 1
	pkey_ti = 1

	
	doneGbD = downloadD_gff+'/completed_gbff'
	if not os.path.exists(doneGbD):
		os.makedirs(doneGbD)

	#count a total # of gz to process

	F = len(os.listdir(downloadD_gff))
	f = 0
	for gbFlatFn in os.listdir(downloadD_gff):
		tick=time()
		
		if gbFlatFn.endswith(gbDnExt[0]) or gbFlatFn.endswith(gbDnExt[1]) or gbFlatFn.endswith(gbDnExt[2]):# or gbFlatFn!='gbcon208.seq.gz': #debug
			cmd='gunzip -c %s/%s > %s\n' % (downloadD_gff,gbFlatFn,gbFlatTmp)
			cmd='%smv %s/%s %s/%s\n' % (cmd,downloadD_gff,gbFlatFn,doneGbD,gbFlatFn)
			os.system(cmd)
			f+=1
		else:
			continue
		print 'processing %s[%d/%d]...' % (gbFlatFn,f,F)
		
		fp = open(gbFlatTmp,'r')
		#skipping header
		header=True
		while header: # Dump the header in the file
			tmp=fp.readline()
			if len(tmp)>5:
				header = (not tmp[:5] == "LOCUS")
		entry = [tmp]
		
		#only focus on the section between "LOCUS" ... "//"
		ti=-1
		for x in fp:
			if re.search(r'^//', x): # Every time we get to a //\n line, we read the current entry and then start collecting a new one.
				gB = parse_ncbi_entry(entry)
				#print gB[0] #debug
				entry=[]
				#.............................................
				#0) check if query gi has multiple sub gis
				has_sub=0
				if len(gB[SUB_GI])>0:
					has_sub=1

				#1) update query gi annotation
				#mysql_update_anno_gi(con,gB[GI],gB[REF_NAME],gB[EDBP],gB[TAXON_ID],gB[PRODUCT],has_sub)
				pkey_anno,fps[ANNO_T] = csv_update_anno_gi(fps[ANNO_T],gB[GI],gB[REF_NAME],gB[EDBP],gB[TAXON_ID],gB[PRODUCT],has_sub,pkey_anno)
				
				#2) update query delimit
				# mysql_update_delim(con,gB[GI],gB[GI],'+',gB[STBP],gB[EDBP])
				
				for s in gB[SUB_GI]:
					#3) update sub_gi annotation
					#mysql_update_anno_sub_gi(con,s[GI],gB[REF_NAME],gB[TAXON_ID],s[PRODUCT_SUB],s[GENE_SUB],s[PROTEIN_ID_SUB])
					#4) update sub_gi delimit
					#mysql_update_delim(con,gB[GI],s[GI_SUB],s[GENE_SUB],s[PROTEIN_ID_SUB],s[STRAND_SUB],s[STBP_SUB],s[EDBP_SUB])
					pkey_delim,fps[DELIM_T]=csv_update_delim(fps[DELIM_T],gB[GI],s[GI_SUB],s[GENE_SUB],s[LOCUS_TAG_SUB],s[PROTEIN_ID_SUB],s[STRAND_SUB],s[STBP_SUB],s[EDBP_SUB],pkey_delim)


				#ti=mysql_update_ti(con,gB[TAXON_ID],gB[ORGANISM],gB[LINEAGE])
				if h_taxLookup.get(gB[TAXON_ID],-1) == -1:
					pkey_ti,fps[TAX_T] = csv_update_ti(fps[TAX_T],gB[TAXON_ID],gB[ORGANISM],gB[LINEAGE],pkey_ti)
					h_taxLookup[gB[TAXON_ID]]=1
					
				#.............................................
			else:
				entry.append(x)
		#for loop end (x in fp)
		fp.close()
		tock=time()
		elapsed=tock-tick
		print 'elasped time:[%g]' % elapsed
	
	fps[0].close()
	fps[1].close()
	fps[2].close()
	#(gbFlatFn) finish for loop
	
	con = dbUtils.init_mysql_innocentive(MySqlConf,0)
	with con:

		print 'loading %s...' % (gi_annoT_fn)
		cur=con.cursor()
		mysql_load_cmd = 'load data local infile \'%s\' into table giAnnoT fields terminated by \'\\t\'' % gi_annoT_fn
		cur.execute(mysql_load_cmd)
		cur=con.cursor()
		mysql_idx_cmd = 'create unique index idx_gi on giAnnoT (gi)'
		cur.execute(mysql_idx_cmd)
		print 'done.'

		print 'loading %s...' % (delimT_fn)
		cur=con.cursor()
		mysql_load_cmd = 'load data local infile \'%s\' into table giDelimT fields terminated by \'\\t\'' % delimT_fn
		cur.execute(mysql_load_cmd)
		cur=con.cursor()
		mysql_idx_cmd = 'create index idx_subgi on giDelimT (gi,stbp,edbp)'
		cur.execute(mysql_idx_cmd)
		print 'done.'
		
		print 'computing database size for each taxon id...'
		if False:
			#collect dbSize for each ti
			h_ti_dbSz = get_ti_db_size(gi_annoT_fn)
			update_taxT_fn(h_ti_dbSz,taxT_fn)
		else:
			add_dbsize2taxonT(tiNtDfn,taxT_fn)
		print 'done.'
		
		print 'loading %s...' % (taxT_fn)
		cur=con.cursor()
		mysql_load_cmd = 'load data local infile \'%s\' into table cj_taxonT fields terminated by \'\\t\'' % taxT_fn
		cur.execute(mysql_load_cmd)
		cur=con.cursor()
		mysql_idx_cmd = 'create unique index idx_taxon on cj_taxonT (taxon)'
		cur.execute(mysql_idx_cmd)
		print 'done.'
		
	dbUtils.mysql_close(con)
	print 'done'

#=====================
def get_dbsize_tiNtFa(tiNtFa):
	print 'computing searching space size for each taxonomy id...'
	h_tiDbSize = {}
	fp = open(tiNtFa,'r')
	ti = 0
	for r in fp:
		if r[0]=='>':
			mObj = re.search(r'>ti\|(\S+)\|gi',r)
			ti = mObj.group(1)
		else:
			h_tiDbSize[ti] = h_tiDbSize.get(ti,0) + len(r)
	fp.close()
	print 'done.'
	return h_tiDbSize
	
#==================
def add_dbsize2taxonT(tiNtDfn,taxT2load_fn):
	#check if tiNtDfn has pyv file of dictionary (ti:dbSize)
	tiNt_pyv = pathoUtilsA.file_tag(tiNtDfn,'dbSize','pyv')
	if not os.path.exists(tiNt_pyv):
		h_tiDbSize = get_dbsize_tiNtFa(tiNtDfn)
		with open(tiNt_pyv,'w') as f:
			pickle.dump(h_tiDbSize,f)
	else:
		with open(tiNt_pyv) as f:
			h_tiDbSize=pickle.load(f)
			
	fp = open(taxT2load_fn,'r')
	taxT2load_fn2 = taxT2load_fn+'.tmp'
	fp2 = open(taxT2load_fn2,'w')
	for i in fp:
		words = i.split('\t')
		tiDbSz = h_tiDbSize.get(words[1],1)
		words[2] = str(tiDbSz)
		j = '\t'.join(words)
		fp2.write('%s' % j)
	fp2.close()
	fp.close()
	os.rename(taxT2load_fn2,taxT2load_fn)
	
#===============
def update_taxT_fn(h_ti_dbSz,taxT_fn):
	fp = open(taxT_fn,'r')
	taxT_fn2 = taxT_fn+'.tmp'
	fp2 = open(taxT_fn2,'w')
	for i in fp:
		ent=i.split('\t')
		fp2.write('%s\t%s\t%d\t%s\t%s' % (ent[0],ent[1],h_ti_dbSz.get(ent[1],1),ent[2],ent[3]))
	fp.close()
	fp2.close()
	os.rename(taxT_fn2,taxT_fn)
	
#====================
def get_ti_db_size(gi_annoT_fn):
	h_ti_dbSz={}
	fp = open(gi_annoT_fn,'r')
	for i in fp:
		ent=i.split('\t')
		ti=ent[4]
		h_ti_dbSz[ti] = h_ti_dbSz.get(ti,0) + int(ent[3])
	fp.close()
	return h_ti_dbSz


#===================================
def group_print_ti_cat_fa(catIdx,h_tax2cat,nt_ti,outD):
	#catIdx=['N','A','B','E','EP','EH','EF','V','U','O','OPT']
	h_cFps={}
	tmpFn=outD+'/online_tax_cat.tmp'
	ncbiLineage=['Fungi','Protozoan','Archaea','Bacteria','Viroids','Viruses','other sequences','unclassified sequences','Eukaryota']
	h_lin2catIdx={'Archaea':'A','Bacteria':'B','Eukaryota':'E','Viroids':'V','Viruses':'V','other sequences':'O','unclassified sequences':'U','Fungi':'EF','Protozoan':'EP','N':'N'}
	
	for c in catIdx:
		cFn=outD+'/'+c+'.fa'
		h_cFps[c]=open(cFn,'w')

	print 'grouping sequence into each cat...'
	
	NOT_AVAIL='X'
	fp = open(nt_ti,'r')
	for r in seqParse.parse(fp,'fasta'):
		mObj=re.search(r'^ti\|(\d+)\|',r.id)
		ti=int(mObj.group(1))
		#print ti #debug
		cat=h_tax2cat.get(ti,NOT_AVAIL)
		if cat==NOT_AVAIL:
			#search in phylogeny tree
			tmp = pathoUtilsA.search_cat_in_online_taxonomy(ti,ncbiLineage,tmpFn)
			cat=h_lin2catIdx.get(tmp)
			h_tax2cat[ti]=cat
		
		fp2=h_cFps.get(cat)
		fp2.write('>%s\n%s\n' % (r.id, r.seq))

	fp.close()
	if os.path.exists(tmpFn):
		os.remove(tmpFn)
	
	for c in catIdx:
		(h_cFps.get(c)).close()
		
	print 'done.'
	return (h_tax2cat)

#=================================
def getCatFromTaxid2(ncbiCatFn):

	loadSuffix='cat5'
	file2sv=ncbiCatFn+'_'+loadSuffix

	if os.path.isfile(file2sv):
		with open(file2sv,'r') as f:
			h_ti2cat=pickle.load(f)
	else:
		h_ti2cat={}
		fp=open(ncbiCatFn,'r')
		for i in fp:
			entry=i.split('\t')
			h_ti2cat[int(entry[1])]=entry[0]
			if entry[1]!=entry[2]:
				h_ti2cat[int(entry[2])]=entry[0]
		fp.close()
		with open(file2sv,'w') as f:
			pickle.dump(h_ti2cat,f)
	print 'done!'
	
	return (h_ti2cat)

# the following module contains mysql

