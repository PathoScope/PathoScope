#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Functions to query from MySQL database

#	Pathoscope - Predicts strains of genomes in Nextgen seq alignment file (sam/bl8)
#	Copyright (C) 2013  Johnson Lab - Boston University
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of thefrom License, or
#	any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################################$
# Which def can be written here?
# anything using mysql connection to pathoDB. But some batch job or application is not welecome to be here. For example two applications using mysql can be found in patholib/pathoLib.py. Here place is good for atomic-level module
#########################################################################################$

import os, pickle
try:
	import MySQLdb as mdb
except ImportError:
	print "Running without mySQLdb library"

#===============================================
'''
objective: This function retrieves the organism name for the given ti(taxon id) from mysql
'''
def findOrganismLineage(con, ti):
	
	with con:
		cur = con.cursor()
		mysql_sel_cmd = 'select organism, lineage from cj_taxonT where taxon = %s' % ti
		cur.execute(mysql_sel_cmd)
		entr = cur.fetchone()
		if entr:
			organism = entr[0]
			lineage = entr[1]
	return organism, lineage

#=================================
def mysql_update_anno_gi2(con,gi,ref_name,seq_len,taxon,product,has_sub):
	with con:
		cur=con.cursor()
		mysql_ins_cmd="insert into giAnnoT(gi,ref_name,seq_len,taxon,product,has_sub) values(%d,\'%s\',%d,%d,\'%s\',%d)" % (gi,ref_name,seq_len,taxon,product,has_sub)
		#print mysql_ins_cmd #debug
		cur.execute(mysql_ins_cmd)
		
		
#=================================
# note: this module to be replaced by mysql_update_anno_gi2 soon
def mysql_update_anno_gi(con,gi,ref_name,seq_len,taxon,product,has_sub):
	NAi=0
	with con:
		cur=con.cursor()
		mysql_sel_cmd="select ref_name,seq_len,taxon,product,has_sub from giAnnoT where gi = %d limit 1" % gi
		cur.execute(mysql_sel_cmd)
		entries = cur.fetchone()
		
		if entries:
			set_cmd=[]
			if entries[0]==None and not ref_name:
				tmp='ref_name = %s' % ref_name; set_cmd.append(tmp)
			if entries[1]==None and seq_len!=NAi:
				tmp='seq_len = %d' % seq_len; set_cmd.append(tmp)
			if entries[2]==None and taxon!=NAi:
				tmp='taxon = %d' % taxon; set_cmd.append(tmp)
			if entries[3]==None and not product:
				tmp='product = %s' % product; set_cmd.append(tmp)
			if entries[4]==None and has_sub!=NAi:
				tmp='has_sub = %d' % has_sub; set_cmd.append(tmp)
			if set_cmd:
				set_cmd_str = ",".join(set_cmd)
				mysql_upd_cmd="update giAnnoT set %s where gi = %d" % (set_cmd_str,gi)
				#print mysql_upd_cmd #debug
				cur.execute(mysql_upd_cmd)
		else:
			mysql_ins_cmd="insert into giAnnoT(gi,ref_name,seq_len,taxon,product,has_sub) values(%d,\'%s\',%d,%d,\'%s\',%d)" % (gi,ref_name,seq_len,taxon,product,has_sub)
			#print mysql_ins_cmd #debug
			cur.execute(mysql_ins_cmd)
			
			
#=================================
def mysql_update_ti(con,ti,organism,lineage):

	#NAs='X'
	
	with con:
		cur = con.cursor()
		mysql_sel_cmd="select lineage from cj_taxonT where taxon = %d limit 1" % ti
		cur.execute(mysql_sel_cmd)
		entries = cur.fetchone()
		if entries: #there is a prev record and note that this case will correspond to the case only single gi
			set_cmd=''
			if entries[0]==None and not lineage:
				set_cmd='lineage = %s' % lineage
			if set_cmd:
				mysql_upd_cmd="update cj_taxonT set %s where taxon = %d" % (set_cmd,ti)
				cur.execute(mysql_upd_cmd)
		else: #first time so it will register
			mysql_ins_cmd="insert into cj_taxonT(taxon,organism,lineage) values(%d,\'%s\',\'%s\')" % (ti,organism,lineage)
			#print mysql_ins_cmd #debug
			cur.execute(mysql_ins_cmd)
	return ti
	
#=================================
def mysql_update_delim(con,gi,sub_gi,gene,locus_tag,protein_id,strand,stbp,edbp):
	#ITSELF=-1
	#NAi=0
	#if gi == sub_gi:
	#	sub_gi = ITSELF

	with con:
		cur = con.cursor()
		#first check if gi is registered prevly,
		#mysql_sel_cmd="select id, strand, stBp, edBp from giDelimT where gi = %d and sub_gi = %d limit 1" % (gi,sub_gi)
		#cur.execute(mysql_sel_cmd)
		#e2 = cur.fetchone()
		#if not e2: #wait one more step, if gi_child is also registered among these gi entries
		mysql_ins_cmd="insert into giDelimT(gi,sub_gi,gene,locus_tag,protein_id,strand,stBp,edBp) values(%d,%d,\'%s\',\'%s\',\'%s\',\'%s\',%d,%d)" % (gi,sub_gi,gene,locus_tag,protein_id,strand,stbp,edbp)
		cur.execute(mysql_ins_cmd)
	return
	
#=================================
def init_mysql_innocentive(MySqlConf,reset_table):
	HOST_NAME,MYSQL_PORT,USER,PASSWORD,DEFAULT_DB = range(5)
	inno_db = 'pathodb'
	if reset_table==1:
		MySqlConf[DEFAULT_DB]='information_schema'
		con = ''

	try:
		con = mdb.connect(host=MySqlConf[HOST_NAME],port=MySqlConf[MYSQL_PORT],user=MySqlConf[USER],passwd=MySqlConf[PASSWORD],db=MySqlConf[DEFAULT_DB])
		#con = mdb.connect(hostname,user,passwd)
	except:
		print 'DB connection error. Make sure that you install pathoDB and your mysql credential is correct.'
		con = ''
		return con
	with con:
		cur = con.cursor()
		if reset_table==0:
			query ='use %s' % inno_db
			print query
			cur.execute(query)
			return con

		query = 'create database IF NOT EXISTS %s' % inno_db
		cur.execute(query)
		query ='use %s' % inno_db
		print query
		cur.execute(query)
		
		#------------------
		#create giDelimT that contains multiple genes in a long chromosome level sequence
		cur.execute("drop table if exists giDelimT")
		cur.execute("CREATE TABLE giDelimT (id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, gi INT NOT NULL, sub_gi INT, gene VARCHAR(64), locus_tag VARCHAR(64), protein_id VARCHAR(64), strand CHAR(1), stbp INT, edbp INT)") #, INDEX idxGis(gi,sub_gi)
		
		#------------------
		#create gi annotation table
		cur.execute("drop table if exists giAnnoT")
		cur.execute("create table giAnnoT (id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, gi INT NOT NULL, ref_name VARCHAR(64), seq_len INT, taxon INT, product VARCHAR(1024), has_sub INT)") #, INDEX idxGi (gi)
		#-------------------
		#create taxon_table
		
		cur.execute("drop table if exists cj_taxonT")
		cur.execute("create table cj_taxonT (id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, taxon INT NOT NULL, dbSize BIGINT UNSIGNED, organism VARCHAR(128), lineage VARCHAR(2048))") #, INDEX idxTaxon (taxon)

		#------------------
		cur.execute("drop table if exists phylo_nodesT")
		cur.execute("create table phylo_nodesT (id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, tax_id INT NOT NULL, parent_tax_id INT NOT NULL, rank VARCHAR(16), embl_code VARCHAR(16), division_id VARCHAR(8), inherited_div_flag INT, genetic_code_id INT, inherited_GC_flag INT, mitochondrial_genetic_code_id INT, inherited_MGC_flag INT, GenBank_hidden_flag INT, hidden_subtree_root_flag INT, comments VARCHAR(128))") #, INDEX idxTax_id (tax_id)

	return con
	
	
#====================================================================
'''
objective: download ncbi phylogeny database and build mysql database
input:
-con:mysql db connection
-downloadD: a temporary space to download ncbi dataset
output: update phylo_nodesT table
'''
def phylo_node2mysql(con,nodesDfn):

	with con:
		cur = con.cursor()
		fp = open(nodesDfn,'r')
		for i in fp:
			j=i.split('\t|\t')
			mysql_ins_cmd="insert into phylo_nodesT(tax_id,parent_tax_id,rank,embl_code,division_id,inherited_div_flag,genetic_code_id,inherited_GC_flag,mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments) values (%d,%d,\'%s\',\'%s\',%d,%d,%d,%d,%d,%d,%d,%d,\'%s\')" % (int(j[0]),int(j[1]),j[2],j[3],int(j[4]),int(j[5]),int(j[6]),int(j[7]),int(j[8]),int(j[9]),int(j[10]),int(j[11]),j[12])
			cur.execute(mysql_ins_cmd)
		fp.close()
		mysql_idx_cmd = 'create index idx_tax on phylo_nodesT (tax_id,parent_tax_id)'
		cur = con.cursor()
		cur.execute(mysql_idx_cmd)

#=================================
def select_ti_db_size(MySqlConf,csv2save):
	#HOST_NAME,MYSQL_PORT,USER,PASSWORD,DEFAULT_DB = range(5)
	
	if not os.path.exists(csv2save):
		con = init_mysql_innocentive(MySqlConf,0) #DON'T REMOVE LAST 0 or CHANGE TO 1
		selCmd='select taxon, dbSize from cj_taxonT'
		with con:
			cur = con.cursor()
			cur.execute(selCmd)
			entr = cur.fetchall()
			fp2 = open(csv2save,'w')
			for j in entr:
				fp2.write('%s\t%s\n' % (j[0],j[1]))
			fp2.close()
		mysql_close(con)

	h_tiDbSz = {}
	fp = open(csv2save,'r')
	for i in fp:
		ent = i.split('\t')
		h_tiDbSz[ent[0]]=int(ent[1].rstrip())
	fp.close()
	return h_tiDbSz
	
	
#=============================
def mysql_close(con):
	if con and con.open:
		cursor=con.cursor()
		cursor.close()
		con.close()
		print 'close db'

		
#=================================
def get_refname_product_mysql(con, gis):
	NAs = 'X'
	with con:
		for g in gis:
			ref_names = []
			products = []
			selCmd = 'select ref_name, product from giAnnoT where gi = %s' % g
			cur = con.cursor()
			cur.execute(selCmd)
			entr = cur.fetchone()
			if entr:
				ref_names.append(entr[0])
				products.append(entr[0])
			else:
				ref_names.append(NAs)
				products.append(NAs)
	return (ref_names,products)
	
	
#=================================
def taxon_rel_in_hash_mysql(con, workD):
	
	print 'storing taxon tree to hash memory...'
	
	TOP_LEVEL=-1
	[NODE, PARENT]=range(2)
	#NOT_VISIT,INCLUDE,EXCLUDE = range(3)
	(NOT_VISIT,_,_) = range(3)
	tax=[-1,0]
	
	node2par_pyv=workD+'/nodes.dmp_hash.pyv'
	if os.path.exists(node2par_pyv):
		with open(node2par_pyv) as f:
			h_tax_node2par=pickle.load(f)
			return h_tax_node2par

	h_tax_node2par={}
	with con:
		cur = con.cursor()
		mysql_sel_cmd="select tax_id, parent_tax_id from phylo_nodesT"
		cur.execute(mysql_sel_cmd)
		entr=cur.fetchall()
		if entr:
			for j in entr:
				tax[NODE]=int(j[0])
				tax[PARENT]=int(j[1])

				if not h_tax_node2par.get(tax[NODE],[-1,False])[1]: #have you seen before
					h_tax_node2par[tax[NODE]]=[-1,NOT_VISIT] #first initialize

				if tax[NODE]!=tax[PARENT]:
					h_tax_node2par[tax[NODE]][0]=tax[PARENT]
				else: #only one case 1:1
					h_tax_node2par[tax[NODE]][0]=TOP_LEVEL

		else:
			print 'check your taxonomy table in pathodb database!'

		with open(node2par_pyv,'w') as f:
			pickle.dump(h_tax_node2par,f)
		
	return h_tax_node2par
	
	
#=================================
def mysql_update_anno_sub_gi(con,gi,ref_name,taxon,product,gene,protein_id):
	#NAs='X'
	NAi=0
	with con:
		cur=con.cursor()
		mysql_sel_cmd="select ref_name,taxon,product,gene,protein_id from giAnnoT where gi = %d" % gi
		cur.execute(mysql_sel_cmd)
		entries = cur.fetchone()
		
		if entries:
			set_cmd=[]
			if entries[0]==None and not ref_name:
				tmp='ref_name = %s' % ref_name; set_cmd.append(tmp)
			if (entries[1]==None or entries[1]=='0') and taxon!=NAi:
				tmp='taxon = %d' % taxon; set_cmd.append(tmp)
			if entries[2]==None and not product:
				tmp='product = %s' % product; set_cmd.append(tmp)
			if entries[3]==None and not gene:
				tmp='gene = %s' % product; set_cmd.append(tmp)
			if entries[4]==None and not protein_id:
				tmp='protein_id = %s' % protein_id; set_cmd.append(tmp)
			if set_cmd:
				set_cmd_str = ",".join(set_cmd)
				mysql_upd_cmd="update giAnnoT set %s where gi = %d" % (set_cmd_str,gi)
				#print mysql_upd_cmd #debug
				cur.execute(mysql_upd_cmd)
		else:
			mysql_ins_cmd="insert into giAnnoT(gi,ref_name,taxon,product,gene,protein_id) values(%d,\'%s\',%d,\'%s\',\'%s\',\'%s\')" % (gi,ref_name,taxon,product,gene,protein_id)
			#print mysql_ins_cmd #debu
			cur.execute(mysql_ins_cmd)
			
