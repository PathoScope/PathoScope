#!/usr/bin/python
# build a ncbi phylogeny in mysql
# build genbank flat(DB: genbank, refSeq, TSA) into mysql
# example: @oligomer, try ./app_gb2mysql.py -g /media/sleepysilver/data/innocentive_id_org/db/genomes/tgt_nhost_9606/blood_tissue/nt_tgt.fa -d /media/sleepysilver/data/innocentive_id_org/db/ncbi -m localhost -u hong -p oligomer -r 0

import os, sys, argparse
import MySQLdb as mdb

pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir)

from pathoscope.patholib import pathoLib
from pathoscope.pathodb import dbUtils

#---------------------------------
parser = argparse.ArgumentParser(description="retrieve mysql from ncbi gene bank")
parser.add_argument('-g', action='store', dest='ti_nt', required=True, help='specify a fasta format reference appended with taxonomy id in the front seq header. Refer to app_build_nt_tgt.py to build such fasta file.')
parser.add_argument('-d', action='store', dest='downloadD', required=True, help='specify a temporary download directory')
parser.add_argument('-m', action='store', dest='hostname', required=False, default='localhost', help='specify hostname running mysql')
parser.add_argument('-P', action='store', dest='port', required=False, type=int, default='3306', help='specify hostname running mysql')
parser.add_argument('-u', action='store', dest='user', required=False, default='root', help='user name to access mysql')
parser.add_argument('-p', action='store', dest='passwd', required=True, default='X', help='provide password associate with user')
parser.add_argument('-r', action='store', dest='reset_table', required=False, type=int, default=0, help='set to 1 if you want to reset mysql table')

args=parser.parse_args()

#####################$
#open mysql connection
#####################$
HOST_NAME,MYSQL_PORT,USER,PASSWORD,DEFAULT_DB = range(5)
MySqlConf=[args.hostname,args.port,args.user,args.passwd,'']
con = dbUtils.init_mysql_innocentive(MySqlConf,args.reset_table)

#####################################$
#create ncbi phylogeny tree into mysql
#####################################$
if True: #debug
	nodesDfn=pathoLib.getNodesDump_online(args.downloadD)
	dbUtils.phylo_node2mysql(con,nodesDfn)

###################################$
#close mysql connection
###################################$
dbUtils.mysql_close(con)

####################################$
#create ncbi genbank flat into mysql
####################################$
#pathoLib.gb2mysql2(con,args.downloadD)
pathoLib.gb2prepare_load_data_file(MySqlConf,args.ti_nt,args.downloadD)

