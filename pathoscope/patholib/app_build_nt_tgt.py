#!/usr/bin/python
#

import os, sys, math, argparse
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir) 
from pathoscope.patholib import pathoLib

#=================================
#example:
#main()
#---------------------------------
parser = argparse.ArgumentParser(description="selecting taxon level genome sequences")
parser.add_argument('-g', action='store', dest='reference', required=True, help='specify reference genome(Download ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz)')
parser.add_argument('-t', action='store', dest='taxon_ids', required=False, default='X', help='specify taxon ids of your interest with comma separated (if you have multiple taxon ids). If you do not specify this option, it will work on all entries in the reference file. For taxonomy id lookup, refer to http://www.ncbi.nlm.nih.gov/taxonomy')
parser.add_argument('--desc', action='store_const', dest='desc', required=False, const=True, default=False, help='to keep an additional description in original fasta seq header. Depending on NGS aligner, a long sequence header may slow down its mapping process.')
parser.add_argument('--subtax', action='store_const', dest='subtax', required=False, const=True, default=False, help='to include all sub taxonomies under the query taxonomy id. e.g., if you set -t 4751 --subtax, it will cover all sub taxonomies under taxon id 4751 (fungi).')
parser.add_argument('--online', action='store_const', dest='online_search', required=False, const=True, default=False, help='to enable online searching in case you cannot find a correct taxonomy id for a given gi. When there are many entries in nt whose gi is invalid, this option may slow down an overall process.')
parser.add_argument('-m', action='store', dest='hostname', required=False, default='oligomer.bumc.bu.edu', help='specify hostname running mysql if you want to use mysql instead of hash method in mapping gi to taxonomy id')
parser.add_argument('-P', action='store', dest='myport', required=False, type=int,  default=3306, help='specify port running mysql')
parser.add_argument('-u', action='store', dest='user', required=False, default='pathoscope', help='user name to access mysql')
parser.add_argument('-p', action='store', dest='passwd', required=False, default='johnsonlab', help='provide password associate with user')

parser.add_argument('-o', action='store', dest='label', required=True, help='provide a short description (w/o any white space) of the target database you are about to build')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args=parser.parse_args()

################################################$
#append taxon id in the front of sequence header
################################################$
NAs = 'X'
if args.user!=NAs and args.user==NAs:
	print 'if you want to use mysql, make sure that you install pathoDB and also specify the corresponding mysql password correclty(Ask to your mysql admin to access the database).'
MysqlConf=(args.hostname,args.myport,args.user,args.passwd,'information_schema')
taxon_ids=pathoLib.parse_input_app_build_nt_tgt(args.taxon_ids)
(ncbiNt_ti,ncbiNt_invalid) = pathoLib.append_ti_into_fasta_app(args.reference,taxon_ids,args.subtax,MysqlConf,args.desc,args.online_search,args.label)
