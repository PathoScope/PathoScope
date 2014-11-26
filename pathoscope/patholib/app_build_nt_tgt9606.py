#!/usr/bin/python
# read alignment file
# adapted from nt2_5cat.py

import os, sys, math, argparse
import psLib

#main()
#---------------------------------
parser = argparse.ArgumentParser(description="generating database for innocentive metagenomics project")
parser.add_argument('-g', action='store', dest='reference', required=True, help='specify reference genome')
parser.add_argument('-t', action='store', dest='gi2tax_dmp', required=True, default='X', help='specify gi to taxon mapping file (Refer to ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)')
parser.add_argument('-c', action='store', dest='cat_dmp', required=True, default='X', help='specify ti to categories (Refer to ftp://ftp.ncbi.nih.gov/pub/taxonomy/categories.dmp.gz)')
parser.add_argument('-n', action='store', dest='nodes_dmp', required=True, default='X', help='specify ncbi phylogeny tree file (Refer to ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip)')
parser.add_argument('-o', action='store', dest='out_directory', required=True, default='X', help='specify output directory')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args=parser.parse_args()
(tgt_fa,unclass_fa,host_fa)=psLib.build_innocentive_hg19_tgt_db(args.reference,args.gi2tax_dmp,args.cat_dmp,args.nodes_dmp,args.out_directory)

print 'check %s' % tgt_fa
print 'check %s' % unclass_fa
print 'check %s' % host_fa
