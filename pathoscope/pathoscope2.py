#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Wrapper file for the following modules:
# patholib: generates host/target genome libraries from ncbi nt database for given taxon IDs 
# pathomap: aligns reads to host/target database independent of read type using Bowtie2
# pathoid: reassigns ambiguous reads to the correct genome using statistical models
# pathoreport: Writes sam files to xml format 

#usage information: pathoscope.py -h

#       Pathoscope 2.0 - Predicts strains of genomes in unassembled Nextgen seq data
#       Copyright (C) 2013  Johnson Lab - Boston University
#
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, sys
pathoscopedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,pathoscopedir) 

import argparse
from pathoscope.patholib import pathoLib
from pathoscope.pathomap import PathoMapA
from pathoscope.pathoid import PathoID
from pathoscope.pathoreport import PathoReportA
from time import time

# ===========================================================
# main ()
parser = argparse.ArgumentParser(description="Pathoscope")

# create the top-level parser
parser.add_argument('--version', action='version', version='%(prog)s 2.0.6')
parser.add_argument('--verbose', action='store_const', dest='verbose',
	required=False, const=True, default=False, help='Prints verbose text while running')
subparsers = parser.add_subparsers(dest='subcommand', help='Select one of the following sub-commands')

# create the parser for the "LIB" command
parser_a = subparsers.add_parser('LIB', help='Pathoscope taxon level reference genome Library creation Module')
parser_a.add_argument('-genomeFile', action='store', dest='lib_reference', required=True, 
	help='Specify reference genome(Download ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz)')
parser_a.add_argument('-taxonIds', action='store', dest='lib_taxon_ids', required=False, default='X', 
	help='Specify taxon ids of your interest with comma separated '
	'(if you have multiple taxon ids). If you do not specify this option, '
	'it will work on all entries in the reference file. For taxonomy id lookup, '
	'refer to http://www.ncbi.nlm.nih.gov/taxonomy')
parser_a.add_argument('-excludeTaxonIds', action='store', dest='lib_exclude_taxon_ids', 
	required=False, default='X', 
	help='Specify taxon ids to exclude with comma separated '
	'(if you have multiple taxon ids to exclude).')
parser_a.add_argument('--noDesc', action='store_const', dest='lib_nodesc', required=False, const=True,
	default=False, help='Do not keep an additional description in original fasta seq header.' 
	'Depending on NGS aligner, a long sequence header may slow down its mapping process.')
parser_a.add_argument('--subTax', action='store_const', dest='lib_subtax', required=False, 
	const=True, default=False, help='To include all sub taxonomies under the query taxonomy id.'
	' e.g., if you set -t 4751 --subtax, it will cover all sub taxonomies under taxon id 4751 '
	'(fungi).')
parser_a.add_argument('--online', action='store_const', dest='lib_online_search', required=False, 
	const=True, default=False, help='To enable online searching in case you cannot find a '
	'correct taxonomy id for a given gi. When there are many entries in nt whose gi is invalid, '
	'this option may slow down the overall process.')
parser_a.add_argument('-dbhost', action='store', dest='lib_dbhost', required=False, 
	default='localhost', help='specify hostname running mysql if you want to use mysql '
	'instead of hash method in mapping gi to taxonomy id')
parser_a.add_argument('-dbport', action='store', dest='lib_dbport', required=False, 
	default=3306, type=int, help='provide mysql server port if different from default (3306)')
parser_a.add_argument('-dbuser', action='store', dest='lib_dbuser', required=False, default='X', 
	help='user name to access mysql')
parser_a.add_argument('-dbpasswd', action='store', dest='lib_dbpasswd', required=False, default='X', 
	help='provide password associate with user')
parser_a.add_argument('-db', action='store', dest='lib_db', required=False, default='pathodb', 
	help='mysql pathoscope database name (default: pathodb)')
parser_a.add_argument('-outDir', action='store', default='.', dest='lib_outdir',
	help='Output Directory (Default=. (current directory))')
parser_a.add_argument('-outPrefix', action='store', dest='lib_outprefix', required=True, 
	help='specify an output prefix to name your target database')

# create the parser for the "MAP" command
parser_b = subparsers.add_parser('MAP', help='Pathoscope MAP Module')
parser_b.add_argument('-U', default='', action='store', dest='map_inputread', required=False, 
	help='Input Read Fastq File (Unpaired/Single-end)')
parser_b.add_argument('-1', default='', action='store', dest='map_inputread1', required=False, 
	help='Input Read Fastq File (Pair 1)')
parser_b.add_argument('-2', default='', action='store', dest='map_inputread2', required=False, 
	help='Input Read Fastq File (Pair 2)')
parser_b.add_argument('-targetRefFiles', default='', action='store', 
	dest='map_targetref', required=False, 
	help='Target Reference Genome Fasta Files Full Path (Comma Separated)')
parser_b.add_argument('-filterRefFiles', default='', action='store', 
	dest='map_filterref', required=False, 
	help='Filter Reference Genome Fasta Files Full Path (Comma Separated)')
parser_b.add_argument('-targetAlignParams', action='store', 
	dest='map_targetalignparams', default=None, required=False, 
	help='Target Mapping Bowtie2 Parameters (Default: Pathoscope chosen best parameters)')
parser_b.add_argument('-filterAlignParams', action='store', 
	dest='map_filteralignparams', default=None, required=False, 
	help='Filter Mapping Bowtie2 Parameters (Default: Use the same Target Mapping Bowtie2 parameters)')
parser_b.add_argument('-outDir', action='store', default='.', 
	dest='map_outdir', required=False, 
	help='Output Directory (Default=. (current directory))')
parser_b.add_argument('-outAlign', action='store', default='outalign.sam', 
	dest='map_outalign', required=False, 
	help='Output Alignment File Name (Default=outalign.sam)')
parser_b.add_argument('-indexDir', action='store', default='.', 
	dest='map_indexdir', required=False, 
	help='Index Directory (Default=. (current directory))')
parser_b.add_argument('-targetIndexPrefixes', default='', action='store', 
	dest='map_targetindex', required=False, 
	help='Target Index Prefixes (Comma Separated)')
parser_b.add_argument('-filterIndexPrefixes', default='', action='store', 
	dest='map_filterindex', required=False, 
	help='Filter Index Prefixes (Comma Separated)')
parser_b.add_argument('-targetAlignFiles', default='', action='store', 
	dest='map_targetalign', required=False, 
	help='Target Alignment Files Full Path (Comma Separated)')
parser_b.add_argument('-filterAlignFiles', default='', action='store', 
	dest='map_filteralign', required=False, 
	help='Filter Alignment Files Full Path (Comma Separated)')
parser_b.add_argument('-btHome', default=None, action='store', 
	dest='map_bthome', required=False, 
	help='Full Path to Bowtie2 binary directory (Default: Uses bowtie2 in system path)')
parser_b.add_argument('-numThreads', action='store', dest='map_numthreads', required=False, 
	default=8, type=int, help='Number of threads to use by aligner (bowtie2) if different from default (8)')
parser_b.add_argument('-expTag', action='store', default='pathomap', dest='map_exp_tag',
	help='Experiment Tag added to files generated for identification (Default: pathomap)')

# create the parser for the "ID" command
parser_c = subparsers.add_parser('ID', help='Pathoscope ID Module')
parser_c.add_argument('--outMatrix', action='store_true', default=False, dest='id_out_matrix',
	help='Output alignment matrix')
parser_c.add_argument('--noUpdatedAlignFile', action='store_true', default=False, 
	dest='id_noalign', help='Do not generate an updated alignment file')
parser_c.add_argument('--noDisplayCutoff', action='store_true', default=False, 
	dest='id_nocutoff', help='Do not cutoff display of genomes, even if it is insignificant')
parser_c.add_argument('-scoreCutoff', action='store', default=0.01, type=float,
	dest='id_score_cutoff', help='Score Cutoff')
parser_c.add_argument('-emEpsilon', action='store', default=1e-7, type=float,
	dest='id_emEpsilon', help='EM Algorithm Epsilon cutoff')
parser_c.add_argument('-maxIter', action='store', default=50, type=int,
	dest='id_maxIter', help='EM Algorithm maximum iterations')
parser_c.add_argument('-piPrior', action='store', default=0, type=int, dest='id_piPrior', 
	help='EM Algorithm Pi Prior equivalent to adding n unique reads (Default: n=0)')
parser_c.add_argument('-thetaPrior', action='store', default=0, type=int, dest='id_thetaPrior', 
	help='EM Algorithm Theta Prior equivalent to adding n non-unique reads (Default: n=0)')
parser_c.add_argument('-expTag', action='store', default='pathoid', dest='id_exp_tag',
	help='Experiment tag added to output file for easy identification (Default: pathoid)')
parser_c.add_argument('-outDir', action='store', default='.', dest='id_outdir',
	help='Output Directory (Default=. (current directory))')
parser_c.add_argument('-fileType', action='store', default='sam', dest='id_ali_format',
	help='Alignment Format: sam/bl8/gnu-sam (Default: sam)')
parser_c.add_argument('-alignFile', action='store', dest='id_ali_file', required=True,
	help='Alignment file path')

# create the parser for the "REPORT" command
parser_d = subparsers.add_parser('REP', help='Pathoscope Report Module')
parser_d.add_argument('-samtoolsHome', default=None, action='store', 
	dest='rep_samtoolshome', required=False, 
	help='Full Path to samtools binary directory (Default: Uses samtools in system path)')
parser_d.add_argument('-dbhost', action='store', dest='rep_dbhost', required=False, 
	default='localhost', help='specify hostname running mysql if you want to use mysql '
	'instead of hash method in mapping gi to taxonomy id')
parser_d.add_argument('-dbport', action='store', dest='rep_dbport', required=False, 
	default=3306, type=int, help='provide mysql server port if different from default (3306)')
parser_d.add_argument('-dbuser', action='store', dest='rep_dbuser', required=False, default='X', 
	help='user name to access mysql')
parser_d.add_argument('-dbpasswd', action='store', dest='rep_dbpasswd', required=False, default='X', 
	help='provide password associate with user')
parser_d.add_argument('-db', action='store', dest='rep_db', required=False, default='pathodb', 
	help='mysql pathoscope database name (default: pathodb)')
parser_d.add_argument('-outDir', action='store', default='.', dest='rep_outdir',
					help='Output Directory')
parser_d.add_argument('--contig', action='store_true', default=False, dest='rep_contig_flag',
	help='Generate Contig Information (Needs samtools package installed)')
parser_d.add_argument('-samfile', action='store', dest='rep_ali_file', required=True,
					help='SAM Alignment file path')
parser_d.add_argument('--noDisplayCutoff', action='store_true', default=False, 
	dest='rep_nocutoff', help='Do not cutoff display of genomes, even if it is insignificant')

parser_e = subparsers.add_parser('QC', help="PathoQC: Quality control program for high throughput sequencing reads")
parser_e.add_argument('-1', action='store', dest='qc_read1', required=True, default='', help='1st pair of read in PER or SER')
parser_e.add_argument('-2', action='store', dest='qc_read2', required=False, default='', help='2nd pair of read in PER')
parser_e.add_argument('-r', action='store', dest='qc_readL', required=False, default=0, type = int, help='let us know a mean read length (0:ignore. [1]:I want to know the range of read length after trimming, i:user_specified_mean_read_length)')
parser_e.add_argument('-t', action='store', dest='qc_phred_offset', required=False, default=33, type = int, help='specify a phred offset used in encoding base quality(0:not sure?, [33]:phred+33, 64:phred+64)')
parser_e.add_argument('-o', action='store', dest='qc_outdir', required=False, default='', help='specify output directory in full path')
parser_e.add_argument('-s', action='store', dest='qc_sequencer', required=False, default='Illumina', help='tell us which sequencer generates the reads ([Illumina], PacBio, Roche454, IonTorrent)')
parser_e.add_argument('-m', action='store', dest='qc_len_cutoff', required=False, type=int, default=35, help='specify the min read length cutoff[35]')
parser_e.add_argument('-a', action='store', dest='qc_adaptor', required=False, default='N', help='specify an adaptor (Y:have pathoQC detect it, [N]:ignore, ACGT:adaptor)')
parser_e.add_argument('-a2', action='store', dest='qc_adaptor2', required=False, default='N', help='specify a second adaptor if it is different from the first one (Y:have pathoQC detect it, [N]:ignore, ACGT:adaptor)')
#parser_e.add_argument('-k', action='store', dest='qc_primer', required=False, default='N', help='specify a primer')
parser_e.add_argument('-q', action='store', dest='qc_qual_cutoff', required=False, type=int, default=0, help='specify a cutoff of base quality value to trim at the end of reads([0]:ignore, 1:let pathoQC take care of it, i:user_specified_cutoff_value)')
parser_e.add_argument('-R', action='store', dest='qc_lower454', required=False, type=int, default=0, help='set to 1 if you want to mask all lower case bases that may correspond to sequence tag or adaptor in Roche454 or IonTorrent')
parser_e.add_argument('-e', action='store', dest='qc_coff_entropy', required=False, type=int, default=30, help='specify a cutoff of entropy of low sequence complexity to filter out[0..100],default:30, set 0 to disable')
parser_e.add_argument('-d', action='store', dest='qc_derep', required=False, type=int, default=1, help='filtering duplicates: [1] (exact duplicate), 2 (5\' duplicate), 3 (3\' duplicate), 4 (reverse complement exact duplicate), 5 (reverse complement 5\'/3\' duplicate)');
parser_e.add_argument('-g', action='store', dest='qc_add_valid_single', required=False, type=int, default=0, help='Set to 1 if you want to add a valid single read into a reduced valid PER set. Note that this option works only with PER')
parser_e.add_argument('--no_prinseq', action='store_const', dest='qc_on_prinseq', required=False, const=False, default=True, help='to force to skip prinseq module')

parser_e.add_argument('-p', action='store', dest='qc_num_threads', required=False, type=int, default=1, help='specify a total number of cpus to use[1]')
parser_e.add_argument('--debug', action='store_const', dest='qc_debugF', required=False, const=True, default=False, help='working on debug mode')
parser_e.add_argument('--version', action='version', version='PathoQC 1.0')

def main():
	# parse some argument lists
	inputArgs = parser.parse_args()
	
	#### PathoID modules ####
	
	start = time();
	
	if (inputArgs.subcommand=='LIB'):
		################################################$
		#append taxon id in the front of sequence header
		################################################$
		NAs = 'X'
		if inputArgs.lib_dbuser!=NAs and inputArgs.lib_dbpasswd==NAs:
			print 'if you want to use mysql, make sure that you install pathoDB and '
			'also specify the corresponding mysql password correctly '
			'(Ask your mysql admin to access the database).'
		MysqlConf=(inputArgs.lib_dbhost,inputArgs.lib_dbport,inputArgs.lib_dbuser,inputArgs.lib_dbpasswd,inputArgs.lib_db)
		taxon_ids=pathoLib.parse_input_app_build_nt_tgt(inputArgs.lib_taxon_ids)
		exclude_taxon_ids=pathoLib.parse_input_app_build_nt_tgt(inputArgs.lib_exclude_taxon_ids)
		(ncbiNt_ti,ncbiNt_invalid) = pathoLib.append_ti_into_fasta_app(inputArgs.lib_reference,
			taxon_ids, exclude_taxon_ids, inputArgs.lib_subtax,MysqlConf, 
			not(inputArgs.lib_nodesc), inputArgs.lib_online_search, inputArgs.lib_outprefix, 
			inputArgs.lib_outdir)
	
	if (inputArgs.subcommand=='MAP'):
		pathoMapOptions = PathoMapA.PathoMapOptions()
		pathoMapOptions.verbose = inputArgs.verbose
		pathoMapOptions.outDir = inputArgs.map_outdir
		pathoMapOptions.indexDir = inputArgs.map_indexdir
		pathoMapOptions.outAlignFile = inputArgs.map_outalign
		pathoMapOptions.inReadFile = inputArgs.map_inputread
		pathoMapOptions.inReadFilePair1 = inputArgs.map_inputread1
		pathoMapOptions.inReadFilePair2 = inputArgs.map_inputread2
		pathoMapOptions.targetAlignParameters = inputArgs.map_targetalignparams
		pathoMapOptions.filterAlignParameters = inputArgs.map_filteralignparams
		if (len(inputArgs.map_targetref)>0):
			pathoMapOptions.targetRefFiles = inputArgs.map_targetref.split(",")
		if (len(inputArgs.map_filterref)>0):
			pathoMapOptions.filterRefFiles = inputArgs.map_filterref.split(",")
		if (len(inputArgs.map_targetindex)>0):
			pathoMapOptions.targetIndexPrefixes = inputArgs.map_targetindex.split(",")
		if (len(inputArgs.map_filterindex)>0):
			pathoMapOptions.filterIndexPrefixes = inputArgs.map_filterindex.split(",")
		if (len(inputArgs.map_targetalign)>0):
			pathoMapOptions.targetAlignFiles = inputArgs.map_targetalign.split(",")
		if (len(inputArgs.map_filteralign)>0):
			pathoMapOptions.filterAlignFiles = inputArgs.map_filteralign.split(",")
		pathoMapOptions.btHome = inputArgs.map_bthome
		pathoMapOptions.numThreads = inputArgs.map_numthreads
		pathoMapOptions.exp_tag = inputArgs.map_exp_tag + "-"
		PathoMapA.processPathoMap(pathoMapOptions)
	
	if (inputArgs.subcommand=='ID'):
		pathoIdOptions = PathoID.PathoIdOptions(inputArgs.id_ali_file)
		pathoIdOptions.ali_format = inputArgs.id_ali_format
		pathoIdOptions.verbose = inputArgs.verbose
		pathoIdOptions.out_matrix_flag = inputArgs.id_out_matrix
		pathoIdOptions.score_cutoff = inputArgs.id_score_cutoff
		pathoIdOptions.exp_tag = inputArgs.id_exp_tag
		pathoIdOptions.outdir = inputArgs.id_outdir
		pathoIdOptions.emEpsilon = inputArgs.id_emEpsilon
		pathoIdOptions.maxIter = inputArgs.id_maxIter
		pathoIdOptions.piPrior = inputArgs.id_piPrior
		pathoIdOptions.thetaPrior = inputArgs.id_thetaPrior
		pathoIdOptions.noalign = inputArgs.id_noalign
		pathoIdOptions.noCutOff = inputArgs.id_nocutoff
		PathoID.pathoscope_reassign(pathoIdOptions)
	
	if (inputArgs.subcommand=='REP'):
		pathoReportOptions = PathoReportA.PathoReportOptions(inputArgs.rep_ali_file)
		pathoReportOptions.verbose = inputArgs.verbose
		pathoReportOptions.contigFlag = inputArgs.rep_contig_flag
		pathoReportOptions.outDir = inputArgs.rep_outdir
		pathoReportOptions.samtoolsHome = inputArgs.rep_samtoolshome
		pathoReportOptions.noCutOff = inputArgs.rep_nocutoff
		mysqlConf=(inputArgs.rep_dbhost,inputArgs.rep_dbport,inputArgs.rep_dbuser,
			inputArgs.rep_dbpasswd,inputArgs.rep_db)
		pathoReportOptions.mysqlConf = mysqlConf
		PathoReportA.processPathoReport(pathoReportOptions)
	
	if (inputArgs.subcommand=='QC'):
		qcargs = sys.argv[2:]
		pathoqcdir = pathoscopedir + os.path.sep + 'pathoscope' + os.path.sep + 'pathoqc'
		pathoqcfile = pathoqcdir + os.path.sep + 'pathoqc.py'
		if os.path.exists(pathoqcfile):
			cmd = sys.executable
			cmd += " " + pathoqcfile + " "
			cmd += " ".join(qcargs)
			print(cmd)
			os.system(cmd)
		else:
			print("PathoQC (" + pathoqcfile + ") not found. Please download pathoqc_vXXX.tar.gz and "
			"install it ("+pathoqcdir+") from http://sourceforge.net/projects/pathoscope/")
	
	elapsed = time() - start;
	if inputArgs.verbose:
		print "Total Elapsed Time: %d" % (elapsed)

if __name__ == "__main__":
	main()
