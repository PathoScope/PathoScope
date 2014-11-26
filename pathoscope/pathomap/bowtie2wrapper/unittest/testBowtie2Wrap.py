#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Unit Test Functions to check the Pathoscope EM algorithm functions

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

import os,sys
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
	os.path.dirname(os.path.abspath(__file__))))))
sys.path.insert(0,pathoscopedir) 

from pathoscope.pathomap.bowtie2wrapper import bowtie2Wrap
from pathoscope.utils import seqParse
import unittest

class TestBowtie2WrapFunctions(unittest.TestCase):

	def setUp(self):
		self.targetAlignFile = "target.sam"
		self.filterAlignFiles = ["filter1.sam", "filter2.sam"]
		self.outAlignFile = "pathomap_filtered.sam"
		self.fastqOutFile = "target_aligned.fq"
		self.verbose = False

	def test_filter_alignment(self):
		bowtie2Wrap.filter_alignment(self.targetAlignFile, self.filterAlignFiles, 
			self.outAlignFile)
		expectedReadId = ["HWI-ST998R:270:H7NJ9ADXX:1:1101:1797:1927:A"]
		with open(self.outAlignFile,'r') as in1:
			count = 0
			for ln in in1:
				if (ln[0] == '@' or ln[0] == '#'):
					continue
				l = ln.split('\t')
				readId=l[0]
				self.assertTrue(count < len(expectedReadId), 
					"Filter Alignment: Expected number of reads Mismatch!")
				self.assertTrue(readId == expectedReadId[count], 
					"Filter Alignment: Expected Reads Mismatch!")
				count += 1
			self.assertTrue(count == len(expectedReadId), 
				"Filter Alignment: Expected number of reads Mismatch!")
					

	def test_extractRead(self):
		bowtie2Wrap.extractRead(self.targetAlignFile, self.fastqOutFile)
		expectedReadId = ["HWI-ST998R:270:H7NJ9ADXX:1:1101:1797:1927", 
			"HWI-ST998R:270:H7NJ9ADXX:1:1101:1797:1927:A",
			"HWI-ST998R:270:H7NJ9ADXX:1:1101:1797:1927:B"]
		with open(self.fastqOutFile,'r') as fp:
			count = 0
			for r in seqParse.parse(fp,'fastq'):
				self.assertTrue(count < len(expectedReadId), 
					"Extract Read: Expected number of reads Mismatch!")
				self.assertTrue(r.id == expectedReadId[count], 
					"Extract Read: Expected Reads Mismatch!")
				count += 1
			self.assertTrue(count == len(expectedReadId), 
				"Extract Read: Expected number of reads Mismatch!")

if __name__ == '__main__':
	unittest.main()

