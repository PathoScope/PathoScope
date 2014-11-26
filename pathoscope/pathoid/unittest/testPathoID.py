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
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0,pathoscopedir) 

from pathoscope.pathoid import PathoID
import unittest

class TestPathoscopeFunctions(unittest.TestCase):

	def setUp(self):
		self.maxIter = 10
		self.emEpsilon = 1e-7
		self.scoreCutoff = 0.01
		self.verbose = False
		self.piPrior = 0
		self.thetaPrior = 0

	## ex1.bl8 ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	def test_conv_align2GRmat_bl8_1(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex1.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		score = 7036818.21 #Expected score from the alignment file
		expectedU = {0: [0, score], 1: [0, score]} 
		self.assertEquals(len(expectedU), len(U), "Failed bl8 Example 1 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed bl8 Example 1 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed bl8 Example 1 Unique Read %d score Assertion" %read)
		expectedNU = {2: [[1, 2, 0], [score, score, score], [0.3333, 0.3333, 0.3333], score], 
					3: [[1, 2, 0], [score, score, score], [0.3333, 0.3333, 0.3333], score], 
					4: [[1, 2, 0], [score, score, score], [0.3333, 0.3333, 0.3333], score]}
		self.assertEquals(len(expectedNU), len(NU), "Failed bl8 Example 1 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed bl8 Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 1 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed bl8 Example 1 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome3', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 1 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_123', 'read4_123', 'read5_123']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 1 Reads Assertion")

	## ex2.bl8 ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	def test_conv_align2GRmat_bl8_2(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex2.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		score = 7036818.21 #Expected score from the alignment file
		expectedU = {}
		self.assertEquals(expectedU, U, "Failed bl8 Example 2 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [score, score], [0.5, 0.5], score], 
			1: [[0, 1], [score, score], [0.5, 0.5], score], 
			2: [[2, 1], [score, score], [0.5, 0.5], score], 
			3: [[2, 1], [score, score], [0.5, 0.5], score], 
			4: [[2, 0, 1], [score, score, score], [0.3333, 0.3333, 0.3333], score], 
			5: [[2, 0, 1], [score, score, score], [0.3333, 0.3333, 0.3333], score], 
			6: [[2, 0, 1], [score, score, score], [0.3333, 0.3333, 0.3333], score]}
		self.assertEquals(len(expectedNU), len(NU), "Failed bl8 Example 1 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed bl8 Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 2 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed bl8 Example 1 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome2', 'genome1', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 2 Genomes Assertion")
		expectedReads = ['read1_12', 'read2_12', 'read3_13', 'read4_13', 'read5_123', 'read6_123', 'read7_123']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 2 Reads Assertion")

	## ex3.bl8 ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	def test_conv_align2GRmat_bl8_3(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex3.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		score = 7036818.21 #Expected score from the alignment file
		expectedU = {0: [0, score], 1: [0, score], 2: [0, score], 
			3: [1, score], 4: [1, score]} 
		self.assertEquals(len(expectedU), len(U), "Failed bl8 Example 3 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed bl8 Example 3 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed bl8 Example 3 Unique Read %d score Assertion" %read)
		expectedNU = {5: [[2, 1, 0], [score, score, score], [0.3333, 0.3333, 0.3333], score], 
			6: [[2, 1, 0], [score, score, score], [0.3333, 0.3333, 0.3333], score], 
			7: [[2, 1, 0], [score, score, score], [0.3333, 0.3333, 0.3333], score]}
		self.assertEquals(len(expectedNU), len(NU), "Failed bl8 Example 3 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed bl8 Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 3 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed bl8 Example 3 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 3 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_1', 'read4_3', 'read5_2', 'read6_123', 'read7_123', 'read8_123']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 3 Reads Assertion")

	## ex4.bl8 ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	def test_conv_align2GRmat_bl8_4(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex4.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		score = 7036818.21 #Expected score from the alignment file
		expectedU = {2: [0, score]}
		self.assertEquals(len(expectedU), len(U), "Failed bl8 Example 4 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed bl8 Example 4 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed bl8 Example 4 Unique Read %d score Assertion" %read)
		expectedNU = {0: [[0, 1], [score, score], [0.5, 0.5], score], 
			1: [[0, 1], [score, score], [0.5, 0.5], score], 
			3: [[2, 3], [score, score], [0.5, 0.5], score], 
			4: [[1, 3], [score, score], [0.5, 0.5], score]}
		self.assertEquals(len(expectedNU), len(NU), "Failed bl8 Example 4 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed bl8 Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 4 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed bl8 Example 1 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome4', 'genome3', 'genome2', 'genome1']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 4 Genomes Assertion")
		expectedReads = ['read1_34', 'read2_34', 'read3_4', 'read4_12', 'read5_13']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 4 Reads Assertion")

	## ex1.g.sam ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	def test_conv_align2GRmat_gsam_1(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex1.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: [0, 1.0], 1: [0, 1.0]} 
		self.assertEquals(len(expectedU), len(U), "Failed gnusam Example 1 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed gnusam Example 1 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed gnusam Example 1 Unique Read %d score Assertion" %read)
		expectedNU = {2: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33], 
					3: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33], 
					4: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33]}
		self.assertEquals(len(expectedNU), len(NU), "Failed gnusam Example 1 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed gnusam Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 1 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed gnusam Example 1 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 1 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_123', 'read4_123', 'read5_123']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 1 Reads Assertion")

	## ex2.g.sam ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	def test_conv_align2GRmat_gsam_2(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex2.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {}
		self.assertEquals(expectedU, U, "Failed gnusam Example 2 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [0.5, 0.5], [0.5, 0.5], 0.5], 1: [[0, 1], [0.5, 0.5], [0.5, 0.5], 0.5], 
			2: [[0, 2], [0.5, 0.5], [0.5, 0.5], 0.5], 3: [[0, 2], [0.5, 0.5], [0.5, 0.5], 0.5], 
			4: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33], 
			5: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33], 
			6: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33]}
		self.assertEquals(len(expectedNU), len(NU), "Failed gnusam Example 2 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed gnusam Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 2 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed gnusam Example 2 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 2 Genomes Assertion")
		expectedReads = ['read1_12', 'read2_12', 'read3_13', 'read4_13', 'read5_123', 'read6_123', 'read7_123']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 2 Reads Assertion")

	## ex3.g.sam ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	def test_conv_align2GRmat_gsam_3(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex3.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: [0, 1.0], 1: [0, 1.0], 2: [0, 1.0], 3: [1, 1.0], 4: [1, 1.0]} 
		self.assertEquals(len(expectedU), len(U), "Failed gnusam Example 3 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed gnusam Example 3 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed gnusam Example 3 Unique Read %d score Assertion" %read)
		expectedNU = {5: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33], 
			6: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33], 
			7: [[0, 1, 2], [0.33, 0.33, 0.33], [0.3333, 0.3333, 0.3333], 0.33]}
		self.assertEquals(len(expectedNU), len(NU), "Failed gnusam Example 3 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed gnusam Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 3 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed gnusam Example 3 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 3 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_1', 'read4_3', 'read5_2', 'read6_123', 'read7_123', 'read8_123']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 3 Reads Assertion")

	## ex4.g.sam ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	def test_conv_align2GRmat_gsam_4(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex4.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {2: [1, 0.99]} 
		self.assertEquals(len(expectedU), len(U), "Failed gnusam Example 4 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed gnusam Example 4 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed gnusam Example 4 Unique Read %d score Assertion" %read)
		expectedNU = {0: [[0, 1], [0.5, 0.5], [0.5, 0.5], 0.5], 1: [[0, 1], [0.5, 0.5], [0.5, 0.5], 0.5], 
			3: [[2, 3], [0.5, 0.5], [0.5, 0.5], 0.5], 4: [[2, 0], [0.5, 0.5], [0.5, 0.5], 0.5]}
		self.assertEquals(len(expectedNU), len(NU), "Failed gnusam Example 4 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed gnusam Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 4 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed gnusam Example 4 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome3', 'genome4', 'genome1', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 4 Genomes Assertion")
		expectedReads = ['read1_34', 'read2_34', 'read3_4', 'read4_12', 'read5_13']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 4 Reads Assertion")

	## ex1.sam ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	def test_conv_align2GRmat_sam_1(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex1.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		expectedU = {0: [0, scaledScore], 1: [0, scaledScore]}
		self.assertEquals(len(expectedU), len(U), "Failed sam Example 1 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed sam Example 1 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed sam Example 1 Unique Read %d score Assertion" %read)
		expectedNU = {2: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
					3: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
					4: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore]}
		self.assertEquals(len(expectedNU), len(NU), "Failed sam Example 1 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed sam Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 1 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed sam Example 1 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome3', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 1 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_123', 'read4_123', 'read5_123']
		self.assertEquals(expectedReads, reads, "Failed sam Example 1 Reads Assertion")

	## ex2.sam ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	def test_conv_align2GRmat_sam_2(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex2.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		expectedU = {}
		self.assertEquals(expectedU, U, "Failed sam Example 2 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			1: [[0, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			2: [[2, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			3: [[2, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			4: [[0, 2, 1], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
			5: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
			6: [[0, 1, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore]}
		self.assertEquals(len(expectedNU), len(NU), "Failed sam Example 2 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed sam Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 2 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed sam Example 2 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome2', 'genome1', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 2 Genomes Assertion")
		expectedReads = ['read1_12', 'read2_12', 'read3_13', 'read4_13', 'read5_123', 'read6_123', 'read7_123']
		self.assertEquals(expectedReads, reads, "Failed sam Example 2 Reads Assertion")

	## ex3.sam ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	def test_conv_align2GRmat_sam_3(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex3.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		expectedU = {0: [0, scaledScore], 1: [0, scaledScore], 2: [0, scaledScore], 3: [1, scaledScore], 4: [1, scaledScore]}
		self.assertEquals(len(expectedU), len(U), "Failed sam Example 3 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed sam Example 3 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed sam Example 3 Unique Read %d score Assertion" %read)
		expectedNU = {5: [[2, 0, 1], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
					6: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
					7: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore]}
		self.assertEquals(len(expectedNU), len(NU), "Failed sam Example 3 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed sam Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 3 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed sam Example 3 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 3 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_1', 'read4_3', 'read5_2', 'read6_123', 'read7_123', 'read8_123']
		self.assertEquals(expectedReads, reads, "Failed sam Example 3 Reads Assertion")

	## ex4.sam ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	def test_conv_align2GRmat_sam_4(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex4.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		expectedU = {2: [0, scaledScore]}
		self.assertEquals(len(expectedU), len(U), "Failed sam Example 4 Unique Reads length Assertion")
		for read in expectedU:
			self.assertEquals(expectedU[read][0], U[read][0], 
				"Failed sam Example 4 Unique Read %d genome mapping Assertion" %read)
			self.assertAlmostEquals(expectedU[read][1], U[read][1], 2, 
				"Failed sam Example 4 Unique Read %d score Assertion" %read)
		expectedNU = {0: [[0, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			1: [[1, 0], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			3: [[2, 3], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			4: [[2, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore]}
		self.assertEquals(len(expectedNU), len(NU), "Failed sam Example 4 Non-Unique Reads length Assertion")
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			for j in range(len(expectedNU[read][1])):
				self.assertAlmostEquals(expectedNU[read][1][j], NU[read][1][j], 2, 
					"Failed sam Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 4 Non-Unique Read %d proportion Assertion" %read)
			self.assertAlmostEquals(expectedNU[read][3], NU[read][3], 2, 
				"Failed sam Example 4 Non-Unique Read %d weight Assertion" %read)
		expectedGenomes = ['genome4', 'genome3', 'genome1', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 4 Genomes Assertion")
		expectedReads = ['read1_34', 'read2_34', 'read3_4', 'read4_12', 'read5_13']
		self.assertEquals(expectedReads, reads, "Failed sam Example 4 Reads Assertion")

	## Example 1 ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	# Answer: all to genome1
	def test_pathoscope_em_1(self):
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		U = {0: [0, scaledScore], 1: [0, scaledScore]}
		NU = {2: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
					3: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
					4: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore]}
		### Genome hash
		genomes = {0:"genome1", 1:"genome3", 2:"genome2"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, 
			self.maxIter, self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.6, 0.2, 0.2]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 1 Initial PI Assertion")
		expectedPi = [1.0, 0.0, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 4, "Failed EM Example 1 PI Assertion")
		expectedTheta = [1.0, 0.0, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 4, "Failed EM Example 1 Theta Assertion")
		expectedNU = {2: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0, 1.0, 0], scaledScore], 
			3: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0, 1.0, 0], scaledScore], 
			4: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0, 0, 1.0], scaledScore]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed EM Example 1 Non-Unique Read %d proportion Assertion" %read)

	## Example 2 ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	# Answer: all to genome1
	def test_pathoscope_em_2(self):
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		U = {}
		NU = {0: [[0, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			1: [[0, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			2: [[2, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			3: [[2, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			4: [[0, 2, 1], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
			5: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
			6: [[0, 1, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore]}
		### Genome hash
		genomes = {0:"genome2", 1:"genome1", 2:"genome3"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, 
			self.maxIter, self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.2857, 0.4286, 0.2857]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 2 Initial PI Assertion j=%d" %j)
		expectedPi = [0.0, 1.0, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 4, "Failed EM Example 2 PI Assertion")
		expectedTheta = [0.0, 1.0, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 4, "Failed EM Example 2 Theta Assertion")
		expectedNU = {0: [[0, 1], [scaledScore, scaledScore], [0.0, 1.0], scaledScore], 
			1: [[0, 1], [scaledScore, scaledScore], [0.0, 1.0], scaledScore], 
			2: [[2, 1], [scaledScore, scaledScore], [0.0, 1.0], scaledScore], 
			3: [[2, 1], [scaledScore, scaledScore], [0.0, 1.0], scaledScore], 
			4: [[0, 2, 1], [scaledScore, scaledScore, scaledScore], [0.0, 0.0, 1.0], scaledScore], 
			5: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.0, 1.0, 0.0], scaledScore], 
			6: [[0, 1, 2], [scaledScore, scaledScore, scaledScore], [0.0, 1.0, 0.0], scaledScore]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed EM Example 2 Non-Unique Read %d proportion Assertion" %read)

	## Example 3 ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	# Answer: 6 to genome1; 2 to genome2
	def test_pathoscope_em_3(self):
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		U = {0: [0, scaledScore], 1: [0, scaledScore], 2: [0, scaledScore], 
			3: [1, scaledScore], 4: [1, scaledScore]}
		NU = {5: [[2, 0, 1], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
			6: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore], 
			7: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.3333, 0.3333, 0.3333], scaledScore]}
		### Genome hash
		genomes = {0:"genome1", 1:"genome2", 2:"genome3"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, 
			self.maxIter, self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.5, 0.375, 0.125]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 3 Initial PI Assertion j=%d" %j)
		expectedPi = [0.75, 0.25, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 2, "Failed EM Example 3 PI Assertion")
		expectedTheta = [1.0, 0.0, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 2, "Failed EM Example 3 Theta Assertion")
		expectedNU = {5: [[2, 0, 1], [scaledScore, scaledScore, scaledScore], [0.0, 1.0, 0.0], scaledScore], 
			6: [[1, 0, 2], [scaledScore, scaledScore, scaledScore], [0.0, 1.0, 0.0], scaledScore], 
			7: [[2, 1, 0], [scaledScore, scaledScore, scaledScore], [0.0, 0.0, 1.0], scaledScore]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 2, 
					"Failed EM Example 3 Non-Unique Read %d proportion Assertion" %read)

	## Example 4 ##
	# 2 reads to 3,4; 1 read to 4; 1 read to 1,2; 1 read to 1,3 
	# Answer: 3 reads to genome4; 2 reads to genome1
	def test_pathoscope_em_4(self):
		scaledScore = 2.6881171418161356e+43 #Expected Re-scaled score
		U = {2: [0, scaledScore]}
		NU = {0: [[0, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			1: [[1, 0], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			3: [[2, 3], [scaledScore, scaledScore], [0.5, 0.5], scaledScore], 
			4: [[2, 1], [scaledScore, scaledScore], [0.5, 0.5], scaledScore]}
		### Genome hash
		genomes = {0:"genome4", 1:"genome3", 2:"genome1", 3:"genome2"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, 
			self.maxIter, self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.4, 0.3, 0.2, 0.1]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 4 Initial PI Assertion j=%d" %j)
		expectedPi = [0.20, 0.57, 0.23, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 2, "Failed EM Example 4 PI Assertion")
		expectedTheta = [0.0, 0.72, 0.28, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 2, "Failed EM Example 4 Theta Assertion")
		expectedNU = {0: [[0, 1], [scaledScore, scaledScore], [0.0, 1.0], scaledScore], 
			1: [[1, 0], [scaledScore, scaledScore], [1.0, 0.0], scaledScore], 
			3: [[2, 3], [scaledScore, scaledScore], [1.0, 0.0], scaledScore], 
			4: [[2, 1], [scaledScore, scaledScore], [0.14, 0.86], scaledScore]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 2, 
					"Failed EM Example 4 Non-Unique Read %d proportion Assertion" %read)

	## Example to test the trivial case with only unique reads: no non-unique reads ##
	def test_pathoscope_em_x1(self):
	### 2 unique reads: 2 reads to genome1
		U = {0: [0, 0.5], 1: [0, 0.5]}
	### 0 non-unique reads: 0 total reads  readnum:[[genomes],[qij],[xij]]
		NU = {}
	### Genome hash
		genomes = {0:"genome1", 1:"genome2", 2:"genome3", 3:"genome4"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, 
			self.maxIter, self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [1.0, 0.0, 0.0, 0.0]
		expectedPi = [1.0, 0.0, 0.0, 0.0]
		expectedTheta = [0.0, 0.0, 0.0, 0.0]
		diffInitPi = [abs(a - b) for a, b in zip(initPi, expectedInitPi)]
		diffPi = [abs(a - b) for a, b in zip(pi, expectedPi)]
		diffTheta = [abs(a - b) for a, b in zip(theta, expectedTheta)]
		delta = [0.0001]*4
		self.assertTrue(diffInitPi < delta, "Failed initPi Assertion")
		self.assertTrue(diffPi < delta, "Failed Pi Assertion")
		self.assertTrue(diffTheta < delta, "Failed Theta Assertion")


if __name__ == '__main__':
	unittest.main()

