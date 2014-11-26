#!/usr/bin/python
# Author: Solaiappan Manimaran
# Utility class and functions to parse fasta/fastq sequence files

#	Pathoscope2
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
#

def parse(sfile, fmt):
	if fmt=='fastq':
		return fastq_parse(sfile)
	else:
		return fasta_parse(sfile)

def fasta_parse(sfile):
	while True:
		line = sfile.readline()
		if line == "":
			return  # Premature end of file, or just empty?
		if line[0] == ">":
			break

	while True:
		if line[0] != ">":
			raise ValueError(
				"Records in Fasta files should start with '>' character")
		idline = line[1:].rstrip()
		lines = []
		line = sfile.readline()
		while True:
			if not line:
				break
			if line[0] == ">":
				break
			lines.append(line.rstrip())
			line = sfile.readline()
		#Remove trailing whitespace, and any internal spaces
		#(and any embedded \r which are possible in mangled files
		#when not opened in universal read lines mode)
		seqid = idline.split(None, 1)[0]
		sequence = Sequence(seqid, idline)
		sequence.seq = "".join(lines).replace(" ", "").replace("\r", "")
		yield sequence
		if not line:
			return  # StopIteration
	assert False, "Should not reach this line"

def fastq_parse(sfile):
	while True:
		line = sfile.readline()
		if line == "":
			return  # Premature end of file, or just empty?
		if line[0] == "@":
			break

	while True:
		qualLine = False
		if line[0] != "@":
			raise ValueError(
				"Records in Fastq files should start with '@' character")
		idline = line[1:].rstrip()
		seqLines = []
		sLen = 0
		qualLines = []
		qLen = 0
		line = sfile.readline()
		while True:
			if not line:
				break
			if line[0] == "@":
				if (qLen >= sLen): # Checking whether quality line has exceeded sequence line
					break
			if line[0] == "+":
				qualLine = True
				line = sfile.readline()
				continue
			strippedLine = line.rstrip().replace(" ", "").replace("\r", "")
			if qualLine == True:
				qLen += len(strippedLine)
				qualLines.append(strippedLine)
			else:
				sLen += len(strippedLine)
				seqLines.append(strippedLine)
			line = sfile.readline()
		seqid = idline.split(None, 1)[0]
		sequence = Sequence(seqid, idline)
		#Remove trailing whitespace, and any internal spaces
		#(and any embedded \r which are possible in mangled files
		#when not opened in universal read lines mode)
		sequence.seq = "".join(seqLines)
		sequence.qual = "".join(qualLines)
		yield sequence
		if not line:
			return  # StopIteration
	assert False, "Should not reach this line"


class Sequence:
	id = ""
	description = ""
	length = 0
	seq = ""
	qual = None
	def __init__(self, sid, desc):
		self.id = sid
		self.description = desc.strip()
	def format(self, fmt):
		if fmt=='fastq':
			line = '@%s %s\n%s\n+\n%s\n' % (self.id, self.description, self.seq, self.qual)
		else:
			line = '>%s %s\n%s\n' % (self.id, self.description, self.seq)
		return line
