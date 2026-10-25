##############################################################################
#	Copyright (c) 2020-2026 Broad Institute of MIT and Harvard
#		All rights reserved
#    CONFIDENTIAL AND PROPRIETARY INFORMATION
##############################################################################

#Filename:	lsmEnsemblUtilities.py
#Purpose:	Manipulate FASTA formatted sequence databases obtained from Ensembl and Gencode
#Authors:	Karl Clauser 
#Created:	June 2020

import gzip
#	KRC 6/27/2021 software engineering thoughts
#	the look-up dictionaries
#	DDtranscriptInfo (created from the Gencode GTF file)
#	DDproteinID_info (created from the BioMart .txt file)
#	are read into objects of my PERL style. A more Pythonic way of doing it would be
#	to read into a DataFrame via Pandas or an ordered dictionary via CSV from an 
#	Encapsulation standpoint those look-up dictionaries should NOT be passed out of
#	this module but remain internal, to mitigate the external impact on subsequent
#	internal changes in the reading and objectifying of these look-ups

##############################################################################
############################ parseGencodeGTFfile #############################
##############################################################################
##A GTF file is a hybrid tab separated variable format containing exactly 9 columns
	##The 9th column consists of semi-colon space delimited key-value pairs, where the value is in quotes
		
##top few lines of gencode.v34.basic.annotation.gtf
#	#description: evidence-based annotation of the human genome (GRCh38), version 34 (Ensembl 100)
#	##provider: GENCODE
#	##contact: gencode-help@ebi.ac.uk
#	##format: gtf
#	##date: 2020-03-24
#	chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.1";
#	chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.1"; havana_transcript "OTTHUMT00000362751.1";
#	chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.1"; havana_transcript "OTTHUMT00000362751.1";
#	chr2	HAVANA	transcript	29192774	29223900	.	-	.	gene_id "ENSG00000171094.18"; transcript_id "ENST00000642122.1"; gene_type "protein_coding"; gene_name "ALK"; transcript_type "protein_coding"; transcript_name "ALK-207"; level 2; protein_id "ENSP00000493203.1"; hgnc_id "HGNC:427"; tag "basic"; tag "CCDS"; ccdsid "CCDS86828.1"; havana_gene "OTTHUMG00000152034.1"; havana_transcript "OTTHUMT00000493449.1";

##https://www.gencodegenes.org/pages/data_format.html
##col num	content					values/format
##1			chromosome name			chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M} or GRC accession a
##2			annotation source		{ENSEMBL,HAVANA}
##3			feature type			{gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
##4			genomic start location	integer-value (1-based)
##5			genomic end location	integer-value
##6			score					(not used)	.
##7			genomic strand			{+,-}
##8			genomic phase			 (for CDS features)	{0,1,2,.}
##9			key-value pairs	see below
	##key name
	##gene_id
	##transcript_id
	##gene_type
	##gene_status
	##gene_name
	##transcript_type
	##transcript_status
	##transcript_name
	##exon_number
	##exon_id
	#protein_id

# A multi-exon gene will will have multiple lines, different feature type & start end, with same ENSG, ENST, ENSP
#	goal is to get genome coordinates from start_codon to stop_codon - 1
"""
chr15	HAVANA	transcript	48120990	48142672	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48120990	48121165	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 1; exon_id "ENSE00001382203.3"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48121045	48121165	.	+	0	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 1; exon_id "ENSE00001382203.3"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	start_codon	48121045	48121047	.	+	0	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 1; exon_id "ENSE00001382203.3"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48121857	48122036	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 2; exon_id "ENSE00001364100.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48121857	48122036	.	+	2	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 2; exon_id "ENSE00001364100.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48134258	48134341	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 3; exon_id "ENSE00003489229.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48134258	48134341	.	+	2	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 3; exon_id "ENSE00003489229.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48134435	48134538	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 4; exon_id "ENSE00003667994.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48134435	48134538	.	+	2	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 4; exon_id "ENSE00003667994.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48134884	48134984	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 5; exon_id "ENSE00003647742.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48134884	48134984	.	+	0	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 5; exon_id "ENSE00003647742.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48136683	48136963	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 6; exon_id "ENSE00001385395.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48136683	48136963	.	+	1	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 6; exon_id "ENSE00001385395.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48138969	48139175	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 7; exon_id "ENSE00001364865.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48138969	48139175	.	+	2	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 7; exon_id "ENSE00001364865.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48141113	48141214	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 8; exon_id "ENSE00001384487.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48141113	48141214	.	+	2	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 8; exon_id "ENSE00001384487.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	exon	48142029	48142672	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 9; exon_id "ENSE00001919044.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	CDS	48142029	48142348	.	+	2	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 9; exon_id "ENSE00001919044.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	stop_codon	48142349	48142351	.	+	0	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 9; exon_id "ENSE00001919044.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	UTR	48120990	48121044	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 1; exon_id "ENSE00001382203.3"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
chr15	HAVANA	UTR	48142349	48142672	.	+	.	gene_id "ENSG00000188467.11"; transcript_id "ENST00000341459.8"; gene_type "protein_coding"; gene_name "SLC24A5"; transcript_type "protein_coding"; transcript_name "SLC24A5-201"; exon_number 9; exon_id "ENSE00001919044.1"; level 2; protein_id "ENSP00000341550.3"; transcript_support_level "1"; hgnc_id "HGNC:20611"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS10128.1"; havana_gene "OTTHUMG00000131494.5"; havana_transcript "OTTHUMT00000254340.3";
"""
def fParseGencodeGTFfile_CDS(sGencodeAnnotationPath):
	#8/30/2025 revision to get genome coordinates that are for the CDS
	Lheaders = ['chromosome_name', 'annotation_source', 'feature_type', 'genomic_start_location', 'genomic_end_location', 'score', 'genomic_strand', 'genomic_phase', 'key_values']
	LDtranscriptInfo = fparseTSVfileNoHeader(sGencodeAnnotationPath, Lheaders)
	#LDtranscriptInfo = [] #initialize empty list
	
	#create a lookup dictionary with transcript_id as the primary key
		# expand the key_values column to get access to 	gene_id, transcript_id, transcript_type, protein_id
	DDtranscriptInfo = {}	#initialize empty dictionary
	iNumTranscriptIDsParsed = 0
	for DtranscriptInfo in LDtranscriptInfo:
		Lkey_values = DtranscriptInfo['key_values'].split('\"; ') #split on double-quote semi-colon space
		for sKeyValuePair in Lkey_values:
			LkeyValuePair = sKeyValuePair.split(' \"') #split on space double-quote
			#LkeyValuePair[1].rstrip('\"') #trim off trailing "
			if(len(LkeyValuePair) > 1):
				DtranscriptInfo[LkeyValuePair[0]] = LkeyValuePair[1]
			else:
				DtranscriptInfo[LkeyValuePair[0]] = '' #empty string
		if('transcript_id' in DtranscriptInfo):
			iNumTranscriptIDsParsed += 1
			#DDtranscriptInfo[DtranscriptInfo['transcript_id']] = DtranscriptInfo #if duplicates present, this will overwrite with the last occurence in the .gtf file
			#KRC 8/30/2025 no longer use last entry for each transcript for genomic coordinates could be just UTR
			if DtranscriptInfo['transcript_id'] in DDtranscriptInfo:
				#includes genomic coordinates for start & stop codon
				#note that for - strand
				# to put genomic coordinates in standard low to high order
				# 'genomic_start_location_protein' 'genomic_end_location_protein' need to be reversed
				# genomic_coordinates = ['chromosome_name']:['genomic_end_location_protein']-['genomic_start_location_protein']:['genomic_strand']
				if DtranscriptInfo['feature_type'] == 'start_codon':
					if DtranscriptInfo['genomic_strand'] == '+':
						DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_start_location_protein'] = DtranscriptInfo['genomic_start_location']
					elif DtranscriptInfo['genomic_strand'] == '-':
						DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_start_location_protein'] = DtranscriptInfo['genomic_end_location']
				elif DtranscriptInfo['feature_type'] == 'stop_codon':
					#term of protein relative to stop codon
					if DtranscriptInfo['genomic_strand'] == '+':
						DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_end_location_protein'] = str(int(DtranscriptInfo['genomic_start_location']) - 1)
					elif DtranscriptInfo['genomic_strand'] == '-':
						DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_end_location_protein'] = str(int(DtranscriptInfo['genomic_start_location']) + 3)
			elif DtranscriptInfo['feature_type'] == 'transcript':
				#initially fills with 'genomic_start_location', 'genomic_end_location' for full transcript
				DDtranscriptInfo[DtranscriptInfo['transcript_id']] = DtranscriptInfo
				#fill defaults in case start or stop is absent
				if DtranscriptInfo['genomic_strand'] == '+':
					DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_start_location_protein']	= DtranscriptInfo['genomic_start_location']
					DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_end_location_protein']		= DtranscriptInfo['genomic_end_location']
				elif DtranscriptInfo['genomic_strand'] == '-':
					#these will be reversed in fUpdateGencodeHeadersCreateCategoriesFile to restore low-high genomic coordinates 
					DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_start_location_protein']	= DtranscriptInfo['genomic_end_location']
					DDtranscriptInfo[DtranscriptInfo['transcript_id']]['genomic_end_location_protein']		= DtranscriptInfo['genomic_start_location']


	#verify how many transcript_ids parsed 
	iNumItemsStored = len(DDtranscriptInfo)
	print('<br>Number of transcript_ids parsed from ', sGencodeAnnotationPath, ': ', iNumTranscriptIDsParsed)
	print('<br>Number of transcript_ids stored from ', sGencodeAnnotationPath, ': ', iNumItemsStored)
#s
	return (DDtranscriptInfo) #Should keep this as an internal object to this module

##############################################################################
############################ parseGencodeGTFfile #############################
##############################################################################
def fParseGencodeGTFfile(sGencodeAnnotationPath):
	Lheaders = ['chromosome_name', 'annotation_source', 'feature_type', 'genomic_start_location', 'genomic_end_location', 'score', 'genomic_strand', 'genomic_phase', 'key_values']
	LDtranscriptInfo = fparseTSVfileNoHeader(sGencodeAnnotationPath, Lheaders)
	#LDtranscriptInfo = [] #initialize empty list
	
	#create a lookup dictionary with transcript_id as the primary key
		# expand the key_values column to get access to 	gene_id, transcript_id, transcript_type, protein_id
	DDtranscriptInfo = {}	#initialize empty dictionary
	iNumTranscriptIDsParsed = 0
	for DtranscriptInfo in LDtranscriptInfo:
		Lkey_values = DtranscriptInfo['key_values'].split('\"; ') #split on double-quote semi-colon space
		for sKeyValuePair in Lkey_values:
			LkeyValuePair = sKeyValuePair.split(' \"') #split on space double-quote
			#LkeyValuePair[1].rstrip('\"') #trim off trailing "
			if(len(LkeyValuePair) > 1):
				DtranscriptInfo[LkeyValuePair[0]] = LkeyValuePair[1]
			else:
				DtranscriptInfo[LkeyValuePair[0]] = '' #empty string
		if('transcript_id' in DtranscriptInfo):
			iNumTranscriptIDsParsed += 1
			DDtranscriptInfo[DtranscriptInfo['transcript_id']] = DtranscriptInfo #if duplicates present, this will overwrite with the last occurence in the .gtf file

	#verify how many transcript_ids parsed 
	iNumItemsStored = len(DDtranscriptInfo)
	print('<br>Number of transcript_ids parsed from ', sGencodeAnnotationPath, ': ', iNumTranscriptIDsParsed)
	print('<br>Number of transcript_ids stored from ', sGencodeAnnotationPath, ': ', iNumItemsStored)
#s
	return (DDtranscriptInfo) #Should keep this as an internal object to this module

##############################################################################
######################### fParseEnsemblBiomartFile ###########################
##############################################################################
def	fParseEnsemblBiomartFile(sBioMartFile, sProteinIDkeyName = '', sProteinNameKeyName = ''):
	
	LDproteinInfo = fParseTSVfileHeader(sBioMartFile)
	
	#make hash for cross referencing by transcriptID
	DDproteinIDinfo = {}	#initialize empty dictionary
	for DproteinInfo in LDproteinInfo:
		#sProteinIDkeyName 		hg19: 'Ensembl Protein ID'			Ensembl-75_GRCh37_ENST_mart_export.txt
		#sProteinNameKeyName	hg19: 'Description'					Ensembl-75_GRCh37_ENST_mart_export.txt
		#sProteinIDkeyName 		hg38: 'Protein stable ID version'	Ensembl-100_GRCh38_mart_export.txt
		#sProteinNameKeyName	hg38: 'Gene description'			Ensembl-100_GRCh38_mart_export.txt
		if(sProteinIDkeyName in DproteinInfo):
			if(sProteinNameKeyName in  DproteinInfo):
				Lfields = DproteinInfo[sProteinNameKeyName].split(' [Source') #trim off superfluous text
				DproteinInfo['proteinName'] = Lfields [0]
			DDproteinIDinfo[ DproteinInfo[sProteinIDkeyName] ] = DproteinInfo
		
		

	iNumInfoEntries = len(DDproteinIDinfo )
	print('Read in gene symbols, protein names, genome coordinates, transcriptIDs, and geneIDS for ', iNumInfoEntries, ' proteins from ', sBioMartFile)

	return(DDproteinIDinfo)  #Should keep this as an internal object to this module
		
# belongs in lsmTableRead.py
##############################################################################
########################### fParseTSVfileNoHeader ############################
##############################################################################

def fparseTSVfileNoHeader(sFile_in, Lheaders):

	iLineNumber = 0
	if not sFile_in.endswith('.gz'):
		FileIn = open(sFile_in, 'r')
		Llines = FileIn.read().splitlines()	#list of lines with \n removed
	else:
		with gzip.open(sFile_in, 'rt') as FileIn:
			Llines = FileIn.read().splitlines()	#list of lines with \n removed

	LDentries = [] #make empty list
	for sLine in Llines:
		if (sLine == ''):				continue	#skip blank line in file
		if (sLine.startswith('#') ):	continue	#skip comment lines that have a leading #
		#prevents timing out
		#if ( ( $lineNumber % $NUM_LINES_READ_BEFORE_SENDING_FEEDBACK == 0 ) && (printf ( "\t" ) != 1 ) ){ &fatalErrorHTML("<P>Server timed out while reading $TSVfilename.\n", "lsmTableRead.pl") } #exits script if browser connection terminated by user
		Lvalues = sLine.split('\t')
		iLineNumber += 1
		DcurrentValues = {} #make empty dictionary
		for sHeader in Lheaders:
			#$value = shift @values;	# always shift-- preserve positional relation with col keys
			if(len(Lvalues) > 0):
				sValue = Lvalues.pop(0) # maintain register left to right with Lheaders
			else:
				sValue = ''	#initialize with empty string
			DcurrentValues[sHeader] = sValue
		LDentries.append(DcurrentValues)

	FileIn.close()
#s
	return(LDentries);

##############################################################################
############################ parseTSVfileNoHeader ############################
##############################################################################

def fParseTSVfileHeader(sFile_in):

	iLineNumber = 0

	FileIn = open(sFile_in, 'r')
	Llines = FileIn.read().splitlines()	#list of lines with \n removed

	LDentries = [] #make empty list
	Lheaders = [] #make empty list
	for sLine in Llines:
		if (sLine == ''):				continue	#skip blank line in file
		if (sLine.startswith('#') ):	continue	#skip comment lines that have a leading #		
		#prevents timing out
		#if ( ( $lineNumber % $NUM_LINES_READ_BEFORE_SENDING_FEEDBACK == 0 ) && (printf ( "\t" ) != 1 ) ){ &fatalErrorHTML("<P>Server timed out while reading $TSVfilename.\n", "lsmTableRead.pl") } #exits script if browser connection terminated by user
		#chomp($line);
		Lvalues = sLine.split('\t')
		iLineNumber += 1
		if(iLineNumber == 1):			#read in headers
			Lheaders = Lvalues
			continue

		DcurrentValues = {} #make empty dictionary
		for sHeader in Lheaders:
			#$value = shift @values;	# always shift-- preserve positional relation with col keys
			if(len(Lvalues) > 0):
				sValue = Lvalues.pop(0) # maintain register left to right with Lheaders
			else:
				sValue = ''	#initialize with empty string
			DcurrentValues[sHeader] = sValue
		LDentries.append(DcurrentValues)

	FileIn.close()

	return(LDentries);
	
##############################################################################
#############################  cleanEnsemblFasta #############################
##############################################################################
#purpose: to filter out proteins that have undesirable sources i.e not protein coding or not mapped to a standard chromosome

def fCleanEnsemblFasta(LDfproteinsAll):
#e
#>ENSP00000328693 pep:novel chromosome:NCBI35:1:904515:910768:1 gene:ENSG00000158815:transcript:ENST00000328693 gene_biotype:protein_coding transcript_biotype:protein_coding
# ^               ^   ^     ^                                   ^                    ^                          ^                            ^ 
# ID              |   |  LOCATION                          GENE:stable gene ID       |                       GENE: gene biotype           TRANSCRIPT: transcript biotype
#                 | STATUS                                           TRANSCRIPT: stable transcript ID
#               SEQTYPE

	iNumProteinsAll = len(LDfproteinsAll)
	print('<br># Proteins before Ensemble cleaning: ', iNumProteinsAll)
	
	LDfproteinsClean = []  #initialize empty list
	LDfproteinsDirty = []  #initialize empty list

	iNumNotProteinCoding = 0
	iNumNotChromosomeAssigned = 0
	for Dfprotein in LDfproteinsAll:
		#(@AspaceDelimitedFileds) = split (/\s+/, Dfprotein['header'])
		if(Dfprotein['header'].find('transcript_biotype:protein_coding') > 0):
			if(Dfprotein['header'].find('chromosome:GRCh37:GL') > 0 
				or Dfprotein['header'].find('chromosome:GRCh37:H') > 0
				or Dfprotein['header'].find('chromosome:GRCh37:LRG') > 0
				or Dfprotein['header'].find('chromosome:GRCh38:CHR_') > 0
				or Dfprotein['header'].find('scaffold:GRCh38') > 0
				):
				iNumNotChromosomeAssigned += 1
				LDfproteinsDirty.append(Dfprotein)
				continue
			else:
				LDfproteinsClean.append(Dfprotein)
		else:
			iNumNotProteinCoding += 1
			LDfproteinsDirty.append(Dfprotein)

	iNumProteinsClean = len(LDfproteinsClean)
	print('<br># Proteins after Ensemble cleaning: ', iNumProteinsClean)
	print('<br>Number entries removed that were not transcript_biotype:protein_coding:\t', iNumNotProteinCoding)
	print('<br>Number entries removed that were transcript_biotype:protein_coding, but not assigned to a chromosome:\t', iNumNotChromosomeAssigned)

	return(LDfproteinsClean, LDfproteinsDirty)

###############################################################################
##############################  cleanGencodeFasta #############################
###############################################################################
##s
def fCleanGencodeFasta(LDfproteinsAll, DDtranscriptInfo, bv19Filter = False, sSubsetTag = ''): #sSubsetTag: 'basic', 'GENCODE_Primary'

	iNumProteinsAll = len(LDfproteinsAll)
	print('<br># Proteins before Gencode cleaning: ', iNumProteinsAll)

	LDfproteinsClean = []  #initialize empty list
	LDfproteinsDirty = []  #initialize empty list
	LDfproteinsUnannotated = []  #initialize empty list
	Ltranscript_typePC = ['protein_coding', 'IG_C_gene', 'TR_C_gene'] #KRC 4/26/2023 14 IG_C and 6 TR_C, readily found in proteomics datasets

	iNumNotProteinCoding = 0
	iNumNotChromosomeAssigned = 0
	iNumTranscriptTypesFound = 0
	iNumTranscriptIDsFound = 0
	iNumProteinsUnannotated = 0
	iNumTranscriptStatusNotKnown = 0
	iCountPrimary = 0
	icountBasic = 0
	iCountBothBasicPrimary = 0
	for Dfprotein in LDfproteinsAll:
		#if(exists($RHRH_transcriptInfo->{Dfprotein['transcriptID']})){
		if(Dfprotein['transcriptID'] in DDtranscriptInfo):
			iNumTranscriptIDsFound += 1
			#if(exists($RHRH_transcriptInfo->{Dfprotein['transcriptID']}->{'transcript_type'])){
			if('transcript_type' in DDtranscriptInfo[Dfprotein['transcriptID']]):
				iNumTranscriptTypesFound += 1
# 				if(DDtranscriptInfo[Dfprotein['transcriptID']]['transcript_type'] == 'protein_coding'
# 					or DDtranscriptInfo[Dfprotein['transcriptID']]['transcript_type'] == 'IG_C_gene'
# 					or DDtranscriptInfo[Dfprotein['transcriptID']]['transcript_type'] == 'TR_C_gene'):
				if(any(DDtranscriptInfo[Dfprotein['transcriptID']]['transcript_type'] in x  for x in Ltranscript_typePC)):
					if(DDtranscriptInfo[Dfprotein['transcriptID']]['chromosome_name'] == 'KarltoDo' ):
						iNumNotChromosomeAssigned += 1
						LDfproteinsDirty.append(Dfprotein)
						continue
					elif bv19Filter:
						if DDtranscriptInfo[Dfprotein['transcriptID']]['transcript_status'] == 'KNOWN' or DDtranscriptInfo[Dfprotein['transcriptID']]['tag'] == 'exp_conf' or DDtranscriptInfo[Dfprotein['transcriptID']]['tag'] == 'basic' or 'level 1; protein_id' in DDtranscriptInfo[Dfprotein['transcriptID']]: #make equivalent of BASIC
						#len(DDtranscriptInfo[Dfprotein['transcriptID']]['level 1; protein_id']) > 0: #make equivalent of BASIC
							LDfproteinsClean.append(Dfprotein)
						else:
							iNumTranscriptStatusNotKnown += 1
							LDfproteinsDirty.append(Dfprotein)
							continue
					elif sSubsetTag != '':
						#tag "basic"; tag "Ensembl_canonical"; tag "GENCODE_Primary"; tag "MANE_Select"; tag "appris_principal_1";

						#GENCODE Comprehensive, which captures all annotated transcripts
						#GENCODE Basic, a smaller subset of transcripts containing only full-length transcripts at protein-coding genes and enriched for alternative splicing at lncRNA genes. The length and depth of the reads generated in ONT and PacBio experiments will eventually overwhelm the GENCODE Basic category, for example by adding many tens or hundreds of full-length alternatively spliced transcripts at protein-coding genes.
						#GENCODE Primary transcripts (v47 +) capture a minimal set of transcripts at protein coding genes and include protein-coding exons with evidence of evolutionary constraint and alternative splicing (including exon skips) with high expression, includes all MANE Select/Ensembl Canonical and MANE Plus Clinical transcripts by default.

						# Tier 1: Highest Priority (Highly Likely to Be Translated)These transcripts are functional alternatives that likely just missed the stringent Primary/MANE cutoffs (which favor a single representative transcript per gene).
							# 	APPRIS Principal Variants (appris_principal): Check the 9th column of your GTF for the appris_principal tag. The APPRIS database uses protein structure and conservation to flag principal functional isoforms. Even if a transcript isn't in GENCODE Primary, an APPRIS Principal tag means it preserves critical protein structures
							# 	Alternative APPRIS Tags (appris_candidate_long / appris_candidate): These represent alternative transcripts that still code for long, functional, and evolutionarily conserved protein domains
							# 	UniProt / Swiss-Prot Match: Cross-reference the Ensembl Protein IDs (ENSP) with Swiss-Prot. If the exact alternative sequence has a reviewed, manually curated entry in Swiss-Prot with "Evidence at protein level," it is actively translated.
						if "basic" in DDtranscriptInfo[Dfprotein['transcriptID']]['key_values']:
							icountBasic += 1
						if "GENCODE_Primary" in DDtranscriptInfo[Dfprotein['transcriptID']]['key_values']:
							iCountPrimary += 1
						if "basic" in DDtranscriptInfo[Dfprotein['transcriptID']]['key_values'] and "GENCODE_Primary" in DDtranscriptInfo[Dfprotein['transcriptID']]['key_values']:
							iCountBothBasicPrimary += 1
						if sSubsetTag in DDtranscriptInfo[Dfprotein['transcriptID']]['key_values']:
							LDfproteinsClean.append(Dfprotein)
					else:	#comprehensive or all in .gtf file
						LDfproteinsClean.append(Dfprotein)

				else:
					iNumNotProteinCoding += 1
					LDfproteinsDirty.append(Dfprotein)
		else:
			LDfproteinsUnannotated.append(Dfprotein)
			iNumProteinsUnannotated += 1

	iNumProteinsClean = len(LDfproteinsClean)
	print('<br># Proteins unannotated:\t', iNumProteinsUnannotated)
	print('<br># Proteins after Gencode cleaning:\t', iNumProteinsClean)
	print('<br>Number of entries with transcriptID mapped in the GTF file:\t',  iNumTranscriptIDsFound)
	print('<br>Number of entries with transcript_type mapped in the GTF file:\t', iNumTranscriptTypesFound)
	print('<br>Number of mapped entries removed that were not transcript_type:protein_coding:\t', iNumNotProteinCoding)
	print('<br>Number of mapped entries removed that were not transcript_status - KNOWN: ', iNumTranscriptStatusNotKnown)
	#print('<br>Number of mapped entries removed that were transcript_type:protein_coding, but not assigned to a chromosome:\t', iNumNotChromosomeAssigned)
	print('<br># Gencode PC proteins tagged GENCODE_Primary: ', iCountPrimary)
	print('<br># Gencode PC proteins tagged basic:', icountBasic)
	print('<br># Gencode PC proteins tagged both basic & primary:', iCountBothBasicPrimary)

	return(LDfproteinsClean, LDfproteinsDirty, LDfproteinsUnannotated)
##############################################################################
############################# parseEnsemblHeaders ############################
##############################################################################
#GRCH38
#>ENSP00000493203.1 pep chromosome:GRCh38:2:29192774:29223900:-1 gene:ENSG00000171094.18 transcript:ENST00000642122.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:ALK description:ALK receptor tyrosine kinase [Source:HGNC SymbolAcc:HGNC:427]

#GRCH37
#>ENSP00000328693 pep:novel chromosome:NCBI35:1:904515:910768:1 gene:ENSG00000158815:transcript:ENST00000328693 gene_biotype:protein_coding transcript_biotype:protein_coding
# ^               ^   ^     ^                                   ^                    ^                          ^                            ^ 
# ID              |   |  LOCATION                          GENE:stable gene ID       |                       GENE: gene biotype           TRANSCRIPT: transcript biotype
#                 | STATUS                                           TRANSCRIPT: stable transcript ID
#               SEQTYPE
def fParseEnsemblHeaders(LDfproteins):
	
	## parse original headers to extract info
	for Dfprotein in LDfproteins:
		LfieldsSpace = Dfprotein['header'].split()	#Python for split on white space, can also try  import re re.split('\s+', Dfprotein['header'])
		sAccession_number = LfieldsSpace[0]
		sAccession_number = sAccession_number.strip('>')	#strip leading >
		#my($status, $geneID, $transcriptID, $chromosome, $genomeCoordinates)
		#my($gene_symbol, $transcript_biotype)
		#for sFieldsSpace in LfieldsSpace:
		Dfprotein['proteinName'] = '' #initialize with empty value
		for i in range(1,len(LfieldsSpace)): 
		#for( my $i =1 $i <= $#AfieldsSpace $i++){
			LfieldsColon = LfieldsSpace[i].split(':')
			if(LfieldsColon[0] == 'pep'):
				if(len(LfieldsColon) > 1):		sStatus = LfieldsColon[1]	#typical of GChr37
				else:								sStatus = ''				#typical of GChr38
			elif(LfieldsColon[0] == 'gene'):		sGeneID = LfieldsColon[1]
			elif(LfieldsColon[0] == 'transcript'):	sTranscriptID = LfieldsColon[1]
			elif(LfieldsColon[0] == 'chromosome' or LfieldsColon[0] == 'scaffold'): #scaffold in GRCH38 but not GRCH37
				sChromosome = LfieldsColon[2]
				sGenomeCoordinates = 'chr' + sChromosome + ':' + LfieldsColon[3] + '-' + LfieldsColon[4] + ':'
				if(LfieldsColon[5]): #for strand convert +1 to +, and -1 to -
					sGenomeCoordinates += '+'
				else:
					sGenomeCoordinates += '-'
			elif(LfieldsColon[0] == 'transcript_biotype'):	Dfprotein['transcript_biotype'] = LfieldsColon[1]
			elif(LfieldsColon[0] == 'gene_symbol'):	#in GRCH38 but not GRCH37
				Dfprotein['gene_symbol'] = LfieldsColon[1]
			elif(LfieldsColon[0] == 'description'):	#in GRCH38 but not GRCH37
				#the last entry in the line is a key:value where value has spaces, unlike all other values, so re-join the space delimited description
				#$RHfprotein{'description'] = $AfieldsColon[1] + ' ' + join(' ', @AfieldsSpace[$i+1 .. $#AfieldsSpace])
				Dfprotein['description'] = LfieldsColon[1] + ' ' +  ' '.join(LfieldsSpace[i+1 : len(LfieldsSpace)])
				Lfields = Dfprotein['description'].split(' [Source') #trim off superfluous text
				Dfprotein['proteinName'] = Lfields [0]
				break
					#break prevents problems with encountering keywords in the description field,
					#especially  (chromosome, scaffold) causing the prior $chromosome and genomeCoordinates to be overwritten,
					#perhaps also (pep, gene, transcript, transcript_biotype, gene_symbol) 
		#end i loop
		Dfprotein['accession_number']	= Dfprotein['proteinID'] = sAccession_number
		Dfprotein['geneID']				= sGeneID
		Dfprotein['transcriptID']		= sTranscriptID
		Dfprotein['status']				= sStatus
		Dfprotein['chromosome']			= sChromosome
		Dfprotein['genomeCoordinates']	= sGenomeCoordinates
	#end Dfprotein loop

	return (LDfproteins)

##############################################################################
############### updateGRCh37EnsemblHeadersCreateCategoriesFile ###############
##############################################################################
#s 
def fUpdateGRCh37EnsemblHeadersCreateCategoriesFile(LDfproteins,  DDproteinID_info, sDatabasePath, iWriteToFile):

	##create new header that includes gene symbol and description
	for Dfprotein in LDfproteins:
		LfieldsDot = Dfprotein['accession_number'].split('.')
		Dfprotein['accNumNoVersion'] = LfieldsDot[0]
		Dfprotein['geneSymbol'] = DDproteinID_info[Dfprotein['accNumNoVersion']]['HGNC symbol']
		if(Dfprotein['geneSymbol'] == ''): #fill in with next best choice, when HGNC is empty
			Dfprotein['geneSymbol'] = DDproteinID_info[Dfprotein['accNumNoVersion']]['Associated Gene Name']
	
		#tidy up portein name, if necessary
		if not ('proteinName' in Dfprotein):
			if('description' in Dfprotein):	#KRC trimming now already done in fParseEnsemblBiomartFile
				Lfields = Dfprotein['description'].split('\[Source') #trim off superfluous text
				Dfprotein['proteinName'] = Lfields [0]

		Dfprotein['headerOriginal'] = Dfprotein['header']
		
		#>ENST00000488147.1_1_1:17367-18318:-|GN=WASH7P Pseudogene|nuORFdb_human
		#>ENSBTAP00000064819|ENSBTAG00000052228|GN=TRBV24-1 T cell receptor beta variable 24-1 [Source:HGNC SymbolAcc:HGNC:12203]
		
		Dfprotein['header'] = '>' + Dfprotein['accession_number'] + '|' + Dfprotein['transcriptID'] + '|' + Dfprotein['geneID']
		Dfprotein['header'] += '|' + Dfprotein['genomeCoordinates'] + '|' + 'GN=' + Dfprotein['geneSymbol'] + ' ' + Dfprotein['proteinName'] ##. '|' + 'Human'
		Dfprotein['header'] += '\n'
	
	##write categories file
	if(iWriteToFile):
		fWriteCategoriesFile_Ensembl(sDatabasePath, LDfproteins)
	
	return (LDfproteins)

##############################################################################
############### updateGRCh38EnsemblHeadersCreateCategoriesFile ###############
##############################################################################
#s 
def fUpdateGRCh38EnsemblHeadersCreateCategoriesFile(LDfproteins, sDatabasePath, iWriteToFile):
#e	
	##create new header that includes gene symbol and description
	for Dfprotein in LDfproteins:
#		LfieldsDot = Dfprotein['accession_number'].split('.') #need 
#		Dfprotein['accNumNoVersion'] = LfieldsDot[0]
		
#		if not ('proteinName' in Dfprotein):
#			Dfprotein['proteinName'] = ''
#			if('description' in Dfprotein):
#				Lfields = Dfprotein['description'].split(' [Source') #trim off superfluous text
#				Dfprotein['proteinName'] = Lfields [0]

		Dfprotein['headerOriginal'] = Dfprotein['header']
		#>ENSP00000493203.1 pep chromosome:GRCh38:2:29192774:29223900:-1 gene:ENSG00000171094.18 transcript:ENST00000642122.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:ALK description:ALK receptor tyrosine kinase [Source:HGNC SymbolAcc:HGNC:427]
		#>ENST00000488147.1_1_1:17367-18318:-|GN=WASH7P Pseudogene|nuORFdb_human
		#>ENSBTAP00000064819|ENSBTAG00000052228|GN=TRBV24-1 T cell receptor beta variable 24-1 [Source:HGNC SymbolAcc:HGNC:12203]
		
		Dfprotein['header'] = '>' + Dfprotein['accession_number'] + '|' + Dfprotein['transcriptID'] + '|' + Dfprotein['geneID']
		Dfprotein['header'] += '|' + Dfprotein['genomeCoordinates'] + '|' + 'GN=' + Dfprotein['gene_symbol'] + '|' + Dfprotein['proteinName'] ##. '|' + 'Human'
		Dfprotein['header'] += '\n'
	
	##write categories file
	if(iWriteToFile):
		fWriteCategoriesFile_Ensembl(sDatabasePath, LDfproteins)
	
	return (LDfproteins)

##############################################################################
######################### fWriteCategoriesFile_Ensembl #######################
##############################################################################
def fWriteCategoriesFile_Ensembl(sDatabasePath, LDfproteins, DDtranscriptIDinfo = {}, bWriteTranscriptDetailsColumns = False):
	
	sOutPath = sDatabasePath + '.categories.tsv'
	FileCategories = open(sOutPath , 'w')
	
	##write headers for output
	##no spaces allowed in header variables for subsequent use as SM category files
	##use speciesC so as not to overwrite variable name obtained in tagSummary file
		#KRC 8/30/3025
	FileCategories.write('accession_number\tgeneSymbol\tsequenceLength\tgeneID\ttranscriptID\ttranscript_biotype\tchromosome\tgenomeCoordinatesProtein\tgenomeCoordinatesTranscript\tredundant_accession_numbers\ttagBasic\ttagGENCODE_primary\tproteinName')
	if 'startNT4' in LDfproteins[0]:
		FileCategories.write('\tstartNT4')

	#KRC 5/20/2024 bring along the transcript info to see what's there
	if bWriteTranscriptDetailsColumns:
		Lkeys= []
		if bool(DDtranscriptIDinfo):
			Lkeys = DDtranscriptIDinfo[LDfproteins[0]['transcriptID']].keys()
			for sKey in Lkeys:
				if sKey == 'key_values': continue
				FileCategories.write('\t' + sKey)
	FileCategories.write('\theaderOriginal\n')
	
	#write lines
	for Dfprotein in LDfproteins:
		FileCategories.write(Dfprotein['accession_number'] + '\t')
		
		if('gene_symbol' in Dfprotein):	##in GRCH38 but not GRCH37
				FileCategories.write(Dfprotein['gene_symbol'])
		
		FileCategories.write('\t' + str(Dfprotein['sequenceLength']) )
		FileCategories.write( '\t' + Dfprotein['geneID']		+ '\t' + Dfprotein['transcriptID'])
		
		FileCategories.write( '\t')
		if('transcript_biotype' in Dfprotein):	#Ensembl
			FileCategories.write( Dfprotein['transcript_biotype'])
		elif('transcript_type' in Dfprotein):	#Gencode
			FileCategories.write( Dfprotein['transcript_type'])					

		FileCategories.write('\t')
		if 'chromosome' in Dfprotein:
			FileCategories.write(Dfprotein['chromosome'])

		#KRC 8/30/3025
		FileCategories.write('\t')
		if 'genomeCoordinatesProtein' in Dfprotein:
			 FileCategories.write(Dfprotein['genomeCoordinatesProtein'])

		FileCategories.write('\t')
		if 'genomeCoordinatesTranscript' in Dfprotein:
			 FileCategories.write(Dfprotein['genomeCoordinatesTranscript'])

		FileCategories.write( '\t')
		if('redundant_accession_numbers' in Dfprotein):	#can be populated only after making a database non-redundant
			FileCategories.write(Dfprotein['redundant_accession_numbers'])

		FileCategories.write( '\t')
		if('proteinName' in Dfprotein):
			FileCategories.write( Dfprotein['proteinName'])

		if 'startNT4' in Dfprotein:
			FileCategories.write( '\t')
			FileCategories.write( Dfprotein['startNT4'])

		if bWriteTranscriptDetailsColumns:
			if len(Lkeys) > 0:
				if Dfprotein['transcriptID'] in DDtranscriptIDinfo:
					for sKey in Lkeys:
						if sKey == 'key_values': continue
						FileCategories.write('\t')
						if sKey in DDtranscriptIDinfo[Dfprotein['transcriptID']]:
							FileCategories.write(DDtranscriptIDinfo[Dfprotein['transcriptID']][sKey])
				else:
					print ('Missing transcript info for: ', Dfprotein['transcriptID'])

		FileCategories.write( '\t')
		if('chromosome' in Dfprotein and Dfprotein['chromosome'] == ''): #sometimes it is empty, not sure why
			FileCategories.write(Dfprotein['headerOriginal'])
		else:
			FileCategories.write('\n')	# should add a '\n'

	FileCategories.close()
	print ('<br>Wrote file ' + sOutPath)
	
	return(0)
###############################################################################
############################## parseGencodeHeaders ############################
###############################################################################
##s 
def fParseGencodeHeaders(LDfproteins, iParseTranscripts = 0, DDtranscriptIDinfoBioMart = {}):
##e
	#Gencode.v34.pc_translations.fasta
	#>ENSP00000493203.1|ENST00000642122.1|ENSG00000171094.18|OTTHUMG00000152034.4|OTTHUMT00000493449.1|ALK-207|ALK|552
	#http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/_README.TXT
		#12. gencode.vX.pc_translations.fa.gz:
		#  Translations of protein-coding transcripts on reference chromosomes Fasta file.
		#  The sequence headers contain the following fields:
		#    transcript-id|
		#    gene-id|
		#    Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
		#    Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
		#    transcript-name|gene-name|sequence-length

	#gencode.v41.lncRNA_transcripts.fa
	#>ENST00000473358.1|ENSG00000243485.5|OTTHUMG00000000959.2|OTTHUMT00000002840.1|MIR1302-2HG-202|MIR1302-2HG|712|
	indexTranscript = 1
	if iParseTranscripts:
		indexTranscript = 0

	## parse original headers to extract info
	for Dfprotein in LDfproteins:
		Lfields = Dfprotein['header'].split('|')
		sAccession_number = Lfields[0]
		sAccession_number = sAccession_number.strip('>')	#strip leading >
		Dfprotein['accession_number']	= Dfprotein['proteinID'] = sAccession_number
		if iParseTranscripts:
			Dfprotein['transcriptID']	= sAccession_number
		else:
			Dfprotein['transcriptID']	= Lfields[indexTranscript]
		Dfprotein['geneID']				= Lfields[indexTranscript + 1]
		#Dfprotein['geneSymbol']			= Lfields[6] #used SM favored variable name:  geneSymbol
		Dfprotein['gene_symbol']			= Lfields[indexTranscript + 5] #used Ensembl favored variable name:  gene__symbol

	return (LDfproteins)
###############################################################################
############################ parseGencodeHeaders_hg19 #########################
###############################################################################
#KRC 5/20/2024 this function handles the absence of proteinIds in the GRCh37 (hg19) releases of Gencode (v <20)
def fParseGencodeHeaders_hg19(LDfproteins, DDcategories_transcriptID = {}):
##e
	#Gencode.v34.pc_translations.fasta
	#>ENSP00000493203.1|ENST00000642122.1|ENSG00000171094.18|OTTHUMG00000152034.4|OTTHUMT00000493449.1|ALK-207|ALK|552
	#Gencode.v19.pc_translations_basic.fasta
	#>ENST00000589042.1|ENSG00000155657.19|OTTHUMG00000154448.19|OTTHUMT00000450680.2|TTN-018|TTN|35991

	#http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/_README.TXT
		#12. gencode.vX.pc_translations.fa.gz:
		#  Translations of protein-coding transcripts on reference chromosomes Fasta file.
		#  The sequence headers contain the following fields:
		#    transcript-id|
		#    gene-id|
		#    Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
		#    Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
		#    transcript-name|gene-name|sequence-length

	## parse original headers to extract info
	iNumMappedToEnsembl = 0
	iNumNotMappedToEnsembl = 0
	LDfproteinsMapped = []
	for Dfprotein in LDfproteins:
		Lfields = Dfprotein['header'].split('|')
		indexTranscript = 0		
		sIdentifier = Lfields[0].strip('>')
		if sIdentifier.startswith('ENST'):
			if sIdentifier in DDcategories_transcriptID:
				Dfprotein['accession_number']	= Dfprotein['proteinID'] = DDcategories_transcriptID [sIdentifier]['accession_number']
				Dfprotein['transcriptID']	= sIdentifier
				Dfprotein['header'] = Dfprotein['header'].replace('>', '>' + Dfprotein['accession_number'] + '|')
				LDfproteinsMapped.append(Dfprotein)
				iNumMappedToEnsembl+=1
			else:
				#print('transcriptID not in categories file: ', sIdentifier)
				Dfprotein['accession_number'] = Dfprotein['transcriptID']	= sIdentifier
				iNumNotMappedToEnsembl+=1
		else:
			print('First identifier is not a transcriptID: ', sIdentifier)
		Dfprotein['geneID']				= Lfields[indexTranscript + 1]
		#Dfprotein['geneSymbol']			= Lfields[6] #used SM favored variable name:  geneSymbol
		Dfprotein['gene_symbol']			= Lfields[indexTranscript + 5] #used Ensembl favored variable name:  gene__symbol

	print('iNumMappedToEnsembl: ', iNumMappedToEnsembl)
	print('iNumNotMappedToEnsembl: ', iNumNotMappedToEnsembl)

	return (LDfproteinsMapped)
###############################################################################
################## fUpdateGencodeHeadersCreateCategoriesFile ##################
###############################################################################
 
def fUpdateGencodeHeadersCreateCategoriesFile(LDfproteins, DDtranscriptIDinfo, DDproteinIDinfo, sDatabasePath, iWriteToCategoriesFile = 0):

	##append to header chromosome, genome coordinates, transcript_type, gene symbol and description
	for Dfprotein in LDfproteins:
		Dfprotein['headerOriginal']		= Dfprotein['header']
		#		my($accNumNoVersion, $version) = split(/\./, Dfprotein['accession_number'])
#		Dfprotein['accNumNoVersion'] = $accNumNoVersion
		if not Dfprotein['transcriptID'] in DDtranscriptIDinfo:
			continue
		Dfprotein['chromosome']			= DDtranscriptIDinfo[Dfprotein['transcriptID']]['chromosome_name']
		Dfprotein['chromosome']			= Dfprotein['chromosome'].lstrip('chr')	#strip leading chr

		#KRC before 8/30/2025 genome coordinates of last DtranscriptInfo['feature_type'] in Gencode .GTF file
		#Dfprotein['genomeCoordinates']	= DDtranscriptIDinfo[Dfprotein['transcriptID']]['chromosome_name'] + ':' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_start_location'] + '-' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_end_location'] + ':' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_strand']

		#KRC 8/30/2025
		Dfprotein['genomeCoordinatesTranscript']	= DDtranscriptIDinfo[Dfprotein['transcriptID']]['chromosome_name'] + ':' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_start_location'] + '-' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_end_location'] + ':' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_strand']
		Dfprotein['genomeCoordinatesProtein']	= DDtranscriptIDinfo[Dfprotein['transcriptID']]['chromosome_name'] + ':'
		#put coordinates in coonventional low to high order
		#'transcript_type' TR_C_gene & IG_C_gene have a stop codon annotation but not a start codon
		if DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_strand']  == '+' and 'genomic_start_location_protein' in DDtranscriptIDinfo[Dfprotein['transcriptID']] and 'genomic_end_location_protein' in DDtranscriptIDinfo[Dfprotein['transcriptID']]:
			Dfprotein['genomeCoordinatesProtein']	+= DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_start_location_protein'] + '-' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_end_location_protein'] + ':' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_strand']
		elif DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_strand']  == '-' and 'genomic_start_location_protein' in DDtranscriptIDinfo[Dfprotein['transcriptID']] and 'genomic_end_location_protein' in DDtranscriptIDinfo[Dfprotein['transcriptID']]:
			Dfprotein['genomeCoordinatesProtein']	+= DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_end_location_protein'] + '-' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_start_location_protein'] + ':' + DDtranscriptIDinfo[Dfprotein['transcriptID']]['genomic_strand']

		Dfprotein['transcript_type']	= DDtranscriptIDinfo[Dfprotein['transcriptID']]['transcript_type'] 

		Dfprotein['header']				= Dfprotein['header'].rstrip('\n')

	#>ENSP00000493203.1|ENST00000642122.1|ENSG00000171094.18|OTTHUMG00000152034.4|OTTHUMT00000493449.1|ALK-207|ALK|552
		
		Dfprotein['header'] += '|' + Dfprotein['genomeCoordinatesProtein']
		Dfprotein['header'] += '|' + 'PE=1'								# PE=1 to aid downstream Fragpipe, SN subgroup specific FDR control
		Dfprotein['header'] += '|' + 'GN=' + Dfprotein['gene_symbol']	# GN= format to aid downstream SM parsing
		if(Dfprotein['accession_number'] in DDproteinIDinfo and 'proteinName' in DDproteinIDinfo[Dfprotein['accession_number']]):
			Dfprotein['proteinName'] = DDproteinIDinfo[Dfprotein['accession_number']]['proteinName']
			Dfprotein['header'] += '|' + Dfprotein['proteinName']
		Dfprotein['header'] += '\n'

	##write categories file
	if(iWriteToCategoriesFile):
		fWriteCategoriesFile_Ensembl(sDatabasePath, LDfproteins)

	return (LDfproteins)
