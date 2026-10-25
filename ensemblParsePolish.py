#####################################################
#	Copyright (c) 2020-2026 Broad Institute of MIT and Harvard
#		All rights reserved
#####################################################

#Filename:	ensemblParsePolish.py
#Purpose:	Parse & Polish headers of FASTA formatted Ensembl/Gencode-derived sequence databases.
#Authors:	Karl Clauser 
#Created:	June 2020

#Python libs
import sys
import re	#for case-insensitive string search
import os
import gzip
from libSpecMill import lsmFastaUtilities	as lFU
from libSpecMill import lsmEnsemblUtilities as lEU

#revise the static variables when updating to new version. ToDo make these values CLI arguments
##############################################################################
########################## Static Global Variables ###########################
##############################################################################

#Update these 4 for new version
#KRC 9/7/2026 now expected to be initialized when running fMain() via call to fSetGlobalVariables()
sREF_PROTEOME_DIR		= 'Gencode47_Ensembl113\\' #from $seqdb_dir
sGENCODE_VERSION		= 'v47'
sENSEMBLVERSION			= 'Ensembl-113'
sSEQDB_DIR				= '../seqdb/'

#sREF_PROTEOME_DIR		= '../../' + sREF_PROTEOME_DIR	#CLI from seqdb/HLA/B721221_germline
#sREF_PROTEOME_DIR		= sSEQDB_DIR + sREF_PROTEOME_DIR	 #CLI SpectrumMill>
#sREF_PROTEOME_DIR		= '../' + sREF_PROTEOME_DIR	 #CLI processing-scripts from \\helium\proteomic_storage_vast\CPTAC4\PCC\proteinSequenceDatabases\processing-scripts
#sREF_PROTEOME_DIR		= os.path.join(r'P:\CPTAC4\PCC\proteinSequenceDatabases\\', sREF_PROTEOME_DIR)	 #CLI processing-scripts from Karl's machine
#sREF_PROTEOME_DIR		= os.path.join( r'D:\helium\vast_CPTAC4_PCC\proteinSequenceDatabases\\', sREF_PROTEOME_DIR)	 #CLI processing-scripts db's local

#Keep Fixed unless processing something other than Gencode
#sGRCVERSION			= '_GRCm39_'		#mouse GRCm39 #human GRCh38
#sGRCVERSION			= '_GRCh37_ENST_' 	#human GRCh37
#sREF_PROTEOME_BASE		= 'Ensembl_Homo_sapiens.GRCh37.pep.all.fasta'
#sREF_PROTEOME_BASE		= 'Ensembl.Homo_sapiens.GRCh38.pep.all.fasta'
sGRCVERSION				= '_GRCh38_'		##human GRCh38
sREF_PROTEOME_BASE			= 'Gencode.' + sGENCODE_VERSION + '.pc_translations.fasta'	#for output change .fa to .fasta
sREF_PROTEOME_BASE_GZ		= 'gencode.' + sGENCODE_VERSION + '.pc_translations.fa.gz'	#for input
sGENCODE_ANNOTATION_BASIC	= 'gencode.' + sGENCODE_VERSION + '.basic.annotation.gtf.gz'
sGENCODE_ANNOTATION_PRIMARY	= 'gencode.' + sGENCODE_VERSION + '.primary_assembly.annotation.gtf.gz'
sGENCODE_ANNOTATION_ALL		= 'gencode.' + sGENCODE_VERSION + '.annotation.gtf.gz'
sGENCODE_TRANSCRIPTS		= 'gencode.' + sGENCODE_VERSION + '.transcripts.fa.gz'

sGENCODE_TRANSCRIPTS_LNCRNA	= 'gencode.' + sGENCODE_VERSION + '.lncRNA_transcripts.fa.gz'
sGENCODE_ANNOTATION_LNCRNA	= 'gencode.' + sGENCODE_VERSION + '.long_noncoding_RNAs.gtf.gz'

##############################################################################
########################### fSetGlobalVariables ##############################
##############################################################################
# KRC 8/18/2026 when changing a global variable inside a function, must declare its global rather than local 
def	fSetGlobalVariables(sRefProteomeDir, sGencodeVersion, sEnsemblVersion, sSeqDBdir = r'../seqdb/'):

	#user supplied
	global sREF_PROTEOME_DIR, sGENCODE_VERSION, sENSEMBLVERSION, sSEQDB_DIR
	sREF_PROTEOME_DIR		= sRefProteomeDir
	sGENCODE_VERSION		= sGencodeVersion
	sENSEMBLVERSION			= sEnsemblVersion
	sSEQDB_DIR				= sSeqDBdir

	# dependent on sGENCODE_VERSION
	global	sREF_PROTEOME_BASE, sREF_PROTEOME_BASE_GZ
	global sGENCODE_ANNOTATION_BASIC, sGENCODE_ANNOTATION_PRIMARY, sGENCODE_ANNOTATION_ALL, sGENCODE_ANNOTATION_LNCRNA
	global sGENCODE_TRANSCRIPTS, sGENCODE_TRANSCRIPTS_LNCRNA
	sREF_PROTEOME_BASE			= 'Gencode.' + sGENCODE_VERSION + '.pc_translations.fasta'	#for output change .fa to .fasta
	sREF_PROTEOME_BASE_GZ		= 'gencode.' + sGENCODE_VERSION + '.pc_translations.fa.gz'	#for input
	#sGENCODE_ANNOTATION			= 'gencode.' + sGENCODE_VERSION + '.basic.annotation.gtf.gz'
	#sGENCODE_ANNOTATION		= 'gencode.' + sGENCODE_VERSION + '.annotation.gtf.gz'	#basic does not exist for v19, introduced with v23
	sGENCODE_ANNOTATION_BASIC	= 'gencode.' + sGENCODE_VERSION + '.basic.annotation.gtf.gz'
	sGENCODE_ANNOTATION_PRIMARY	= 'gencode.' + sGENCODE_VERSION + '.primary_assembly.annotation.gtf.gz' # introduced with v47?
	sGENCODE_ANNOTATION_ALL		= 'gencode.' + sGENCODE_VERSION + '.annotation.gtf.gz'
	sGENCODE_TRANSCRIPTS		= 'gencode.' + sGENCODE_VERSION + '.transcripts.fa.gz'

	sGENCODE_TRANSCRIPTS_LNCRNA	= 'gencode.' + sGENCODE_VERSION + '.lncRNA_transcripts.fa.gz'
	sGENCODE_ANNOTATION_LNCRNA	= 'gencode.' + sGENCODE_VERSION + '.long_noncoding_RNAs.gtf.gz'

	return(0)

##############################################################################
################################### main #####################################
##############################################################################
#no script flow outside of functions other than __main__ at end of file
def fMain(sSubsetTagGencode = 'basic', sRefProteomeDir = '', sGencodeVersion = '', sEnsemblVersion = '', sSeqDBdir = ''): #sSubsetTag: 'basic', 'GENCODE_Primary'

	fSetGlobalVariables(sRefProteomeDir, sGencodeVersion, sEnsemblVersion, sSeqDBdir)
	print( sGENCODE_VERSION, sENSEMBLVERSION, sSEQDB_DIR, '\n', sREF_PROTEOME_DIR)

	iCGIparams = 0 #set to 0 to run w/o CGI parameters
	if(iCGIparams):
		#&initializeHTML("Parse & Polish Ensembl", 3, "outputPane");
		#sREF_PROTEOME_DIR = sSeqdb_dir = fReadCGIparameters()
		os.chdir(sREF_PROTEOME_DIR)
	
	#read in the .fasta file to parse & polish
	sFastaPathIn	= os.path.join(sREF_PROTEOME_DIR, sREF_PROTEOME_BASE_GZ)
	#sFastaPathIn	= os.path.join(sREF_PROTEOME_DIR, sREF_PROTEOME_BASE)
	sFastaPathOut	= os.path.join(sREF_PROTEOME_DIR, sREF_PROTEOME_BASE)
	LDfproteins = lFU.fReadFastaFile(sFastaPathIn)
	
	LDfproteins = lFU.filterOutProteinsByMinSequenceLength(LDfproteins, 7);

	#sFastaFilenameClean = lFile.fRemoveFilenameSuffix(sFastaFilename)

	#for GRCh37, base headers lack info retreived from Biomart
	if(sREF_PROTEOME_BASE.find('GRCh37') > 0):
		#KRC 5/20/204 Make a acategories file of all entries
		LDfproteins = lEU.fParseEnsemblHeaders(LDfproteins)
		lEU.fWriteCategoriesFile_Ensembl(sFastaPathOut, LDfproteins)
		return(0)
				
		print('Cleaning up GRCh37')
		#ToDo KRC port from Perl
		LDfproteinsClean, LDfproteinsDirty = lEU.fCleanEnsemblFasta(LDfproteins) #keep only transcript_biotype:protein_coding
		sFastaFilenameClean = sFastaPathOut.replace('.fasta', '_clean2PY.fasta')
		lFU.fWriteFastaFile2(LDfproteinsClean, sFastaFilenameClean);

		#read in file containing geneSymbols and names
		#http://feb2014.archive.ensembl.org/biomart/martview/615bfdb5e5b8cc68780ffb3a2f247c89
		sBioMartFile = sREF_PROTEOME_DIR + 'Ensembl-75_GRCh37_ENST_mart_export.txt'
		DDproteinIDinfo = lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Ensembl Protein ID', 'Description')

		##parse header to extract info, simplify headers, and update with gene Symbols and names
		##write category file
		sFastaFilenameClean = sFastaPathOut.replace('.fasta', '_clean3PY.fasta')
		LDfproteinsClean = lEU.fParseEnsemblHeaders(LDfproteinsClean)
		LDfproteinsClean = lEU.fUpdateGRCh37EnsemblHeadersCreateCategoriesFile(LDfproteinsClean, DDproteinIDinfo, sFastaFilenameClean, 1);
		lFU.fWriteFastaFile2(LDfproteinsClean, sFastaFilenameClean)
	elif(sREF_PROTEOME_BASE.find('GRCh38') > 0 ):
		print('Cleaning up GRCh38')
		LDfproteinsClean, LDfproteinsDirty = lEU.fCleanEnsemblFasta(LDfproteins) #keep only transcript_biotype:protein_coding
		sFastaFilenameClean = sFastaPathOut.replace('.fasta', '_clean2PY.fasta')
		lFU.fWriteFastaFile2(LDfproteinsClean, sFastaFilenameClean);

		##parse header to extract info, simplify headers, and update with gene symbols and names
		##write category file
		sFastaFilenameClean = sFastaPathOut.replace('.fasta', '_clean3PY.fasta')
		LDfproteinsClean = lEU.fParseEnsemblHeaders(LDfproteinsClean)
		LDfproteinsClean = lEU.fUpdateGRCh38EnsemblHeadersCreateCategoriesFile(LDfproteinsClean, sFastaFilenameClean, 1)
		lFU.fWriteFastaFile2(LDfproteinsClean, sFastaFilenameClean);
		
		#write category file of removed proteins to help understand what they were
		LDfproteinsDirty = lEU.fParseEnsemblHeaders(LDfproteinsDirty)
		sFastaFilenameDirty= sFastaPathOut.replace('.fasta', '_dirty2.fasta')
		lEU.fWriteCategoriesFile_Ensembl(sFastaFilenameDirty, LDfproteinsDirty)
	elif(re.search('Gencode', sREF_PROTEOME_BASE, re.IGNORECASE)): #case-insensitive
		#read in file containing protein names, genome coordinates, transcript_type
		#https://www.gencodegenes.org/human/

		#select specific file to use
		#GENCODE Comprehensive, which captures all annotated transcripts
		#GENCODE Basic, a smaller subset of transcripts containing only full-length transcripts at protein-coding genes and enriched for alternative splicing at lncRNA genes. The length and depth of the reads generated in ONT and PacBio experiments will eventually overwhelm the GENCODE Basic category, for example by adding many tens or hundreds of full-length alternatively spliced transcripts at protein-coding genes.
		#GENCODE Primary transcripts (v47 +) capture a minimal set of transcripts at protein coding genes and include protein-coding exons with evidence of evolutionary constraint and alternative splicing (including exon skips) with high expression, includes all MANE Select/Ensembl Canonical and MANE Plus Clinical transcripts by default.
		sGenocodeAnnotationFile = sGENCODE_ANNOTATION_ALL				#Comprehensive
		if   sSubsetTagGencode == 'basic':
			sGenocodeAnnotationFile = sGENCODE_ANNOTATION_BASIC
		elif sSubsetTagGencode == 'GENCODE_Primary':
			sGenocodeAnnotationFile = sGENCODE_ANNOTATION_PRIMARY
		elif sSubsetTagGencode == 'basicAndPrimary':
			sSubsetTagGencode = 'GENCODE_Primary'
			sFastaPathOut = sFastaPathOut.replace('.fasta', '_basic.fasta')
			sGenocodeAnnotationFile = sGENCODE_ANNOTATION_BASIC

		sGencodeAnnotationPath = os.path.join(sREF_PROTEOME_DIR, sGenocodeAnnotationFile)
		#DDtranscriptIDinfo = lEU.fParseGencodeGTFfile(sGencodeAnnotationPath)
		DDtranscriptIDinfo = lEU.fParseGencodeGTFfile_CDS(sGencodeAnnotationPath)

		#read in file containing protein names
		#https://www.ensembl.org/biomart/
		#archive.ensembl.org to find whatever version you need
		#http://apr2020.archive.ensembl.org/biomart/martview/a166902a60e7b484b5b8e01e8fc4ab44 #Ensembl 100/Gencode 34
		#sBioMartFile = sREF_PROTEOME_DIR + 'Ensembl-104_GRCh38_2mart_export.txt'	
		sBioMartFile = os.path.join(sREF_PROTEOME_DIR, sENSEMBLVERSION + sGRCVERSION + 'mart_export.txt')	#KRC 6/19/2021 for Gencode polishing
		DDtranscriptIDinfoBioMart = {}	#emptydictionary
		bv19Filter = False
		if 'v19' in sREF_PROTEOME_BASE: #hg19 versions of Gencode lacked protein IDs, so fill that in from the Biomart file
			DDproteinIDinfo				= lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Ensembl Protein ID', 'Description')
			DDtranscriptIDinfoBioMart	= lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Ensembl Transcript ID', 'Description')
			bv19Filter = True
		else:
			DDproteinIDinfo = lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Protein stable ID version', 'Gene description')

		#parse existing headers
		if 'v19' in sREF_PROTEOME_BASE: #hg19 versions of Gencode lacked protein IDs, so fill that in from the Ensembl file
			LDcategories = lEU.fParseTSVfileHeader(sREF_PROTEOME_DIR + 'Ensembl_Homo_sapiens.GRCh37.pep.all.fasta.categories.tsv')	#
			DDcategories_transcriptID = lFU.fListOfDictionariesToDD(LDcategories, 'transcriptID')
			LDfproteins = lEU.fParseGencodeHeaders_hg19(LDfproteins, DDcategories_transcriptID)
			
		else:
			LDfproteins = lEU.fParseGencodeHeaders(LDfproteins, 0);
			
		lEU.fWriteCategoriesFile_Ensembl(sFastaPathOut, LDfproteins, DDtranscriptIDinfo)

		#clean to keep only transcript_biotype:protein_coding, #sSubsetTag: 'basic', 'GENCODE_Primary'
		LDfproteinsClean, LDfproteinsDirty, LDfproteinsUnannotated = lEU.fCleanGencodeFasta(LDfproteins, DDtranscriptIDinfo, bv19Filter, sSubsetTag = sSubsetTagGencode) 
		#for output filenames strip the Gencode part, if present, and capitalize 1st letter
		sSubsetTagGencode = sSubsetTagGencode.replace('GENCODE_', '').capitalize()
		sFastaFilenameClean = sFastaPathOut.replace('.fasta', '_' + sSubsetTagGencode + '.fasta')	 #clean2
		lFU.fWriteFastaFile2(LDfproteins, sFastaFilenameClean);

		##parse header to extract info, simplify headers, and update with genome coordinates and names
		##write category file
#		sFastaFilenameClean = lFile.fRemoveFilenameSuffix(sFastaFilename)
#		sFastaFilenameClean += "_clean3PY.fasta";
		sFastaFilenameClean = sFastaPathOut.replace('.fasta', '_' + sSubsetTagGencode + 'PC.fasta')	 #clean3
		LDfproteinsClean = lEU.fUpdateGencodeHeadersCreateCategoriesFile(LDfproteinsClean, DDtranscriptIDinfo, DDproteinIDinfo, sFastaFilenameClean, 1)
		lEU.fWriteCategoriesFile_Ensembl(sFastaFilenameClean, LDfproteinsClean, DDtranscriptIDinfo)
		lFU.fWriteFastaFile2(LDfproteinsClean, sFastaFilenameClean)
		
		#write category file of removed proteins to help understand what they were
		#lEU.fWriteCategoriesFile_Ensembl(sFastaFilenameDirty, LDfproteinsDirty)
		sFastaFilenameDirty= sFastaPathOut.replace('.fasta', '_' + sSubsetTagGencode + 'NotPC.fasta')	 #dirty2
		LDfproteinsDirty = lEU.fUpdateGencodeHeadersCreateCategoriesFile(LDfproteinsDirty, DDtranscriptIDinfo, DDproteinIDinfo, sFastaFilenameDirty, 1)

		#unannotated .fasta and .category file
		#DDtranscriptIDinfoAll = lEU.fParseGencodeGTFfile(sREF_PROTEOME_DIR + sGENCODE_ANNOTATION_ALL)
		DDtranscriptIDinfoAll = lEU.fParseGencodeGTFfile_CDS(os.path.join(sREF_PROTEOME_DIR, sGENCODE_ANNOTATION_ALL))
		sFastaFilenameUnannotated= sFastaPathOut.replace('.fasta', '_not' + sSubsetTagGencode + '.fasta')
		LDfproteinsUnannotated = lEU.fUpdateGencodeHeadersCreateCategoriesFile(LDfproteinsUnannotated, DDtranscriptIDinfoAll, DDproteinIDinfo, sFastaFilenameUnannotated, 1)
		lFU.fWriteFastaFile2(LDfproteinsUnannotated, sFastaFilenameUnannotated)

	## Make a non-redundant version	of clean, along with ta list of the removed accession number matched to their retained counterparts
	sFastaFilenameCleanNR = sFastaFilenameClean.replace('.fasta', 'nr.fasta')		#clean3nr
	sMatchesFile = sFastaFilenameCleanNR + '.matches.txt'	#.txt instead of .tsv so less clutter for categories files
	LDfproteinsCleanNR, LDfproteinsCleanRonly = lFU.fMakeFastaNRdatabase(LDfproteinsClean, sFastaFilenameCleanNR, sMatchesFile)
	lEU.fWriteCategoriesFile_Ensembl(sFastaFilenameCleanNR, LDfproteinsCleanNR, DDtranscriptIDinfo, bWriteTranscriptDetailsColumns = False)
	
	return(0)
##############################################################################
############################## fProcess_lncRNA ###############################
##############################################################################
def fProcess_lncRNA():
	#read in the .fasta file to parse & polish
	sFastaFilename = sREF_PROTEOME_DIR + sGENCODE_TRANSCRIPTS_LNCRNA
	LDfTranscripts = lFU.fReadFastaFile(sFastaFilename)

	if(re.search('Gencode', sREF_PROTEOME_BASE, re.IGNORECASE)): #case-insensitive

		sBioMartFile = sREF_PROTEOME_DIR + sENSEMBLVERSION + '_GRCh38_mart_export.txt'	#KRC 6/19/2021 for Gencode polishing
		#DDproteinIDinfo = lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Protein stable ID version', 'Gene description')
		DDproteinIDinfo = lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Transcript stable ID version', 'Gene description')

		DDtranscriptIDinfo = lEU.fParseGencodeGTFfile(sREF_PROTEOME_DIR +  sGENCODE_ANNOTATION_LNCRNA)
		LDfTranscripts = lEU.fParseGencodeHeaders(LDfTranscripts, 1)
		sFastaFilenameLncRNA= sFastaFilename.replace('.fa.gz', '_polished.fasta')
		LDfTranscripts = lEU.fUpdateGencodeHeadersCreateCategoriesFile(LDfTranscripts, DDtranscriptIDinfo, DDproteinIDinfo, sFastaFilenameLncRNA, 1)
		lFU.fWriteFastaFile2(LDfTranscripts, sFastaFilenameLncRNA);

	return(0)
##############################################################################
########################### fProcess_Transcripts #############################
##############################################################################
def fProcess_Transcripts(bFilterToProteinCoding = False):
	#read in the .fasta file to parse & polish
	sFastaFilename = sREF_PROTEOME_DIR + sGENCODE_TRANSCRIPTS
	LDfTranscripts = lFU.fReadFastaFile(sFastaFilename)

	sFastaSuffix = ''
	if(re.search('Gencode', sREF_PROTEOME_BASE, re.IGNORECASE)): #case-insensitive

		sBioMartFile = sREF_PROTEOME_DIR + sENSEMBLVERSION + '_GRCh38_mart_export.txt'	#KRC 6/19/2021 for Gencode polishing
		#DDproteinIDinfo = lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Protein stable ID version', 'Gene description')
		DDproteinIDinfo = lEU.fParseEnsemblBiomartFile(sBioMartFile, 'Transcript stable ID version', 'Gene description')

		DDtranscriptIDinfo = lEU.fParseGencodeGTFfile(sREF_PROTEOME_DIR +  sGENCODE_ANNOTATION_ALL)
		LDfTranscripts = lEU.fParseGencodeHeaders(LDfTranscripts, 1)

		if bFilterToProteinCoding:
			for DfTranscript in LDfTranscripts:
				DfTranscript['startNT4'] = DfTranscript['sequence'][: 4]
				DfTranscript['transcriptID'] =  DfTranscript['accession_number'] #needed for clean
			LDfTranscriptsClean, LDfTranscriptsDirty, LDfTranscriptsUnannotated = lEU.fCleanGencodeFasta(LDfTranscripts, DDtranscriptIDinfo, bv19Filter = False) #keep only transcript_biotype:protein_coding
			LDfTranscripts = LDfTranscriptsClean
			sFastaSuffix += '_basic'

		sFastaFilenameTranscripts= sFastaFilename.replace('.fa.gz', sFastaSuffix + '_polished.fasta')
		LDfTranscripts = lEU.fUpdateGencodeHeadersCreateCategoriesFile(LDfTranscripts, DDtranscriptIDinfo, DDproteinIDinfo, sFastaFilenameTranscripts, 1)
		lFU.fWriteFastaFile2(LDfTranscripts, sFastaFilenameTranscripts);

	return(0)

##############################################################################
############################## fAnnotationSummary ############################
##############################################################################
def fAnnotationSummary(sGencodeAnnotationPath):

	#read in Gencode .gtf file
	#sGencodeAnnotationPath = sREF_PROTEOME_DIR + sGENCODE_ANNOTATION
	#dfAnnotation = pd.read_table(sGencodeAnnotationPath, delimiter= '\t', comment='#', header=0 )
	N = 500
	sTopNpath = sGencodeAnnotationPath.replace('gtf.gz', str(N) + '.txt')
	FilepTopN = open(sTopNpath, 'w')

	#reading
	with gzip.open(sGencodeAnnotationPath, 'rt') as f:
		#head = [next(f) for x in range(N)]
		#head=f.readlines()[0:N] #put here the interval you want
		n = 0
		for sLine in f:
			if 'protein_coding' not in sLine:
				continue
			if N == n:
				break
			FilepTopN.write(sLine)
			#print (line)
			n = n+1
	#print (head)
	FilepTopN.close()
	return(0)

##############################################################################
#no script flow outside of functions other than main at end of file
##############################################################################
################################## __main__ ##################################
##############################################################################
if __name__ == "__main__":	#allow this script to be loaded as a library w/o executing this code block

	#fMain()
	#fProcess_lncRNA()
	#fProcess_Transcripts(bFilterToProteinCoding = True)
	#fAnnotationSummary(sREF_PROTEOME_DIR + sGENCODE_ANNOTATION_LNCRNA)
	#fAnnotationSummary(os.path.join(sREF_PROTEOME_DIR,sGENCODE_ANNOTATION))
	#fAnnotationSummary(r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116\gencode.v50.primary_assembly.basic.annotation.gtf.gz')
	#fAnnotationSummary(r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116\gencode.v50.primary_assembly.annotation.gtf.gz')
	#fAnnotationSummary(r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116\gencode.v50.basic.annotation.gtf.gz')

	#KRC 8/19/2026 1st try at primary set
	#fMain(sSubsetTagGencode = 'GENCODE_Primary', sRefProteomeDir = r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116', sGencodeVersion = 'v50', sEnsemblVersion = 'Ensembl-116', sSeqDBdir = '')	  #sSubsetTag: 'basic', 'GENCODE_Primary', 'basicAndPrimary'

	#KRC 8/24/2026 primary subset of basic
	#fMain(sSubsetTagGencode = 'basicAndPrimary', sRefProteomeDir = r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116', sGencodeVersion = 'v50', sEnsemblVersion = 'Ensembl-116', sSeqDBdir = '')	  #sSubsetTag: 'basic', 'GENCODE_Primary', 'basicAndPrimary'

	#KRC 9/7/2026
	#fMain(sSubsetTagGencode = 'basicAndPrimary', sRefProteomeDir = r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116', sGencodeVersion = 'v50', sEnsemblVersion = 'Ensembl-116', sSeqDBdir = '')	  #sSubsetTag: 'basic', 'GENCODE_Primary', 'basicAndPrimary'
	fMain(sSubsetTagGencode = 'GENCODE_Primary', sRefProteomeDir = r'P:\CPTAC4\PCC\proteinSequenceDatabases\Gencode50_Ensembl116', sGencodeVersion = 'v50', sEnsemblVersion = 'Ensembl-116', sSeqDBdir = '')	  #sSubsetTag: 'basic', 'GENCODE_Primary', 'basicAndPrimary'

	sys.exit(0)
