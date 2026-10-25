##############################################################################
#	Copyright (c) 2020-2024 Broad Institute of MIT and Harvard  
#		All rights reserved
##############################################################################

#Filename:	lsmFastaUtilities.py
#Purpose:	Utilities for reading, writing, filtering content, and retreiving entries from FASTA formatted sequence databases.
#Authors:	Karl Clauser 
#Created:	June 2020
import gzip
import os
import sys
from libSpecMill import lsmFastaHeaders as lFH

##############################################################################
############################## Static Variables ##############################
##############################################################################
sInitializedDataBase = ''
DaccToIndex = {}
LDindexOffsets = []
DDallProteinsByAccessionKey= {}	#initialize to empty dictionary

##############################################################################
############################### fReadFastaFile ###############################
##############################################################################
def fReadFastaFile(sFile_in, bVerbose = True):

	if not sFile_in.endswith('.gz'):
		FileIn = open(sFile_in, 'r')
		Llines = FileIn.read().splitlines()	#list of lines with \n removed
	else:
		with gzip.open(sFile_in, 'rt') as FileIn:
			Llines = FileIn.read().splitlines()	#list of lines with \n removed

	iEntryCounter = 0
	Dfprotein = fMakeNewFastaProteinDictionary() #dummy 1st protein
	LDfproteins = [] #initialize empty list
	for sLine in Llines:
		if (sLine == ''):	 #blank line in file
			continue

		if (sLine.startswith('>') ):	#starts with >
			#store previous protein
			if(Dfprotein['sequenceLength'] > 0):
				LDfproteins.append(Dfprotein)
			Dfprotein = fMakeNewFastaProteinDictionary()
			Dfprotein['header'] = sLine
			iEntryCounter += 1	#no ++ in Python
			if ( iEntryCounter % 50000 == 0 ) and bVerbose:		print ('<br>', iEntryCounter, ' entries read')
		else:
			#append lines to the sequence so that sequence is a single continuous string
			#sLine = sLine.rstrip()	#probably not necessary as readlines already removes \n
			Dfprotein['sequence'] += sLine
			Dfprotein['sequenceLength'] = len(Dfprotein['sequence'])
	
	#include last one
	if(Dfprotein['header'] != "" and Dfprotein['sequence'] != ""): #skip blank last line
		LDfproteins.append(Dfprotein)

	print ('<br>', len(LDfproteins), ' entries read ', sFile_in)
	FileIn.close()

	return(LDfproteins)
##############################################################################
############################# fReadFastaFileAcc ##############################
##############################################################################
#non-verbose version that also parses header for accession number
def fReadFastaFileAcc(sFastaPath):

	if not sFastaPath.endswith('.gz'):
		FileIn = open(sFastaPath, 'r')
		Llines = FileIn.read().splitlines()	#list of lines with \n removed
	else:
		with gzip.open(sFastaPath, 'rt') as FileIn:
			Llines = FileIn.read().splitlines()	#list of lines with \n removed

	iEntryCounter = 0
	Dfprotein = fMakeNewFastaProteinDictionary() #dummy 1st protein
	LDfproteins = [] #initialize empty list
	for sLine in Llines:
		if (sLine == ''):	 #blank line in file
			continue

		if (sLine.startswith('>') ):	#starts with >
			#store previous protein
			if(Dfprotein['sequenceLength'] > 0):
				LDfproteins.append(Dfprotein)
			Dfprotein = fMakeNewFastaProteinDictionary()
			Dfprotein['header'] = sLine
			iEntryCounter += 1	#no ++ in Python
			#if ( iEntryCounter % 50000 == 0 ) and bVerbose:		print ('<br>', iEntryCounter, ' entries read')
		else:
			#append lines to the sequence so that sequence is a single continuous string
			#sLine = sLine.rstrip()	#probably not necessary as readlines already removes \n
			Dfprotein['sequence'] += sLine
			Dfprotein['sequenceLength'] = len(Dfprotein['sequence'])
	
	#include last one
	if(Dfprotein['header'] != "" and Dfprotein['sequence'] != ""): #skip blank last line
		LDfproteins.append(Dfprotein)

	sFastaFile = os.path.basename(sFastaPath)
	print ('<br>', len(LDfproteins), ' entries read ', sFastaFile)
	FileIn.close()

	#parse header for accession number
	print('Parsing headers from: ', sFastaFile)
	for Dfprotein in LDfproteins:
		Dfprotein['accession_number'] = lFH.fAccessionNumberFromHeader(Dfprotein['header'], sFastaFile)

	return(LDfproteins)

##############################################################################
######################## fMakeNewFastaProteinDictionary ######################
##############################################################################
def fMakeNewFastaProteinDictionary():

#	Dfprotein = {}
#	Dfprotein['header'] = ""
#	Dfprotein['sequence'] = ""
#	Dfprotein['sequenceLength'] = 0
#	Dfprotein['accession_number'] = 0
	Dfprotein = {
		'header'			: '',
		'sequence'			: '',
		'sequenceLength'	: 0,
		'accession_number'	: 0
	}

	return(Dfprotein)

##############################################################################
############################### fWriteFastaFile2 #############################
##############################################################################
def fWriteFastaFile2(LDfproteins, sFilePath, iLineLength = 80):

	#open ($outputFileHandle, ">$databasePath") or &fatalErrorHTML("Can't open $databasePath\n", "lsmFastaParser.pl");
	FileFasta = open(sFilePath, 'w')

	iNumProteinsFiltered = 0
	for Dfprotein in LDfproteins:
		fAppendEntryToOpenFastaFile(FileFasta, Dfprotein['header'], Dfprotein['sequence'], iLineLength);
		iNumProteinsFiltered += 1	#no ++ in Python

	FileFasta.close()
	print ("\nFinished writing file containing ", iNumProteinsFiltered, " proteins: ", sFilePath);

	return (0)

##############################################################################
######################## fAppendEntryToOpenFastaFile #########################
##############################################################################
def fAppendEntryToOpenFastaFile(FileFasta, sHeader, sSequence, iLineLength = 80):
  
	if not(sHeader.startswith('>') ): #add > if missing
		sHeader = '>' + sHeader

	FileFasta.write(sHeader)	#sHeader should already have a \n
	if not sHeader.endswith('\n'):
		FileFasta.write('\n')
	#split sequence into lines of length 80. The last one should be less	
	#new_list = [expression for member in iterable] 	list comprehension syntax

	if (len(sSequence) < 1):
		print ('\n<br><b>Error: </b>Sequence is empty. \n<br>', sHeader, '\n')
		sys.exit(-1)

	elif (iLineLength > 0 ):
		#LsequenceChunks = [sSequence[j:j+iLineLength] + '\n' for j in range(0, len(sSequence), iLineLength)]
		LsequenceChunks = [sSequence[j:j+iLineLength] for j in range(0, len(sSequence), iLineLength)]
		#FileFasta.writelines(LsequenceChunks)		#writelines does not add \n
		FileFasta.write('\n'.join(LsequenceChunks)) #write while adding \n
		FileFasta.write('\n') #\n after the last, too
	else:
		FileFasta.write(sSequence)
		FileFasta.write('\n') #\n after each entry

	return(0)

##############################################################################
############################# fCompare2FastaFiles ############################
##############################################################################
def fCompare2FastaFiles(sDatabase1Path, sDatabase2Path, bOneContainedInTwo = False):

# 	my($numFragments1) = &countFragments($RARHfproteins1);
# 	print "\n<BR>num Entries with fragment in header: $numFragments1\n";

	#read in sequence databases to lists of dictionaries
	print('\n<BR>Processing Entries  in ', sDatabase1Path, '....')
	LDfproteins1 = fReadFastaFile(sDatabase1Path, False)
	#print('\n<BR>a total of ', str(len(LDfproteins1)), ' entries read.')

	print('\n<BR>Processing Entries  in ', sDatabase2Path, '....')
	LDfproteins2 = fReadFastaFile(sDatabase2Path, False)
	#print('\n<BR>a total of ', str(len(LDfproteins2)), ' entries read.')

# 	my($numFragments2) = &countFragments($RARHfproteins2);
# 	print "\n<BR>num Entries with fragment in header: $numFragments2\n";

	#sort in descending order of sequenceLength
	LDfproteins1  = sorted(LDfproteins1 , key = lambda i: (i['sequenceLength']), reverse=True) #sort descending by sequence length
	LDfproteins2  = sorted(LDfproteins2 , key = lambda i: (i['sequenceLength']), reverse=True) #sort descending by sequence length

	print('<BR>Sorted proteins by sequence length.')

	#create output files
	FileUpdates = open(sDatabase2Path + '.updates.txt', 'w')
	FileDeletes2 = open(sDatabase2Path + '.deletes.txt', 'w')
	FileDeletes1 = open(sDatabase1Path + '.deletes.txt', 'w')

	sCrossReferencesFile = sDatabase2Path + '.crossRefs.' + os.path.basename(sDatabase1Path) + '.categories.tsv'	#create a categories file
	FileCrossReferences = open(sCrossReferencesFile, 'w')
	sCrossReferencesMisMatchedGSFile =sDatabase2Path + '.crossRefs.GSmismatched.' + os.path.basename(sDatabase1Path) + '.txt'	#create a categories files
	FileCrossReferencesMisMatchedGS = open(sCrossReferencesMisMatchedGSFile, 'w')

	#find exact matches
	fMatchProteinSequences(LDfproteins1, LDfproteins2, os.path.basename(sDatabase1Path), os.path.basename(sDatabase2Path), FileUpdates, FileDeletes2, FileDeletes1, FileCrossReferences, FileCrossReferencesMisMatchedGS, sCrossReferencesFile, bOneContainedInTwo = bOneContainedInTwo)

	#find contained
	#FileMatches = open(sProtFASTAdir + sFastaFile + '.matches.txt', 'w')
	#LDfproteins = fMatchProteinFragmentsLeadingTrailingSpeedyTri(LDfproteins, sFastaFile, FileMatches);
	#FileMatches.close()

	FileUpdates.close()
	FileDeletes2.close()
	FileDeletes1.close()
	FileCrossReferences.close()
	FileCrossReferencesMisMatchedGS.close()
	
# 	&makeSequenceLengthHistograms(\@ARHsortedFproteins1, $database1)
# 	&makeSequenceLengthHistograms(\@ARHsortedFproteins2, $database2)

	return(0)

##############################################################################
############################ fMatchProteinSequences ##########################
##############################################################################

def fMatchProteinSequences(LDfproteins1, LDfproteins2, sDatabase1, sDatabase2, FileUpdates, FileDeletes2, FileDeletes1, FileCrossReferences, FileCrossReferencesMisMatchedGS, sCrossReferencesFile, bOneContainedInTwo = False):

	# $LDfproteins1, $LDfproteins2 expected to be sorted descending by length
	
	iNumProteins1 = len(LDfproteins1)
	iLongestLength1 = LDfproteins1[0]['sequenceLength']

	#some have a sequence 'n/a'
	for Dprotein1 in LDfproteins1:
		if Dprotein1['sequenceLength'] < 5:
			print ('\n<BR>', str(Dprotein1['sequenceLength']), Dprotein1['header'])
			print ('\n<BR>', Dprotein1['sequence'])
			break
		iShortestLength1 = Dprotein1['sequenceLength']

	print('<BR>DB1 longest protein: ', str(iLongestLength1), ' shortest protein: ', str(iShortestLength1))

	#find identical matches
	iNumMismatchedGeneSymbols = iNumMatchedGeneSymbols = iNumMissingGeneSymbols1 = iNumMissingGeneSymbols2 = 0
	DgeneSymbols1 = {}
	DgeneSymbols2 = {}
	for Dprotein1 in LDfproteins1:
		Dprotein1['geneSymbol'] = lFH.fParseGeneSymbolFromHeader(Dprotein1['header'])
		if Dprotein1['geneSymbol']:
			if Dprotein1['geneSymbol'] in DgeneSymbols1:
				DgeneSymbols1[Dprotein1['geneSymbol']] += 1
			else:
				DgeneSymbols1[Dprotein1['geneSymbol']] = 1
		else:
			iNumMissingGeneSymbols1 += 1

	FileDeletes2.write('entryCounter\taccession_number\tLength\tNterm 50 AA\theader\n')
	FileDeletes1.write('entryCounter\taccession_number\tLength\tNterm 50 AA\theader\n')
	FileCrossReferences.write('accession_number\tID_' + sDatabase1 + '\tgeneSymbol\tgeneSymbol_' + sDatabase1 + '\tSequenceMatchType\n')
	FileCrossReferencesMisMatchedGS.write('ID_' + sDatabase2 + '\tID_' + sDatabase1 + '\tgeneSymbol_' + sDatabase2 + '\tgeneSymbol_' + sDatabase1 + '\tmatchType' +'\n')

	#match the sequences
	print('bOneContainedInTwo: ', bOneContainedInTwo)
	iEntryCounter = 0
	iStartIndex = 0
	iNumIDmatches = 0
	iNumUpdated = 0
	iNumDeleted = 0
	iNumContainedMatches = 0
	iCurrentLongestLength = iLongestLength1
	for Dprotein2 in LDfproteins2:
		iEntryCounter += 1
		if iEntryCounter % 10000 == 0 :
			print('<BR>' , str(iEntryCounter), ' entries from ', sDatabase2 ,'compared, identical sequences: ', str(iNumIDmatches), ' contained sequences: ', str(iNumContainedMatches))
		if Dprotein2['sequenceLength'] > iCurrentLongestLength:
			continue
		if Dprotein2['sequenceLength'] < iShortestLength1:
			iStartIndex = iNumProteins1 + 1 #skip j loop below
			#break
		bMatched = False
		sAccession_number2 = lFH.fAccessionNumberFromHeader(Dprotein2['header'], sDatabase2)
		#sAccNum2_polished = sAccession_number2.split('-')[0]
		Dprotein2['geneSymbol'] = lFH.fParseGeneSymbolFromHeader(Dprotein2['header'])
		if Dprotein2['geneSymbol'] :
			if Dprotein2['geneSymbol'] in DgeneSymbols2:
				DgeneSymbols2[Dprotein2['geneSymbol']] += 1
			else:
				DgeneSymbols2[Dprotein2['geneSymbol']] = 1
		else:
			iNumMissingGeneSymbols2 += 1

		for j in range(iStartIndex, iNumProteins1):
			if LDfproteins1[j]['sequenceLength'] < Dprotein2['sequenceLength'] and bOneContainedInTwo == False:
				continue
			elif LDfproteins1[j]['sequenceLength'] > Dprotein2['sequenceLength']:
				#accelerate subsequent passes through this loop
				iStartIndex = j
				iCurrentLongestLength = LDfproteins1[j]['sequenceLength']

			#exact matches
			elif LDfproteins1[j]['sequenceLength'] == Dprotein2['sequenceLength'] and LDfproteins1[j]['sequence'] == Dprotein2['sequence']:
				iNumIDmatches += 1
				bMatched = True
				LDfproteins1[j]['dbMatch']= True
				sAccession_number1 = lFH.fAccessionNumberFromHeader(LDfproteins1[j]['header'], sDatabase1)
				#sAccNum1_polished = sAccession_number1.split('-')[0]
				#print "\n<BR>$entryCounter &nbsp $accession_number2 &nbsp $accession_number1";
				if sAccession_number2 != sAccession_number1: #count # with revised accession number
					iNumUpdated += 1
					#print "\n<BR>$database2: $entryCounter &nbsp $accession_number2 &nbsp sDatabase1: $j $accession_number1";
					#header has carriage return			
					FileUpdates.write(str(iEntryCounter) + '\t' + sAccession_number2 + '\t' + Dprotein2['header'] + '\n')
					FileUpdates.write('\t' + sAccession_number1 + '\t' + LDfproteins1[j]['header'] + '\n')
					#print $crossReferencesFileHandle "$accession_number2\t$accession_number1";
					#print $crossReferencesFileHandle "sAccNum2_polished\tsAccNum1_polished";
					#if($accNum1_polished eq ""){
					#	print $crossReferencesFileHandle "\t" . "$$LDfproteins1[$j]->{'header'}";
					#}else{
					#	print $crossReferencesFileHandle "\n";
					#}
				if Dprotein2['geneSymbol'] != LDfproteins1[j]['geneSymbol']: #count # with revised accession number
					iNumMismatchedGeneSymbols += 1
					FileCrossReferencesMisMatchedGS.write(sAccession_number2 + '\t' + sAccession_number1 + '\t' + Dprotein2['geneSymbol'] + '\t' + LDfproteins1[j]['geneSymbol'] + '\n')
				elif Dprotein2['geneSymbol']  and  LDfproteins1[j]['geneSymbol'] : #values != -1
					iNumMatchedGeneSymbols += 1
				FileCrossReferences.write(sAccession_number2 + '\t' + sAccession_number1 + '\t' + Dprotein2['geneSymbol'] + '\t' + LDfproteins1[j]['geneSymbol'] + '\tExact' + '\n')
				break #assumes each db is non-redundant
				##if there is redundancy, then db1 vs db2 will have different $numIDmatches than db2 vs db1
			
			#DB1 entry is substring of DB2 entry
			elif bOneContainedInTwo== True and LDfproteins1[j]['sequence'] in Dprotein2['sequence']:
				iNumContainedMatches += 1
				bMatched = True 
				LDfproteins1[j]['dbMatch']= True
				sAccession_number1 = lFH.fAccessionNumberFromHeader(LDfproteins1[j]['header'], sDatabase1)
				if sAccession_number2 != sAccession_number1: #count # with revised accession number
					iNumUpdated += 1

					#header has carriage return			
					FileUpdates.write(str(iEntryCounter) + '\t' + sAccession_number2 + '\t' + Dprotein2['header'] + '\n')
					FileUpdates.write('\t' + sAccession_number1 + '\t' + LDfproteins1[j]['header'] + '\n')
				if Dprotein2['geneSymbol'] != LDfproteins1[j]['geneSymbol']: #count # with revised accession number
					iNumMismatchedGeneSymbols += 1
					FileCrossReferencesMisMatchedGS.write(sAccession_number2 + '\t' + sAccession_number1 + '\t' + Dprotein2['geneSymbol'] + '\t' + LDfproteins1[j]['geneSymbol'] + '\n')
				elif Dprotein2['geneSymbol']  and  LDfproteins1[j]['geneSymbol'] : #values != -1
					iNumMatchedGeneSymbols += 1
				FileCrossReferences.write(sAccession_number2 + '\t' + sAccession_number1 + '\t' + Dprotein2['geneSymbol'] + '\t' + LDfproteins1[j]['geneSymbol'] + '\tContained' + '\n')
				break #assumes each db is non-redundant
				##if there is redundancy, then db1 vs db2 will have different $numIDmatches than db2 vs db1

		# end $j loop
		if not bMatched: # count # unmatched entries
			iNumDeleted += 1
			sPrefix = Dprotein2['sequence'][ 0: 50] #get 1st 50 residues of sequence
			FileDeletes2.write( str(iEntryCounter) + '\t' + sAccession_number2 + '\t' + str(Dprotein2['sequenceLength']) + '\t' + sPrefix + '\t' + Dprotein2['header'] + '\n')

	#write the unmatched proteins from db 1
	iNumNewSequences1 = 0
	for Dprotein1 in LDfproteins1:
		if 'dbMatch' not in Dprotein1:
			iNumNewSequences1 += 1
			sPrefix = Dprotein1['sequence'][ 0: 50] #get 1st 50 residues of sequence
			sAccession_number1 = lFH.fAccessionNumberFromHeader(Dprotein1['header'], sDatabase1)
			FileDeletes1.write( str(iNumNewSequences1) + '\t' + sAccession_number1 + '\t' + str(Dprotein1['sequenceLength']) + '\t' + sPrefix + '\t' + Dprotein1['header'] + '\n')

	iNumAccUnchanged = iNumIDmatches - iNumUpdated
	#iNumNewSequences = iNumProteins1 - iNumIDmatches - iNumContainedMatches

	print('<BR># entries from ' + sDatabase1 + ' found in ' + sDatabase2 + ' with sequences as exact matches:\t' + str(iNumIDmatches))
	print('<BR># entries from ' + sDatabase1 + ' found in ' + sDatabase2 + ' with sequences as contained matches:\t' + str(iNumContainedMatches))
	print('<BR># matching sequences in both databases with different accession numbers :\t' + str(iNumUpdated))
	print('<BR># matching sequences in both databases with same accession numbers :\t' + str(iNumAccUnchanged))
	print('<BR># matching sequences in both databases with different GeneSymbols :\t' + str(iNumMismatchedGeneSymbols))
	print('<BR># matching sequences in both databases with same GeneSymbols :\t' + str(iNumMatchedGeneSymbols))
	print('<BR># sequences from ' + sDatabase2 + ' absent from ' + sDatabase1 + ':\t' + str(iNumDeleted))

	iNumProteins2 = len(LDfproteins2)
	dPercentDeleted2 = 100 * iNumDeleted / iNumProteins2;
	print('\t%0.1f'%dPercentDeleted2 + '%')

	dPercentNew1 = 100 * iNumNewSequences1 / len(LDfproteins1);
	print('<BR># sequences from ' + sDatabase1 + ' absent from ' + sDatabase2 + ':\t' + str(iNumNewSequences1))
	print('\t%0.1f'%dPercentNew1 + '%')

	print('<BR># entries missing GeneSymbols from ' + sDatabase1 + ':' + str(iNumMissingGeneSymbols1))
	print('<BR># entries missing GeneSymbols from ' + sDatabase2 + ':' + str(iNumMissingGeneSymbols2))
	#fcalculateGeneSymbolOverlap(DgeneSymbols1, DgeneSymbols2, sDatabase1, sDatabase2, sCrossReferencesFile)

	return (0)

##############################################################################
######################### fCalculateGeneSymbolOverlap ########################
##############################################################################
#KRC 6/2/2023 from millscripts/fastaManipulator.pl
def fCalculateGeneSymbolOverlap(DgeneSymbols1, DgeneSymbols2, sDatabase1, sDatabase2, sCrossReferencesFile):

# 	my(@AgeneSymbols1) = keys(%$RHgeneSymbols1);
# 	my(@AgeneSymbols2) = keys(%$RHgeneSymbols2);
# 	my($numGeneSymbols1) = scalar(@AgeneSymbols1);
# 	my($numGeneSymbols2) = scalar(@AgeneSymbols2);
# 	my($databaseSmall);
# 	print "\n<BR>\# GeneSymbols from $database1 :\t$numGeneSymbols1";
# 	print "\n<BR>\# GeneSymbols from $database2 :\t$numGeneSymbols2";
# 	if($numGeneSymbols1 > 0 && $numGeneSymbols2){
# 		my(@AuniqueGeneSymbols1);
# 		my(@AuniqueGeneSymbols2);

# 		my($numOverlappingGeneSymbols12)=0;
# 		foreach my $gs (@AgeneSymbols1){
# 			if(exists($RHgeneSymbols2->{$gs}) ){
# 				$numOverlappingGeneSymbols12++;
# 			}else{
# 				push (@AuniqueGeneSymbols1, $gs);
# 			}
# 		}
# 		print "\n<BR>\# GeneSymbols in both 1 and 2:\t$numOverlappingGeneSymbols12";
# 		
# 		&writeUniqueGeneSymbolsFiles(\@AuniqueGeneSymbols1, $database1, $crossReferencesFile);

# 		my($numOverlappingGeneSymbols21)=0;
# 		foreach my $gs (@AgeneSymbols2){
# 			if(exists($RHgeneSymbols1->{$gs}) ){
# 				$numOverlappingGeneSymbols21++;
# 			}else{
# 				push (@AuniqueGeneSymbols2, $gs);
# 			}
# 		}
# 		print "\n<BR>\# GeneSymbols in both 2 and 1:\t$numOverlappingGeneSymbols21";
# 		
# 		&writeUniqueGeneSymbolsFiles(\@AuniqueGeneSymbols2, $database2, $crossReferencesFile);

	return (0)
##############################################################################
####################### fMatchProteinSequencesRedundant ######################
##############################################################################
# Make a non-redundant fasta of the input protein list, along with a list of the removed accession numbers matched to their retained counterparts
def fMatchProteinSequencesRedundant(LDfproteins, sDatabaseFilename, sMatchesFile = ''):

	iNumProteins = len(LDfproteins)
	
	if('accession_number' in LDfproteins[0]):
		LDfproteins  = sorted(LDfproteins , key = lambda i: (i['accession_number'])) #sort ascending by accession_number

	#enforce some database specific preferences for which one to keep
	if(sDatabaseFilename.startswith('UniProt') or sDatabaseFilename.startswith('RefSeq')):
		if(sDatabaseFilename.startswith('UniProt')):
			for Dfprotein in LDfproteins:
				if(Dfprotein['header'].startswith('>sp')):
					Dfprotein['db_source'] = "sp"
				elif(Dfprotein['header'].startswith('>tr')):
					Dfprotein['db_source'] = "tr"				
		elif(sDatabaseFilename.startswith('RefSeq')):
			##NP accession prefixes are preferred
			for Dfprotein in LDfproteins:
				if(Dfprotein['header'].startswith('>NP_')):
					Dfprotein['db_source'] = "NP_"		
		
		## sort descending by length with sp entries 1st
		LDfproteins  = sorted(LDfproteins , key = lambda i: (i['db_source'])) #sort ascending by 'db_source'

	LDfproteins  = sorted(LDfproteins , key = lambda i: (i['sequenceLength']), reverse=True) #sort descending by sequence length
	
	iNumCompleteMatches = 0
	iStartIndex = 0
	iNum_sp_tr_reundundant_pair = 0
	iNumTotalMultipleRedundancies = 0
	LmatchesStrings = []
	for Dfprotein in LDfproteins:
		iMultipleRedundancies = 0
		iStartIndex += 1
		if ( iStartIndex % 20000 == 0 ):	print (iStartIndex, " entries compared")
		#skip over entries already marked as redundant
		if('dbMatch' in Dfprotein and Dfprotein['dbMatch'] == 'R'):	continue
		else:								Dfprotein['dbMatch'] = 0	#set this so all others will be filtered out in fMakeDatabaseUnmatched
		
		#$accession_number_long = &accesionNumberFromHeader(Dfprotein['header'], sDatabaseFilename)	
		#sAccession_number_long = Dfprotein['accession_number']
		
		#for( my $j = $startIndex $j < $numProteins $j++):
		for j in range(iStartIndex, iNumProteins):
			##skip ones that are too short
			if( LDfproteins[j]['sequenceLength'] < Dfprotein['sequenceLength']):
				break
			
			if(LDfproteins[j]['sequence'] == Dfprotein['sequence']):
				iNumCompleteMatches += 1
				LDfproteins[j]['dbMatch'] = "R"	#mark as redundant
				iMultipleRedundancies += 1
				
				#keep redundant_accession_numbers assoicated with the kept one for output into a categories file
				if('redundant_accession_numbers' in Dfprotein):
					Dfprotein['redundant_accession_numbers'] += '|' + str(LDfproteins[j]['accession_number'])
				else:
					Dfprotein['redundant_accession_numbers'] = str(LDfproteins[j]['accession_number'])
				
				##determine if this redundant pair has 1 sp and 1 tr entry
				iMatched_sp_tr_reundundant_pair = 0
				
				#enforce some database specific preferences for which one to keep
				if(sDatabaseFilename.startswith('UniProt')):
					if(Dfprotein['header'].startswith('>sp')):
						if(LDfproteins[j]['header'].startswith('>tr')):
							iNum_sp_tr_reundundant_pair += 1
							iMatched_sp_tr_reundundant_pair = 1
						
					elif(Dfprotein['header'].startswith('>tr')):
						if(LDfproteins[j]['header'].startswith('>sp')):
							iNum_sp_tr_reundundant_pair += 1
							iMatched_sp_tr_reundundant_pair = 1
							print('<br> Entries not sorted with SP before TR')
				
				# append to a list of the removed accession numbers matched to their retained counterparts
				if not(iMatched_sp_tr_reundundant_pair):

					#include transcriptID and gene_symbol if present
					keep_transcriptID = ''
					redundant_transcriptID = ''
					keep_gene_symbol = ''
					redundant_gene_symbol = ''
					if('transcriptID' in Dfprotein): keep_transcriptID				= Dfprotein['transcriptID']
					if('transcriptID' in LDfproteins[j]): redundant_transcriptID	= LDfproteins[j]['transcriptID']
					if('gene_symbol' in Dfprotein): keep_gene_symbol				= Dfprotein['gene_symbol']
					if('gene_symbol' in LDfproteins[j]): redundant_gene_symbol		= LDfproteins[j]['gene_symbol']
					LmatchesStrings.append(str(Dfprotein['accession_number']) + '\t' + str(keep_transcriptID) + '\t' + str(keep_gene_symbol) + '\t' + str(Dfprotein['sequenceLength']) + '\t' + str(LDfproteins[j]['accession_number']) + '\t' + redundant_transcriptID + '\t' + redundant_gene_symbol + '\n')

				continue	
		# end $j loop
		if(iMultipleRedundancies > 1):
			#print "\n<BR>$multipleRedundancies: Dfprotein['header']"
			iNumTotalMultipleRedundancies += 1	

	#KRC 5/26/2024 only write matches file if something to write
	if len(LmatchesStrings):
		FileMatches = open(sMatchesFile, 'w')
		##write headers for output
		FileMatches.write('accession_number\ttranscript_id\tgene_symbol\tlength\tredundant_accession_number\tredundant_transcript_id\tredundant_gene_symbol\n')
		for sMatches in LmatchesStrings:
			FileMatches.write(sMatches) 
		FileMatches.close()
	
	print('<br># of entries with redundant sequences: ', iNumCompleteMatches)	
	if(sDatabaseFilename.startswith('UniProt')):
		print('<br># of redundancies with sp tr pair: ', iNum_sp_tr_reundundant_pair)	
	
	print('<br># of entries with multiple redundancies: ', iNumTotalMultipleRedundancies)
	
	return (LDfproteins)

##############################################################################
############### fMatchProteinFragmentsLeadingTrailingSpeedyTri ###############
##############################################################################
def fMatchProteinFragmentsLeadingTrailingSpeedyTri(LDfproteins, sDatabaseFile, FileMatches):

	inumProteins = len(LDfproteins)
	
	#Divide the total batch of proteins into 2 sets of ~ 8000 batches according to the tripeptides found at the N and C terminal
	#then each protein only has to undergo 2 x 1/8000 as many string compare operations as it would if the subsets were not created
	# in order to find shorter proteins contained entirely within a longer protein

	#initialize Tripeptide accelerating objects
	DLDleadingTripeptideArrays = {}
	DLDtrailingTripeptideArrays = {}
	DleadingTripeptideArrayStartIndices = {}
	DtrailingTripeptideArrayStartIndices = {}
	#allow for any alphabetic character eventhough the database should have only the std 20 and X, there may be strays
	Lalphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
	for sAminoAcid1 in Lalphabet:
		for sAminoAcid2 in Lalphabet:
			for sAminoAcid3 in Lalphabet:
				sTripeptide = sAminoAcid1 + sAminoAcid2 + sAminoAcid3
				DLDleadingTripeptideArrays[sTripeptide] = []
				DLDtrailingTripeptideArrays[sTripeptide] = []
				DleadingTripeptideArrayStartIndices[sTripeptide] = 0
				DtrailingTripeptideArrayStartIndices[sTripeptide] = 0

	LDleadingTripeptideArray = []
	LDtrailingTripeptideArray = []
	#fill protein object with Tripeptide accelerating info
	for Dprotein in LDfproteins:
		#skip previously matched ones
		if 'dbMatch' in Dprotein:
			continue
		Dprotein['dbMatch'] = 0

		Dprotein['leadingTripeptide']	= Dprotein['sequence'][1:4] #get residues 2-4 of sequence skip 1 it is usually Met
		LDleadingTripeptideArray = DLDleadingTripeptideArrays[Dprotein['leadingTripeptide']]
		LDleadingTripeptideArray.append(Dprotein )

		Dprotein['trailingTripeptide']	= Dprotein['sequence'][-4:-1]  #get last 3 residues of sequence
		LDtrailingTripeptideArray = DLDtrailingTripeptideArrays[Dprotein['trailingTripeptide']]
		LDtrailingTripeptideArray.append(Dprotein)

	#write headers for output
	FileMatches.write("Matched_terminal\tnum_matches\tlengthAA\taccession_number\tlengthType\tNterm_10AA\theader\n")

	iEntryCounter = iNumCtermMatches = iNumNtermMatches = iNumPreviousMatches = iNumCompleteMatches = 0
	for Dprotein in LDfproteins:
		iEntryCounter+=1
		if iEntryCounter % 10000 == 0 :
			print(iEntryCounter, ' entries compared')
			print("\t\tnumCompleteMatches: ", iNumCompleteMatches, 'iNumCtermMatches: ', iNumCtermMatches, ' numNtermMatches: ', iNumNtermMatches, ' numPreviousMatches: ', iNumPreviousMatches)

		#skip previously matched ones
		if Dprotein['dbMatch']:
			iNumPreviousMatches+=1
			continue

		sAccession_number_long = lFH.fAccessionNumberFromHeader(Dprotein['header'], sDatabaseFile)

		#Retrieve the subset of proteins with a N-terminal Tripeptide match to examine for full N-terminal exact matches
		LDtrailingTripeptideArray = DLDtrailingTripeptideArrays[Dprotein['trailingTripeptide']]
		iNumProteins = len(LDtrailingTripeptideArray)
		iStartIndex = DtrailingTripeptideArrayStartIndices[Dprotein['trailingTripeptide']]

		#match C-term
		for j in range(iStartIndex, iNumProteins):
			#skip ones that are too long
			if LDtrailingTripeptideArray[j]['sequenceLength'] > Dprotein['sequenceLength']:
				DtrailingTripeptideArrayStartIndices[Dprotein['trailingTripeptide']] = j + 1
				continue
			if LDtrailingTripeptideArray[j]['sequenceLength'] == Dprotein['sequenceLength']:
				DtrailingTripeptideArrayStartIndices[Dprotein['trailingTripeptide']] = j
				##count complete matches
				if j != (iEntryCounter - 1) and LDfproteins[j]['sequence'] == Dprotein['sequence']:
					iNumCompleteMatches+=1
				continue

			#skip previously matched ones
			if LDtrailingTripeptideArray[j]['dbMatch']:
				continue
	
			if Dprotein['sequence'].endswith(LDtrailingTripeptideArray[j]['sequence']):
				iNumCtermMatches+=1
				LDtrailingTripeptideArray[j]['dbMatch'] = 'C'
				sAccession_number_short = lFH.fAccessionNumberFromHeader(LDtrailingTripeptideArray[j]['header'], sDatabaseFile)
#				print "<P>C-term fragment $numCtermMatches long: $RHprotein->['header']";
#				print "<BR>&nbsp&nbsp&nbsp short: $$RARHtrailingTripeptideArray[$j]->['header']";
				#sPrefix10 = Dprotein['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = Dprotein['sequence'] #all sequence
				FileMatches.write('C\t' + str(iNumCtermMatches) + '\t' + str(Dprotein['sequenceLength']) + '\t' + sAccession_number_long + '\tlong\t' + sPrefix10 + '\t' + Dprotein['header'] + '\n')
				#sPrefix10 = LDtrailingTripeptideArray[j]['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = LDtrailingTripeptideArray[j]['sequence'] #all sequence
				FileMatches.write('C\t' + str(iNumCtermMatches) +'\t' + str(LDtrailingTripeptideArray[j]['sequenceLength']) + '\t' +  sAccession_number_short + '\tshort\t' + sPrefix10 + '\t' + LDtrailingTripeptideArray[j]['header'] + '\n')

		#Retrieve the subset of proteins with an N-terminal Tripeptide match to examine for full N-terminal exact matches
		LDleadingTripeptideArray = DLDleadingTripeptideArrays[Dprotein['leadingTripeptide']]
		iNumProteins = len(LDleadingTripeptideArray)
		iStartIndex = DleadingTripeptideArrayStartIndices[Dprotein['leadingTripeptide']]

		#match N-term
		for j in range(iStartIndex, iNumProteins):
			#skip ones that are too long
			if LDleadingTripeptideArray[j]['sequenceLength'] > Dprotein['sequenceLength']:
				DleadingTripeptideArrayStartIndices[Dprotein['leadingTripeptide']] = j + 1
				continue
			if LDleadingTripeptideArray[j]['sequenceLength'] == Dprotein['sequenceLength']:
				DleadingTripeptideArrayStartIndices[Dprotein['leadingTripeptide']] = j
				continue

			#skip previously matched ones
			if LDleadingTripeptideArray[j]['dbMatch']:
				continue
			
			if Dprotein['sequence'].startswith(LDleadingTripeptideArray[j]['sequence']):
				iNumNtermMatches+=1
				LDleadingTripeptideArray[j]['dbMatch'] = 'N'
				sAccession_number_short = lFH.fAccessionNumberFromHeader(LDleadingTripeptideArray[j]['header'], sDatabaseFile)
#				print "<P>N-term fragment $numNtermMatches long: $RHprotein->{'header'}";
#				print "<BR>&nbsp&nbsp&nbsp short: $$RARHleadingTripeptideArray[$j]->{'header'}";
				#sPrefix10 = Dprotein['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = Dprotein['sequence'] #all sequence
				FileMatches.write('N\t' + str(iNumNtermMatches) + '\t' + str(Dprotein['sequenceLength']) + '\t' + sAccession_number_long + '\tlong\t' + sPrefix10 + '\t' + Dprotein['header'] + '\n')
				#sPrefix10 = LDleadingTripeptideArray[j]['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = LDleadingTripeptideArray[j]['sequence'] #all sequence
				FileMatches.write('N\t' + str(iNumNtermMatches) +'\t' + str(LDleadingTripeptideArray[j]['sequenceLength']) + '\t' + sAccession_number_short + '\tshort\t' + sPrefix10 + '\t' + LDleadingTripeptideArray[j]['header'] + '\n')

	print('numCompleteMatches: ',	iNumCompleteMatches)
	print('numCtermMatches ',		iNumCtermMatches)
	print('numNtermMatches: ',		iNumNtermMatches)

	return (LDfproteins)

##############################################################################
############### fMatchProteinFragmentsLeadingTrailingSpeedyTri2 ###############
##############################################################################
#unfinished KRC 10/9/2025
def fMatchProteinFragmentsLeadingTrailingSpeedyTri2(LDfproteins1, LDfproteins2, sDatabaseFile, FileMatches):

	inumProteins = len(LDfproteins)
	
	#Divide the total batch of proteins into 2 sets of ~ 8000 batches according to the tripeptides found at the N and C terminal
	#then each protein only has to undergo 2 x 1/8000 as many string compare operations as it would if the subsets were not created
	# in order to find shorter proteins contained entirely within a longer protein

	#initialize Tripeptide accelerating objects
	DLDleadingTripeptideArrays = {}
	DLDtrailingTripeptideArrays = {}
	DleadingTripeptideArrayStartIndices = {}
	DtrailingTripeptideArrayStartIndices = {}
	#allow for any alphabetic character eventhough the database should have only the std 20 and X, there may be strays
	Lalphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
	for sAminoAcid1 in Lalphabet:
		for sAminoAcid2 in Lalphabet:
			for sAminoAcid3 in Lalphabet:
				sTripeptide = sAminoAcid1 + sAminoAcid2 + sAminoAcid3
				DLDleadingTripeptideArrays[sTripeptide] = []
				DLDtrailingTripeptideArrays[sTripeptide] = []
				DleadingTripeptideArrayStartIndices[sTripeptide] = 0
				DtrailingTripeptideArrayStartIndices[sTripeptide] = 0

	#fill protein object in DB2 with Tripeptide accelerating info
	LDleadingTripeptideArray = []
	LDtrailingTripeptideArray = []
	for Dprotein in LDfproteins2:
		#skip previously matched ones
		if 'dbMatch' in Dprotein:
			continue
		Dprotein['dbMatch'] = 0

		Dprotein['leadingTripeptide']	= Dprotein['sequence'][1:4] #get residues 2-4 of sequence skip 1 it is usually Met
		LDleadingTripeptideArray = DLDleadingTripeptideArrays[Dprotein['leadingTripeptide']]
		LDleadingTripeptideArray.append(Dprotein )

		Dprotein['trailingTripeptide']	= Dprotein['sequence'][-4:-1]  #get last 3 residues of sequence
		LDtrailingTripeptideArray = DLDtrailingTripeptideArrays[Dprotein['trailingTripeptide']]
		LDtrailingTripeptideArray.append(Dprotein)

	#write headers for output
	FileMatches.write("Matched_terminal\tnum_matches\tlengthAA\taccession_number\tlengthType\tNterm_10AA\theader\n")

	iEntryCounter = iNumCtermMatches = iNumNtermMatches = iNumPreviousMatches = iNumCompleteMatches = 0
	for Dprotein in LDfproteins:
		iEntryCounter+=1
		if iEntryCounter % 10000 == 0 :
			print(iEntryCounter, ' entries compared')
			print("\t\tnumCompleteMatches: ", iNumCompleteMatches, 'iNumCtermMatches: ', iNumCtermMatches, ' numNtermMatches: ', iNumNtermMatches, ' numPreviousMatches: ', iNumPreviousMatches)

		#skip previously matched ones
		if Dprotein['dbMatch']:
			iNumPreviousMatches+=1
			continue

		sAccession_number_long = lFH.fAccessionNumberFromHeader(Dprotein['header'], sDatabaseFile)

		#Retrieve the subset of proteins with a N-terminal Tripeptide match to examine for full N-terminal exact matches
		LDtrailingTripeptideArray = DLDtrailingTripeptideArrays[Dprotein['trailingTripeptide']]
		iNumProteins = len(LDtrailingTripeptideArray)
		iStartIndex = DtrailingTripeptideArrayStartIndices[Dprotein['trailingTripeptide']]

		#match C-term
		for j in range(iStartIndex, iNumProteins):
			#skip ones that are too long
			if LDtrailingTripeptideArray[j]['sequenceLength'] > Dprotein['sequenceLength']:
				DtrailingTripeptideArrayStartIndices[Dprotein['trailingTripeptide']] = j + 1
				continue
			if LDtrailingTripeptideArray[j]['sequenceLength'] == Dprotein['sequenceLength']:
				DtrailingTripeptideArrayStartIndices[Dprotein['trailingTripeptide']] = j
				##count complete matches
				if j != (iEntryCounter - 1) and LDfproteins[j]['sequence'] == Dprotein['sequence']:
					iNumCompleteMatches+=1
				continue

			#skip previously matched ones
			if LDtrailingTripeptideArray[j]['dbMatch']:
				continue
	
			if Dprotein['sequence'].endswith(LDtrailingTripeptideArray[j]['sequence']):
				iNumCtermMatches+=1
				LDtrailingTripeptideArray[j]['dbMatch'] = 'C'
				sAccession_number_short = lFH.fAccessionNumberFromHeader(LDtrailingTripeptideArray[j]['header'], sDatabaseFile)
#				print "<P>C-term fragment $numCtermMatches long: $RHprotein->['header']";
#				print "<BR>&nbsp&nbsp&nbsp short: $$RARHtrailingTripeptideArray[$j]->['header']";
				#sPrefix10 = Dprotein['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = Dprotein['sequence'] #all sequence
				FileMatches.write('C\t' + str(iNumCtermMatches) + '\t' + str(Dprotein['sequenceLength']) + '\t' + sAccession_number_long + '\tlong\t' + sPrefix10 + '\t' + Dprotein['header'] + '\n')
				#sPrefix10 = LDtrailingTripeptideArray[j]['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = LDtrailingTripeptideArray[j]['sequence'] #all sequence
				FileMatches.write('C\t' + str(iNumCtermMatches) +'\t' + str(LDtrailingTripeptideArray[j]['sequenceLength']) + '\t' +  sAccession_number_short + '\tshort\t' + sPrefix10 + '\t' + LDtrailingTripeptideArray[j]['header'] + '\n')

		#Retrieve the subset of proteins with an N-terminal Tripeptide match to examine for full N-terminal exact matches
		LDleadingTripeptideArray = DLDleadingTripeptideArrays[Dprotein['leadingTripeptide']]
		iNumProteins = len(LDleadingTripeptideArray)
		iStartIndex = DleadingTripeptideArrayStartIndices[Dprotein['leadingTripeptide']]

		#match N-term
		for j in range(iStartIndex, iNumProteins):
			#skip ones that are too long
			if LDleadingTripeptideArray[j]['sequenceLength'] > Dprotein['sequenceLength']:
				DleadingTripeptideArrayStartIndices[Dprotein['leadingTripeptide']] = j + 1
				continue
			if LDleadingTripeptideArray[j]['sequenceLength'] == Dprotein['sequenceLength']:
				DleadingTripeptideArrayStartIndices[Dprotein['leadingTripeptide']] = j
				continue

			#skip previously matched ones
			if LDleadingTripeptideArray[j]['dbMatch']:
				continue
			
			if Dprotein['sequence'].startswith(LDleadingTripeptideArray[j]['sequence']):
				iNumNtermMatches+=1
				LDleadingTripeptideArray[j]['dbMatch'] = 'N'
				sAccession_number_short = lFH.fAccessionNumberFromHeader(LDleadingTripeptideArray[j]['header'], sDatabaseFile)
#				print "<P>N-term fragment $numNtermMatches long: $RHprotein->{'header'}";
#				print "<BR>&nbsp&nbsp&nbsp short: $$RARHleadingTripeptideArray[$j]->{'header'}";
				#sPrefix10 = Dprotein['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = Dprotein['sequence'] #all sequence
				FileMatches.write('N\t' + str(iNumNtermMatches) + '\t' + str(Dprotein['sequenceLength']) + '\t' + sAccession_number_long + '\tlong\t' + sPrefix10 + '\t' + Dprotein['header'] + '\n')
				#sPrefix10 = LDleadingTripeptideArray[j]['sequence'][0:10] #get 1st 10 residues of sequence
				sPrefix10 = LDleadingTripeptideArray[j]['sequence'] #all sequence
				FileMatches.write('N\t' + str(iNumNtermMatches) +'\t' + str(LDleadingTripeptideArray[j]['sequenceLength']) + '\t' + sAccession_number_short + '\tshort\t' + sPrefix10 + '\t' + LDleadingTripeptideArray[j]['header'] + '\n')

	print('numCompleteMatches: ',	iNumCompleteMatches)
	print('numCtermMatches ',		iNumCtermMatches)
	print('numNtermMatches: ',		iNumNtermMatches)

	return (LDfproteins)

##############################################################################
############################# fMakeTSVfromMatched ############################
##############################################################################
def fMakeTSVfromMatched(LDfproteins, sDatabasePath, sDatabaseFile):

	if len(LDfproteins) == 0:
		return (0)

	FileTSV = open(sDatabasePath + '.deletes.txt', 'w')

	##write headers for output
	if 'fragment' in LDfproteins[0]:
		FileTSV.write("accession_number\tlength\tDeleteReason\tFragment\tprefix10\tHeader\n")
	else:
		FileTSV.write("accession_number\tlength\tRedundantGeneSymbol\tHeader\n")
		

	#write each entry that was not found to be matched
	for Dprotein in LDfproteins:
		if 'dbMatch' in Dprotein and Dprotein['dbMatch'] != 0:
			if 'accession_number' not in Dprotein:
				sAccession_number = lFH.fAccesionNumberFromHeader(Dprotein['header'], sDatabaseFile)
			else:
				sAccession_number = Dprotein['accession_number']
			
			if 'fragment' in LDfproteins[0]:
				sPrefix10 = Dprotein['sequence'][0:10] #get 1st 10 residues of sequence
				
				sFragment = "N";
				if '(Fragment)' in Dprotein['header']:
					sFragment = "Y";
				FileTSV.write(sAccession_number + '\t' + str(Dprotein['sequenceLength']) + '\t' + Dprotein['dbMatch'] + '\t' + sFragment + '\t' + sPrefix10 + '\t' + Dprotein['header'] + '\n')
			else:
				FileTSV.write(sAccession_number + '\t' + str(Dprotein['sequenceLength']) + '\t' + Dprotein['dbMatch'] + '\t' + Dprotein['header'] + '\n')


	FileTSV.close()

	return (0)

##############################################################################
################# fMatchProteinSequencesRedundantKeepPreferred ###############
##############################################################################
def fMatchProteinSequencesRedundantKeepPreferred(LDfproteins, DaccPreferred, sMatchesFile = ''):

	##starting sort by sequence length
	LDfproteins  = sorted(LDfproteins , key = lambda i: (i['sequenceLength']), reverse=True) #sort descending by sequence length
	
	#print('\n<br>Num Preferred protein acc: ', len(DaccPreferred), ' ',DaccPreferred.keys() )
	print('\n<br>Num Preferred protein acc: ', len(DaccPreferred), ' ' )
	iNumCompleteMatches = 0
	iNumPreferredKept = 0
	iStartIndex = 0
	iNumTotalMultipleRedundancies = 0
	iNumProteins = len(LDfproteins)
	LmatchesStrings = []
	for Dfprotein in LDfproteins:
		iMultipleRedundancies = 0
		if ( iStartIndex % 20000== 0 ):
			print('\n<br>', iStartIndex,  ' entries compared')
		iStartIndex+=1
		
		#skip over entries already marked as redundant
		if 'dbMatch' in Dfprotein and Dfprotein['dbMatch'] == "R":
			continue
		else:
			Dfprotein['dbMatch'] = 0
		
		bPreferDfprotein = False
		#if(exists($RHaccPreferred->{Dfprotein['accession_number'])):
		if Dfprotein['accession_number'] in DaccPreferred:
			bPreferDfprotein = True
			if not 'preferredKept' in Dfprotein:
				Dfprotein['preferredKept'] = 1
				iNumPreferredKept+=1

		#loop through proteins of same length to find redundant members
		for j in range(iStartIndex, iNumProteins):
			##skip ones that are too short
			if( LDfproteins[j]['sequenceLength'] < Dfprotein['sequenceLength']):
				break
			
			if LDfproteins[j]['sequence'] == Dfprotein['sequence']:
				iNumCompleteMatches+=1
				iMultipleRedundancies+=1
				#check for preferred status, eliminate those lacking it
				
				if LDfproteins[j]['accession_number'] in DaccPreferred:
					#print('\nkeeping preferred ', LDfproteins[j]['accession_number'])
					LDfproteins[j]['preferredKept'] = 1
					iNumPreferredKept+=1
					if not bPreferDfprotein:
						Dfprotein['dbMatch'] = 'R'	#mark as redundant
						LmatchesStrings.append('Preferred\t' + LDfproteins[j]['accession_number'] + '\t' + str(LDfproteins[j]['sequenceLength']) + '\t' + Dfprotein['accession_number'] + '\t' + str(Dfprotein['sequenceLength']) + '\n')
						LmatchesStrings.append(Dfprotein['accession_number'] + '\t' + str(Dfprotein['sequenceLength']) + '\t' + LDfproteins[j]['accession_number'] + '\t' + str(LDfproteins[j]['sequenceLength']) + '\n')
						continue #try more j's
					break #exit j loop and move on to next protein
				else:
					LDfproteins[j]['dbMatch'] = 'R'	#mark as redundant
					LmatchesStrings.append('First\t' + Dfprotein['accession_number'] +'\t' + str(Dfprotein['sequenceLength']) + '\t' + LDfproteins[j]['accession_number'] + '\t' + str(LDfproteins[j]['sequenceLength']) + '\n')
					continue
			
		# end $j loop

		if(iMultipleRedundancies > 1):
			#print "\n<BR>$multipleRedundancies: Dfprotein['header']"
			iNumTotalMultipleRedundancies+=1

	#KRC 5/26/2024 only write matches file if something to write
	if len(LmatchesStrings):
		FileMatches = open(sMatchesFile, 'w')
		##write headers for output
		FileMatches.write("KeepType\tAccession_number\tlength\tRedundant Accession Number\tlength\n")
		for sMatches in LmatchesStrings:
			FileMatches.write(sMatches) 
		FileMatches.close()

	print('\n<br># of entries with redundant sequences: ',		iNumCompleteMatches)
	print('\n<br># of entries with multiple redundancies: ',	iNumTotalMultipleRedundancies)
	print('\n<br># of preferred entries kept: ',				iNumPreferredKept)


	return (LDfproteins)

##############################################################################
########################### fRetainEntriesUnmatched ##########################
##############################################################################
def fRetainEntriesUnmatched(LDfproteins, FileLog = sys.stdout):
  
	LDfproteinsRetained = []	#initialize empty

	##retain each entry that was not found to be matched
	iNumProteinsUnFiltered = len(LDfproteins)
	iNumProteinsFiltered = 0
	for Dfprotein in LDfproteins:
		if ('dbMatch' in Dfprotein.keys() and Dfprotein['dbMatch'] == 0): #need eq instead of == because non-zero values are strings not numbers
			LDfproteinsRetained.append(Dfprotein)
			iNumProteinsFiltered += 1			
	
	iPercentRetained = 100
	if(iNumProteinsUnFiltered > 0):
		iPercentRetained = int(100 * iNumProteinsFiltered / iNumProteinsUnFiltered)
	
	FileLog.write('<br>Num proteins retained: ' + str(iNumProteinsFiltered) + ' (' "{0:0.1f}".format(iPercentRetained,)+ '%)')

	return (LDfproteinsRetained)

##############################################################################
#################### filterOutProteinsByMinSequenceLength ####################
##############################################################################
def filterOutProteinsByMinSequenceLength(LDfproteins, iMinLength):

	LDfproteins  = sorted(LDfproteins , key = lambda i: (i['sequenceLength']), reverse=True) #sort descending by sequence length
	iNumProteinsBefore = len(LDfproteins)

	LDfproteinsFiltered = []
	for Dfprotein in LDfproteins:
		if(Dfprotein['sequenceLength'] < iMinLength):	break	
		LDfproteinsFiltered.append(Dfprotein)
	
	iNumProteinsAfter = len(LDfproteinsFiltered)
	iNumProteinsDeleted = iNumProteinsBefore - iNumProteinsAfter
	print('<br># Proteins before length filtering: ', iNumProteinsBefore)
	print('<br># Proteins after length filtering: ', iNumProteinsAfter)
	print('<br># Proteins deleted < length ', iMinLength, ':', iNumProteinsDeleted)

	return (LDfproteinsFiltered)

##############################################################################
########################## filterByMinSequenceLength #########################
##############################################################################
def filterByMinSequenceLength(LDfSequences, iMinLength):

	LDfSequences  = sorted(LDfSequences , key = lambda i: (i['sequenceLength']), reverse=True) #sort descending by sequence length
	iNumSequencesBefore = len(LDfSequences)

	LDfSequencesRetained = []
	LremovedAcc = []
	for DfSequence in LDfSequences:
		if(DfSequence['sequenceLength'] < iMinLength): #KRC 6/18/2023 potential to accelerate this step by just taking all remaining
			LremovedAcc.append(DfSequence['accession_number'])
		else:
			LDfSequencesRetained.append(DfSequence)
	
	iNumSequencesAfter = len(LDfSequencesRetained)
	iNumSequencesDeleted = iNumSequencesBefore - iNumSequencesAfter
	print('<br># Sequences before length filtering: ', iNumSequencesBefore)
	print('<br># Sequences after length filtering: ', iNumSequencesAfter)
	print('<br># Sequences deleted < length ', iMinLength, ':', iNumSequencesDeleted)

	return (LDfSequencesRetained, LremovedAcc)

##############################################################################
################################ fFilterByAcc ################################
##############################################################################
def fFilterByAcc(LDfSequences, LremovedAcc):

	LDfSequencesRetained = []
	for DfSequence in LDfSequences:
		if DfSequence['accession_number'] in LremovedAcc:
			continue
		LDfSequencesRetained.append(DfSequence)

	return (LDfSequencesRetained)

##############################################################################
######################### fMakeFastaDatabaseUnmatched ########################
##############################################################################
def fMakeFastaDatabaseUnmatched(LDfproteins, sFilePath, sSuffix = ''):

	if not (sSuffix == ''):
		if(sFilePath.endswith('.fasta')):
			#$databasePath =~ s/\.fasta/$suffix\.fasta/g #insert $suffix before .fasta			#ToDo KRC use re. function
			sFilePath.replace('.fasta', sSuffix + '.fasta')
		else:								sFilePath += sSuffix

	FileFasta = open(sFilePath, 'w')

	#write each entry that was not found to be matched
	iNumProteinsFiltered = 0
	LDfproteinsNR=[]
	LDfproteinsRonly=[]
	for Dfprotein in LDfproteins:
		if (Dfprotein['dbMatch'] == 0): #need eq instead of == because non-zero values are strings not numbers
			fAppendEntryToOpenFastaFile(FileFasta, Dfprotein['header'], Dfprotein['sequence'])
			iNumProteinsFiltered += 1
			LDfproteinsNR.append(Dfprotein)
		else:
			LDfproteinsRonly.append(Dfprotein)
			
	FileFasta.close()
	print('<br>Finished writing file containing ', iNumProteinsFiltered, ' proteins: ', sFilePath)

	return (LDfproteinsNR, LDfproteinsRonly)
##############################################################################
############################  fMakeFastaNRdatabase ###########################
##############################################################################
# Make a non-redundant fasta of the input protein list, along with a list of the removed accession numbers matched to their retained counterparts
def fMakeFastaNRdatabase(LDfproteins, sFastaFilenameNR, sMatchesFile):

	LDfproteins = fMatchProteinSequencesRedundant(LDfproteins, sFastaFilenameNR, sMatchesFile)
	LDfproteinsNR, LDfproteinsRonly = fMakeFastaDatabaseUnmatched(LDfproteins, sFastaFilenameNR, '')

	return(LDfproteinsNR, LDfproteinsRonly) #returns ONLY the proteins with the marking of the attribute $RHfprotein->{'dbMatch'} = 0

##############################################################################
########################## fListOfDictionariesToDD ###########################
##############################################################################
def	fListOfDictionariesToDD(LDentries, sPrimaryKey, bTruncatePrimaryVersionID = False, sTruncateChar = '.'):

	DDentries = {}	#initialize empty dictionary
	for Dentry in LDentries:
		#if(sPrimaryKey in Dentry.keys()):
		if bTruncatePrimaryVersionID:
			if sPrimaryKey =='redundant_accession_numbers':
				if Dentry[sPrimaryKey] == '':
					continue
				LredAccc = Dentry[sPrimaryKey].split('|')
				for sRedAcc in LredAccc:
					DDentries[sRedAcc.split(sTruncateChar)[0]] = Dentry
				continue
			else:
				Dentry[sPrimaryKey] = Dentry[sPrimaryKey].split(sTruncateChar)[0]
			
		DDentries[Dentry[sPrimaryKey]] = Dentry

	print ('number xrefs configured for ', sPrimaryKey, ':', str(len(DDentries.keys())))

	return(DDentries)

#####################################################################
############## fFastaRetrieveSequenceViaAccessionNumber ############
#####################################################################
def fFastaRetrieveSequenceViaAccessionNumber(sDBpath, sAccNum):

	LaccNums = [] #initialize empty list
	if sDBpath == sInitializedDataBase: 	#return previously requested/retrieved info
		bVerboseProgressOff = True
		DDallProteinsByAccessionKey = fFastaFileToAccessionNumberOrientedDDprotein(sDBpath, LaccNums, bVerboseProgressOff) #relies on .idx, and .acn files for accession number

	return(DDallProteinsByAccessionKey[sAccNum]['sequence'])

#####################################################################
############ fFastaFileToAccessionNumberOrientedDDprotein ###########
#####################################################################
def fFastaFileToAccessionNumberOrientedDDprotein(sDBpath, LaccessionKeys, bVerboseProgressOff):

	# JR Note for >4.2Gb databases 03/11/09
	# faindex now creates '.idx2' files (__int64), and will only create 'idx' files if size < 4.2Gb
	# That's because unpack("Q") isn't supported with this version of Perl, so it must access idx files,
	# and I think seek() will fail also.
	# So, I moved the idx check here an it exits with error it doesn't exist
	#KRC 7/28/2021 now supporting only idx2 present. 64-bit PERL supports unpack("Q")
	bOffsets64Bit = False
	if os.path.isfile(sDBpath+ '.idx'):
		DDproteinsByAccessionKey = fFastaFileToAccessionNumberOrientedDDprotein64(sDBpath, LaccessionKeys, bVerboseProgressOff, bOffsets64Bit);
	elif os.path.isfile(sDBpath + '.idx2'):
		bOffsets64Bit = True
		DDproteinsByAccessionKey = fFastaFileToAccessionNumberOrientedDDprotein64(sDBpath, LaccessionKeys, bVerboseProgressOff, bOffsets64Bit);
	else:	# this really shouldn't happen unless someone deleted it
		print ('<br><b>Error: </b>The ' + sDBpath + '.idx2' + 'does not exist - you must index the database.')
		sys.exit(-1)

	return(DDproteinsByAccessionKey)

#####################################################################
########### fFastaFileToAccessionNumberOrientedDDprotein64 ##########
#####################################################################
#KRC 7/9/2023 - did not finish porting this from perl to python because found fReadFastaFileAcc already existed
def fFastaFileToAccessionNumberOrientedDDprotein64(sDBpath, LaccessionKeys, bVerboseProgressOff, bOffsets64Bit):

# 	my($RHRHproteinsByAccessionKey);

# 	if($db eq sInitializedDataBase){ 	#return previously requested/retrieved info
# 		$RHRHproteinsByAccessionKey = &parseFastaFileUsingIndicesToAccessionNumberOrientedRHRHprotein(sDBpath, LDaccessionKeys, $RAindexOffsets, $RHaccToIndex);
# 		return($RHRHproteinsByAccessionKey);

# 	# this call accesses entry offsets using .ac[cn] and .idx files
# 	# and only for those accKeys passed in (called by aligner.pl)
# 	# JR AGT 09/12/05

# 	sAcnFile = "$db.acn";	# try number-based file first
# 	my($isAlphaAccKey) = 0;
# 	if ( !(-e sAcnFile) ) {
# 		$acnFile = "$db.acc";		# assume string-based accession keys file
# 		$isAlphaAccKey = 1;
# 	}

# 	if not bVerboseProgressOff: print('\n<BR> Processing accession key file: $acnFile \n')
# 	# First read in all accessionKey <-> Index# maps from .ac[cn] file

# 	# while not eof, read in <acc#><index> pair
# 	# look at optimizing this.. for now just one field at a time
# 	# create a hash mapping acc# to index
# 	# TODO: Rethink how this might get handled in 64-bit windows (future)
# 	# For all SM3.2 and earlier acn and idx files, and int is 4 bytes
# 	my($sizeOfInt) = $Config{intsize};
# 	my($record);
# 	my($recordSize) = 2*$sizeOfInt;
	DaccToIndex= {}	#=&makeNewHash();
# 	my($entryCounter) = 0;
# 	if ($isAlphaAccKey) {
# 		my($acnLine) = "";
# 		open(ACN, $acnFile) or die "Can't open $acnFile : $!\n";
# 		# both values are text separated by space with \n terminator
# 		while ( $acnLine = <ACN> ) {
# 			chomp $acnLine;
# 			my($accnum, $index) = split / /, $acnLine;
# 			$RHaccToIndex->{$accnum} = $index;
# 			#print "Acc# $accnum - Index# $index\n";
# 			$entryCounter++;
# 			unless($verboseProgressOff){
# 				if ( $entryCounter % 10000 == 0 ) {
# 					print "."; #visual progress indicator to screen
# 				}
# 			}
# 		}
# 		close(ACN);
# 	} else {
# 		# open file in binary mode
# 		open(ACN, $acnFile) or die "Can't open $acnFile : $!\n";
# 		binmode(ACN);
# 		until ( eof(ACN) ) {
# 			read(ACN, $record, $recordSize);
# 			my($accnum, $index) = unpack("I I", $record);
# 			$RHaccToIndex->{$accnum} = $index;
# 			#print "Acc# $accnum - Index# $index\n";
# 			$entryCounter++;
# 			unless($verboseProgressOff){
# 				if ( $entryCounter % 10000 == 0 ) {
# 					print "."; #visual progress indicator to screen
# 				}
# 			}
# 		}
# 		close(ACN);
# 	}

# 	#Now read in the idx or .idx2 file, and create a map of index's to fasta file offets
# 	#KRC 7/28/2021 revision to support idx2 file reading with 64-bit PERL
# 	#while the offsets for .idx2 are 64-bit ints, the initial records for num_entries, max_size remain 32-bit
# 	my($idxFile) = $db . ".idx";
# 	if($offsets64Bit){
# 		$idxFile = $db . ".idx2";
# 	}
# 	unless($verboseProgressOff){
# 		print "\n<BR> Processing index file: $idxFile \n";
# 	}
# 	open(IDX, $idxFile) or die "Can't open $idxFile : $!\n";
# 	binmode(IDX);
# 	# while not eof, read in <idx#><offset> pair
# 	# look at optimizing this.. for now just one field at a time
# 	# create a hash mapping idx# to file offset
# 	# (I think the first two values, though are:
# 	#   num_entries, max_size, then array of offsets )
# 	my($firstRecordSize) = 2*$sizeOfInt;
# 	# read in num_entries and max_size
# 	read(IDX, $record, $firstRecordSize);
# 	my($num_entries, $max_size) = unpack("I I", $record);	#.idx & .idx2 the initial records for num_entries, max_size remain 32-bit
# 	#my($num_entries, $max_size) = unpack("Q Q", $record);	#.idx & .idx2 the initial records for num_entries, max_size remain 32-bit

# 	#KRC 7/29/2021 offsets for idx2 are int64
# 	if($offsets64Bit){
# 		$sizeOfInt = $Config{i64size};
# 	}
# 	my($offsetSize) = $sizeOfInt;

# 	# now read in array of offsets
	LindexOffsets = []
# 	my($idx)=0;
# 	until ( eof(IDX) ) {
# 		read(IDX, $record, $offsetSize);
# 		my($offset);
# 		if($offsets64Bit){
# 			$offset = unpack("Q", $record);
# 		}else{
# 			$offset = unpack("I", $record);
# 		}
# 		push @AindexOffsets, $offset;
# 		$idx++;
# 		#print "Index#$idx  - offset: $offset\n";
# 		unless($verboseProgressOff){
# 			if ( $idx % 10000 == 0 ) {
# 				print "."; #visual progress indicator to screen
# 			}
# 		}
# 	}
# 	close(IDX);
# 	
# 	#store in case requested again
# 	$initializedDataBase = $db;
# 	$RAindexOffsets = \@AindexOffsets;
# 	
# 	# now access the database file to get the entry informations for each accession key
# 	unless($verboseProgressOff){ print "\n<BR> Parsing entry information \n"; }
# 	my($numEntriesRequested) = $#$RAaccessionKeys + 1;
	if(len(LaccessionKeys) > 0):
		DDproteinsByAccessionKey = fParseFastaFileUsingIndicesToAccessionNumberOrientedRHRHprotein(sDBpath, LaccessionKeys, LindexOffsets, DaccToIndex)
	else:
# 		#get them all
		DDproteinsByAccessionKey = fParseFastaFileUsingIndicesToAccessionNumberOrientedRHRHprotein(sDBpath, DaccToIndex.keys(), LindexOffsets, DaccToIndex);

	return(DDproteinsByAccessionKey)

#####################################################################
## fParseFastaFileUsingIndicesToAccessionNumberOrientedRHRHprotein ##
#####################################################################
#s
def fParseFastaFileUsingIndicesToAccessionNumberOrientedRHRHprotein(sDB, LaccessionKeys, LindexOffsets, LDaccToIndex):

# 	open(DB, $db) or die "Can't open $db : $!\n";
# 	# text file

# 	my($accKey, $RHfprotein);
	DDproteinsByAccessionKey= {}	# = &makeNewHash(); #dummy 1st protein
# 	my($numidx) = 0;
# 	foreach $accKey (@$RAaccessionKeys) {
# 		$numidx++;
# 		my($idx) = $RHaccToIndex->{$accKey};	
# 		if ($idx) {
# 			my($seekOffset) = $$RAindexOffsets[$idx - 1];	# 0-based access
# 			#print "\n\nReading entry for Acc# $accKey - Idx# $idx at Offset $seekOffset\n";
# 			seek(DB, $seekOffset, 0) or die "Failed to seek to offset $seekOffset : $!\n";
# 			my($line);
# 			# note: offset is always past the leading '>'
# 			$line = <DB>;	# read in header line
# 			$RHfprotein = &makeNewFastaProteinHash();
# 			$RHfprotein->{'header'} = $line;
# 			# now read in sequence
# 			my($sequence) = "";
# 			while ($line = <DB>) {
# 				chomp $line;
# 				if ($line =~ m/\>/) {
# 					last;
# 				} 
# 				$sequence .= $line;
# 			}
# 			# TODO:? replace octal 001 (ncbi indicates next gi) with \n.
# 			$RHfprotein->{'sequence'} = $sequence;
# 			# aligner doesn't use the sequenceLength and accession_number
# 			# but I assigned them to keep this consistent with earlier code - JR
# 			$RHfprotein->{'sequenceLength'} = length($sequence);
# 			$RHfprotein->{'accession_number'} = $accKey;
# 		} else {
# 			# this can happen if first hit's search was say .stdmix and second's was .mammals
# 			# Could always use the 'super set db' (NCBInr), but that takes a long time.
# 			$RHfprotein = &makeNewFastaProteinHash();
# 			$RHfprotein->{'header'} = "Accession# not found in $db";
# 			$RHfprotein->{'sequence'} = "";
# 		}
# 		$RHRHproteinsByAccessionKey->{$accKey} = $RHfprotein;
# 		if ( $numidx % 10 == 0 ) {
# 			#when called to fill protein names for all hits
# 			#this print happens for each row in the protein table but is between cells, so it gets shown before the entire table
# 			#without being on the same line, so don't print it KRC 3/12/2008
# #			print "."; #visual progress indicator to screen
# 		}
# 	}
# 	# Now, for each accessionKey in list, get the entry header and sequence

# 	close (DB);	
# #s
	return(DDproteinsByAccessionKey)
