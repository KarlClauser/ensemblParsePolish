##############################################################################
#	Copyright (c) 2022-2023 Broad Institute of MIT and Harvard
#		All rights reserved
##############################################################################

#Filename:	lsmFastaHeaders.py
#Purpose:	Functions for manipulating headers in FASTA formatted sequence databases.
#Authors:	Karl Clauser 
#Created:	March 2022
#Adapted from SM perl file lsmFastaHeaders.pl

import re # split on regex
import pandas as pd

##############################################################################
############# fParseHeadersTo_accession_species_geneSymbol_name ##############
##############################################################################
def fParseHeadersTo_accession_species_geneSymbol_name(LDfproteins, sDatabaseFilename):

	#fill category dataframe
	for Dfprotein in LDfproteins:
		#if 'accession_number' not in Dfprotein:
		Dfprotein['accession_number']	= fAccessionNumberFromHeader(Dfprotein['header'], sDatabaseFilename)
		Dfprotein['entryName'], Dfprotein['species'] = fEntryNameFromHeader(Dfprotein['header'], sDatabaseFilename)
		if Dfprotein['species'].startswith('nuORF') and len (Dfprotein['entryName'].split(' ')) >= 2:
			Dfprotein['transcriptType'] = ' '.join(Dfprotein['entryName'].split(' ')[1:]) #all except first
			Dfprotein['genomeCoordinates'] = Dfprotein['accession_number'].split('_')[-1] #last item
			Dfprotein['chr'] = Dfprotein['genomeCoordinates'].split(':')[0] #first item
		if 'OS=' in Dfprotein['header']:
			Dfprotein['speciesOS']			= fParseSpeciesOSfromHeader(Dfprotein['header'])
		Dfprotein['geneSymbol']				= fParseGeneSymbolFromHeader(Dfprotein['header'])

	return (LDfproteins)
##############################################################################
####### fParseHeadersTo_accession_species_geneSymbol_name_CategoriesFile #####
##############################################################################
def fParseHeadersTo_accession_species_geneSymbol_name_CategoriesFile(LDfproteins, sDatabasePath, sDatabaseFilename, bSkipOutput = False):

	LDfproteins = fParseHeadersTo_accession_species_geneSymbol_name(LDfproteins, sDatabaseFilename)

	print ('<BR>First accession_number: ', LDfproteins[0]['accession_number'], ' Last accession_number: ', LDfproteins[len(LDfproteins) - 1]['accession_number'])

	#convert list of dictionaries to dataframe without the header and sequence
	dfProteins = pd.DataFrame(LDfproteins)
	#dfProteins.drop(columns=['sequence', 'header'], inplace=True)
	dfProteins.drop(columns=['sequence'], inplace=True)
	dfProteins = dfProteins[ ['accession_number'] + [ col for col in dfProteins.columns if col != 'accession_number' ] ] # move accession_number to left of table
	if 'header' in dfProteins.columns:
		dfProteins = dfProteins[[ col for col in dfProteins.columns if col != 'header' ] +  ['header']] # move header to right of table

	if not bSkipOutput:
		#create the category filename
		if 'transcriptType' in dfProteins.columns and 'genomeCoordinates' in dfProteins.columns:
			sCategoryFilePath = sDatabasePath + ".SGSNTGC.categories.tsv";
		else:
			sCategoryFilePath = sDatabasePath + ".SGSN.categories.tsv";

		#write the dataframe to a file
		dfProteins.to_csv(sCategoryFilePath , sep='\t', index=False)
		print ('<BR>Wrote file ', sCategoryFilePath)

	return (dfProteins)
#####################################################################
##################### fAccessionNumberFromHeader ####################
#####################################################################

def fAccessionNumberFromHeader(sHeader, sDatabaseFilename):

	sDatabaseFilename = sDatabaseFilename.lower() 	#enable below startswith statements to become case insensitive
	sHeader = sHeader.lstrip('\>') #remove leading >
	sAccessionNumber = ''

	if sHeader.startswith(('sp|', 'tr|')): #use a tuple - double parens for multiple
		sHeader = sHeader.replace('sp|', '')	#strip sp|
		sHeader = sHeader.replace('tr|', '')	#strip tr|
		sAccessionNumber = fAccessionNumberFromHeaderPreceedingDelimiter(sHeader, '\|')
	elif sHeader.startswith('ENS') or sDatabaseFilename.startswith('ensembl') or sDatabaseFilename.startswith('gencode'):
		sAccessionNumber = fAccessionNumberFromHeaderPreceedingDelimiter(sHeader, " |\|"); #space or |
	else:
		sAccessionNumber = fAccessionNumberFromHeaderPreceedingDelimiter(sHeader, " |\|"); #space or |

	return(sAccessionNumber)

#####################################################################
########### fAccessionNumberFromHeaderPreceedingDelimiter ###########
#####################################################################
def fAccessionNumberFromHeaderPreceedingDelimiter(sHeader, sDelimiter_re):

	#Lentries = sHeader.split(sDelimiter)
	Lentries = re.split(sDelimiter_re, sHeader) #split on regex
	sAccessionNumber = Lentries[0]
	sAccessionNumber = sAccessionNumber.strip() #remove spaces, leading or trailing

	return(sAccessionNumber)

#####################################################################
#################### fParseGeneSymbolFromHeader #####################
#####################################################################
def fParseGeneSymbolFromHeader(sHeader):

	return(fExtractString(sHeader,'GN='))
#####################################################################
#################### fParseSpeciesOSfromHeader ######################
#####################################################################
def fParseSpeciesOSfromHeader(sHeader):
	sSpeciesOS = ''
	sTagNotRegex = 'OS='
	if sTagNotRegex in sHeader:
		sValuePlusTrailers = sHeader.split(sTagNotRegex)[1] #trim off everything preceding tag

		Lvalues = sValuePlusTrailers.split()
		sSpeciesOS = Lvalues[0]				#genus
		if len(Lvalues) >= 2:
			sSpeciesOS += ' ' + Lvalues[1]	#genus + species separated by space

	return(sSpeciesOS)
#####################################################################
########################## fExtractString ###########################
#####################################################################
def fExtractString(sHeader,sTagNotRegex):
	sValue = ''
	if sTagNotRegex in sHeader:
		sValuePlusTrailers = sHeader.split(sTagNotRegex)[1] #trim off everything preceding tag
		sValue = sValuePlusTrailers.split()[0] # trim off whitespace and everything after, empty split args and split will take into account multiple spaces, newlines, tabs and all other forms of whitespace.
		sValue = sValue.split('|')[0] # trim off everything after |.

	return(sValue)
#####################################################################
######################## fEntryNameFromHeader ########################
#####################################################################
def fEntryNameFromHeader(sHeader, sDatabaseFilename):

	sEntryName = -1 #return value for unsupported databases
	sSpecies = -1
	
	sHeader = sHeader.lstrip('\>') #remove leading >
	sAccessionNumber = ''

	if sHeader.startswith(('sp|', 'tr|')):
		sEntryName, sSpecies = fEntryName_species_FromHeaderUniProt(sHeader)
# 	elif($databaseFile =~ /^NCBInr/i
# 		|| ($databaseFile =~ /^ipi/i && $header =~ /^gi/i)  #NCBInr entries added to IPI
# 		){
# 		my($delimiter) = "gi\\|";
# 		my($junk, @entries) = split(/$delimiter/,$header);
# 		my($numEntries) = $#entries + 1;
# 		chop($entries[0]);#removes \n or control character joining entries
# 		($species, $entryName) = &speciesAndNameFromSplitNCBInrHeader($entries[0]);
	elif sHeader.startswith('ENS'):
		##Public 7/2013 Homo_sapiens.GRCh37.72.pep.all.fa.gz from http://useast.ensembl.org/info/data/ftp/index.html
		##>ENSP00000468826 pep:known chromosome:GRCh37:HG27_PATCH:26597180:26600278:1 gene:ENSG00000268612 transcript:ENST00000593546 gene_biotype:protein_coding transcript_biotype:protein_coding
		##9/3/2019 from Sarah Calvo
		##>ENSBTAP00000064819|ENSBTAG00000052228|GN= TRBV24 - 1 T cell receptor beta variable 24 - 1[Source:HGNC Symbol; Acc:HGNC : 12203]
		##>ENSBTAP00000058331|ENSBTAG00000049985|GN=
		#>ENSBTAP00000002299|ENSBTAG00000001753|GN= chloride intracellular channel 6 [Source:NCBI gene;Acc:507243]
		##10/30/2020 from Karl Ensembl.human.hg19.clean3nr
		#>ENSP00000467141.1|ENST00000589042.1|ENSG00000155657.19|chr2:179390716-179672150:+|GN=TTN titin									Ensembl.human.hg19.clean3nr
		#>ENSP00000342313.3-W343R|GN=SAMD11|transcriptID=ENST00000342066.3|chr=1|type=G|qual=0.01|SAAP=W343R|SAMD11 p.W343R c.1027T>C		DFCI-5485-P01_fullpeptide_variant_proteomenr.fasta
		#>ENST00000488147.1_1_1:17367-18318:-|GN=WASH7P Pseudogene|nuORFdb_human															nuORFdb derived
		#>A1BG_19:58858248-58858388:-|GN=A1BG|smORF_v4|dORF|																				smORF_v4
		#print "\n<br>Ensembl or Gencode";

		Lfields = sHeader.split('|'); #part after 1st space
		if sHeader.startswith('ENSP'):
			sSpecies = 'Human'
			sEntryName = Lfields[len(Lfields)-1] #last field
		elif  sHeader.startswith('ENSBTA'):
			sSpecies = 'Bovine'
			sEntryName = Lfields[len(Lfields)-1] #last field
		else: #parse as if it is SM generic format >acc|name|species
			sEntryName, sSpecies = fEntryName_species_FromHeaderSMgeneric(sHeader) #use this which tidies up leading trailing spaces

		#$entryName =~ s/^GN=//g; #strip off leading indicator of Gene symbol
	else: #parse as if it is SM generic format >acc|name|species
		sEntryName, sSpecies = fEntryName_species_FromHeaderSMgeneric(sHeader) #use this which tidies up leading trailing spaces
# 	}elsif ($databaseFile =~ /^RefSeq/i) { #KRC 3/11/2020 revised version
# 		#from David Fenyo
# 		##>NP_001254479 titin isoform IC [Homo sapiens]. (MAP:chr2:179391738- 
# 		#from REfSeq website
# 		##>gi|341915075|ref|XP_003403591.1| PREDICTED: hypothetical protein LOC100652771 [Homo sapiens]
# 		
# 		#KRC 3/11/2020 also handle lines with  multiple items in brackets
# 		#>NP_060024.2 poly [ADP-ribose] polymerase 14 GN=PARP14 [Homo sapiens]
# 		
# 		#my($plusEntryName, $speciesPlus) = split(/\[/,$header);
# 			#$species = $speciesPlus;
# 		#my($junk2);
# 		#($species,$junk2) = split(/\]/,$speciesPlus);
# 		#keep everything after 1st space
# 		my(@AbracketedFields) = split(/\[/,$header);
# 		my($numBracketedIndices) = $#AbracketedFields;
# 		
# 		if($numBracketedIndices == 0
# 			&& $header =~ /\|/
# 			){
# 			($entryName, $species) = &entryName_species_FromHeaderSMgeneric($header);
# 		}else{
# 			#get species
# 			my($junk2);
# 			($species,$junk2) = split(/\]/, $AbracketedFields[$numBracketedIndices]); #use last bracketed field
# 			#keep everything after 1st space
# 			my($junk, @AspaceDelimited) = split(/ /, $AbracketedFields[0]);	#take out spaces
# 			my($proteinName0)= join(" ", @AspaceDelimited);					#put spaces back in
# 			
# 			if($numBracketedIndices > 1){
# 				$AbracketedFields[0] =$proteinName0; #exclude everything before the name
# 				$entryName = join(" [", @AbracketedFields[0 .. $numBracketedIndices-1]); #exclude the last bracketed field
# 			}else{
# 				$entryName = $proteinName0; 	#removes \n if present
# 			}
# 		}
	sEntryName = sEntryName.lstrip('GN=')#strip off leading indicator of Gene symbol)'

	return(sEntryName, sSpecies)
#####################################################################
############### fEntryName_species_FromHeaderSMgeneric ##############
#####################################################################
def fEntryName_species_FromHeaderSMgeneric(sHeader):
#e
	##>ENSMUSP00000027035|Sox17|MOUSE|ENSMUST00000027035|ENSMUSG00000025902|Transcription factor SOX-17. [Source:Uniprot/SWISSPROT;Acc:Q61473]
	sEntryName =sSpecies = '' #initialize with empty values
	Lentries = sHeader.split('|')

	iNumEntries = len(Lentries)
	if iNumEntries > 1:
		sEntryName = Lentries[1].strip() #remove leading & trailing spaces
		if iNumEntries > 2:
			sSpecies = Lentries[2].strip() #remove leading & trailing spaces

	return(sEntryName, sSpecies)

#####################################################################
################ fEntryName_species_FromHeaderUniProt ###############
#####################################################################
def fEntryName_species_FromHeaderUniProt(sHeader):

	sEntryName =sSpecies = '' #initialize with empty values
	# new format - split at spaces following first '_'
	##>sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens GN=YWHAB PE=1 SV=3
	##>sp|A2ASS6|TITIN_MOUSE Titin OS=Mus musculus GN=Ttn PE=1 SV=1
	Lunderscores = sHeader.split('_', 1) #split on first underscore
	if len(Lunderscores) == 2:
		sSpeciesName = Lunderscores[1].split(' OS=')[0]  # remove all text following name
		LspeciesName = sSpeciesName.split(' ', 1) #split on first space
		if len(LspeciesName) == 2:
			return(LspeciesName[1], LspeciesName[0])

	return(sEntryName, sSpecies)
