# ensemblParsePolish
Python script to create a base reference protein sequence database derived from a human Gencode release (ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/) consisting of all non-redundant protein coding transcript types mapped to the human reference genome GRCh38 with length > 6 amino acids.

Function of ensemblParsePolish.py script
The protein sequence file that one can download from Gencode, gencode.vXX.pc_translations.fa.gz, has a few undesirable characteristics for subsequent use in proteomics for database searching of LC-MS/MS.  
These include:
  - All transcripts not just the ones tagged by Gencode as 'basic' or 'GENCODE_Primary'.
  - Length < 7
  - Redundant sequences present
  - Headers lack a protein name (gene description) and genome coordinates
 
The script ensemblParsePolish.py addresses each of the above issues to yield a reference proteome better suited for database searching of LC-MS/MS data.  
For troubleshooting files of removed and retained entries are kept at each stage.

A single running of the script can make 1 of 4 possible subsets (the 3 below and one requiring both basic and primary tags):    
  - GENCODE Comprehensive, which captures all annotated transcripts  
  - GENCODE Basic, a smaller subset of transcripts containing only full-length transcripts at protein-coding genes and enriched for alternative splicing at lncRNA genes. The length and depth of the reads generated in ONT and PacBio experiments will eventually overwhelm the GENCODE Basic category, for example by adding many tens or hundreds of full-length alternatively spliced transcripts at protein-coding genes.  (average 6.1 transcripts/gene)  
  - GENCODE Primary transcripts (v47 +) capture a minimal set of transcripts at protein coding genes and include protein-coding exons with evidence of evolutionary constraint and alternative splicing (including exon skips) with high expression, includes all MANE Select/Ensembl Canonical and MANE Plus Clinical transcripts by default. (average 1.9 transcripts/gene)

The script was written and tested on Windows operating systems. Minor revision may be required for Linux.

Example output of 3 possible subsets for Gencode v50 can be downloaded from:  
sSubsetTagGencode supported values 'GENCODE_Primary', 'basic', 'basicAndPrimary' (requires both tags)
  - https://proteomics.broadinstitute.org/ExampleDB/Gencode50_Ensembl116.zip  

Prerequisites to run:
---------------------
Create a working directory for inputs/outputs  
&emsp;C:\Gencode50_Ensembl116

Download required inputs from ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/
&emsp;(copy/paste ftp URL to address bar of a Windows Explorer window):

&emsp;gencode.v50.pc_translations.fa.gz\
&emsp;gencode.v50.basic.annotation.gtf.gz\
&emsp;gencode.v50.annotation.gtf.gz

Obtain annotations from BioMart website  
&emsp;https://www.ensembl.org/biomart/martview  
&emsp;Set the Dataset to Human genes (GRCh38.p14)

&emsp;Click the  checkboxes for the following attributes in the BioMart UI:  
&emsp;&emsp;Protein stable ID version, Transcript stable ID version, Gene stable ID version,  
&emsp;&emsp;Chromosome/scaffold name, Transcript start (bp), Transcript end (bp),  
&emsp;&emsp;Strand, Transcript type, HGNC symbol, Gene description

&emsp;Export as TSV format and save the file to the working directory with the following name  
&emsp;C:\Gencode50_Ensembl116\Ensembl-v50_GRCh38_mart_export.txt

Run the script from the command line:
-------------------------------------
With directory structure  
C:\ensemblParsePolish\  
C:\Gencode50_Ensembl116\  

Run  
C:\>cd ensemblParsePolish  
C:\ensemblParsePolish>python ensemblParsePolish.py  

Incorporate into a pipeline with another python script that treats ensemblParsePolish as a library:
---------------------------------------------------------------------------------------------------
>import sys\
>sys.path.append(r'C:\\ensemblParsePolish')\
>import ensemblParsePolish as EPP  
>
>#Revise versions to run with a later release  
>#Designate a subset to make  
>#sSubsetTagGencode supported values: '' (all), 'GENCODE_Primary', 'basic', 'basicAndPrimary' (requires both tags)   
>#GENCODE Comprehensive, which captures all annotated transcripts  
>#GENCODE Basic, a smaller subset of transcripts containing only full-length transcripts at protein-coding genes and enriched for alternative splicing at lncRNA genes. The length and depth of the reads generated in ONT and PacBio experiments will eventually overwhelm the GENCODE Basic category, for example by adding many tens or hundreds of full-length alternatively spliced transcripts at protein-coding genes.  (average 6.1 transcripts/gene)  
>#GENCODE Primary transcripts (v47 +) capture a minimal set of transcripts at protein coding genes and include protein-coding exons with evidence of evolutionary constraint and alternative splicing (including exon skips) with high expression, includes all MANE Select/Ensembl Canonical and MANE Plus Clinical transcripts by default. (average 1.9 transcripts/gene)  
>
>EPP.fMain(sSubsetTagGencode = 'GENCODE_Primary', sRefProteomeDir = r'C:\Gencode50_Ensembl116', sGencodeVersion = 'v50', sEnsemblVersion = 'Ensembl-116')\	  
