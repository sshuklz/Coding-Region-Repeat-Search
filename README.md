# Coding-Region-Repeat-Search
Project dedicated towards generating data for proteins with repeating elements (repeats of interest specified by user input). Data pertains currently to cancer mutation data from cBioPortal, but the project will be expanded to pull data from other databases such as ProteinAtlas, StringDB, etc.

## Code workflow

1. Static data set import function from github page
2. User defined inputs for specified search parameters as follows:
* Reps: Give a list of amino acids in quotes i.e: 'CAG', 'GAA', 'GCGCGTT' inside square brackets [], note if multiple repeats entered at once i.e: ['CAG', 'GTA'] they must have same length
* Size: Give repeat stretch minimum size i.e: for 'CAG' and Size = 4 repeats the following are returned: CAGCAGCAGCAG, CAGCAGCAGCAGCAGCAG but not CAGCAG
* Gap_allowance: Strech of nucleic bases allowed inbetween repeats ie. if Gap_allowance = 1 and 'CAG' repeats, CAG - 'TTT' - CAGCAGCAGCAG - 'CATGGG' - CAGCAGCAGCAGCAG output is CAGTTTCAGCAGCAGCAG
3. Search for genes with repetitive domains and accompanying genomic coordinates presented in tabular form
4. Fetch mutation data from online databases such as CbioPortal, CCLE, and ICGC
5. Perform calculations for mutation frequency and abundance in repitive regions compared to other regions
6. Present results in tabular / graph form

## Working static dataset explainers on Github page

1. CDSlibrary = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/Human_full_CDS_libraryV2.csv'

CSV file that contains the following three columns: Gene name, Ensemble Transcript ID, DNA sequence of transcript obtained from grch37 Ensemble biomart (https://grch37.ensembl.org/biomart/martview/ac6be9fa83be44339e854dfa8eef10e7) 

2. RefTable = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/Exon2SeqDictionary.csv'

CSV file that contains the following columns: Gene name, Transcript stable ID version, Chromosome/scaffold name, Exon region start (bp), Exon region end (bp), Exon rank in transcript, CDS start, CDS end, 5' UTR end, 3' UTR start
obtained from grch37 Ensemble biomart (https://grch37.ensembl.org/biomart/martview/ac6be9fa83be44339e854dfa8eef10e7) 

3. CancerStudy = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/cancer_study_id.csv'

CSV file that contains a list of PanAtlas Cancer studies used in CbioPortal search

4. GeneName2EnsembleID = 'https://raw.githubusercontent.com/sshuklz/Coding-Region-Repeat-Search/master/Datasets/GeneName2EnsembleID.csv'

CSV file that contains two columns: gene name and accompanying Ensemble ID (multiple IDs for same gene due to halpotypes)

## Calculations and positioning consideration (in progress)

## Project updates:

Jul 24 
- Implemented dynamic ICGC dynamic data table creation on Collab page. Unfortunately due to memory problems sessions often times run out before tabulations are completed.
Possible options include: upgrading to Colloab pro for $10 a month subscription (https://colab.research.google.com/signup). Trim data oon import, only query neccasary data points
and cut out the fluff. Perform calculations offline and move away from google collab. 
- wrote code in HTML to decipher JSON objects produced by ICGC API
- Rewrote / organized code in accordance to PEP 8 style guide (https://www.python.org/dev/peps/pep-0008/), so code is understandable by a large audience
- Updates from meeting (in progress)
