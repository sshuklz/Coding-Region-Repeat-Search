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

        DomainLen = 0
        ResNum = 0
        GapNum = 0
        
        for CodonPos in range (round(len(row[2])/3)):
            DomainLen += GapLength
            RepLen += GapLength
            if row[2][DomainLen-3:DomainLen] in Reps : # repeat residues
                ResNum += 1
                GapNum = 0

            else:
                GapNum += GapLength

                if GapNum >= Gap_allowance*GapLength:

                    if ResNum >= Size: 
                        
                        Gene_ref = row[0] + ': ' + row[1]; Up_down_stream_adj = 0
                        Coord = [((CodonPos*GapLength) - RepLen + 3),(CodonPos*GapLength)+(GapLength*2) - GapNum - 3]
                        
                        with closing(requests.get(RefTable, stream=True)) as r:
                            reader = csv.reader(codecs.iterdecode(r.iter_lines(), 'utf-8'), delimiter=',') 
                            for ref in reader:
                               
                                if ref[1] == row[1]:
                                    stop = 0; Rep_ST = 0; Rep_ED = 0
                                    
                                    if int(ref[8]) == 1:
                                            
                                        if Coord[0] > int(ref[7]) or Coord[1] < int(ref[6]):
                                            stop = 1; Up_down_stream_adj = 0
                                        
                                        Rep_ST = int(ref[3])
                                        Rep_ED = int(ref[4]) 
                                        
                                        if Coord[0] > int(ref[6]):
                                            Rep_ST = int(ref[3]) + (Coord[0] - int(ref[6])) - 1
                                            
                                        if Coord[1] < int(ref[7]):
                                            Rep_ED = int(ref[4]) - (int(ref[7]) - Coord[1]) - 1 + Up_down_stream_adj
                                               
                                    if int(ref[8]) == -1:
                                            
                                        if Coord[0] > int(ref[7]) or Coord[1] < int(ref[6]):
                                            stop = 1; Up_down_stream_adj = 0
                                            
                                        Rep_ST = int(ref[3])
                                        Rep_ED = int(ref[4])
                                        
                                        if Coord[0] > int(ref[6]):
                                            Rep_ED = int(ref[4]) - (Coord[0] - int(ref[6]))
                                            
                                        if Coord[1] < int(ref[7]):
                                            Rep_ST = int(ref[3]) + (int(ref[7]) - Coord[1]) - Up_down_stream_adj
                                    
                                    Elen = int(ref[4]) - int(ref[3]); Rlen = Rep_ED - Rep_ST;
                                    
                                    if stop != 1:       
                                        Gene_dict = Gene_dict.append({'Gene': ref[0],'Transcript': ref[1],'CHR': ref[2],'Exon_num':ref[5],'Exon_ST': int(ref[3]),'Exon_ED': int(ref[4]),'Rep_ST': Rep_ST,'Rep_ED': Rep_ED,'Coding_LEN': len(row[2]),'Exon_LEN': Elen,'Rep_LEN': Rlen,'Sense': ref[8]}, ignore_index=True)
                                        Up_down_stream_adj = 1
                                    else:
                                        Up_down_stream_adj = 0
                                    
                        DNAseq = row[2][(CodonPos*GapLength)-RepLen:(CodonPos*GapLength)+(GapLength*2) - GapNum]
                        AAseq = Seq(DNAseq, generic_dna).translate()
                        
                        if row[0] not in Gene_list or Gene_ref in Gene_ref_list:
                            if row[0] not in Gene_list:
                                Gene_list = Gene_list + [row[0]]


## Project updates:

Jul 24 
- Implemented dynamic ICGC dynamic data table creation on Collab page. Unfortunately due to memory problems sessions often times run out before tabulations are completed.
Possible options include: upgrading to Colloab pro for $10 a month subscription (https://colab.research.google.com/signup). Trim data oon import, only query neccasary data points
and cut out the fluff. Perform calculations offline and move away from google collab. 
- wrote code in HTML to decipher JSON objects produced by ICGC API
- Rewrote / organized code in accordance to PEP 8 style guide (https://www.python.org/dev/peps/pep-0008/), so code is understandable by a large audience
- Updates from meeting: Try to minimize data size by removing enseble ID. Wait for subscription to get approved and also recieve ideal model table. Work on searching for intronic regions as opposed to exonic regions (will have to create new dataset and reorganize code slightly). READ literature.
