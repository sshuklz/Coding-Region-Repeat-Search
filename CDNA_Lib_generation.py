import csv
import pandas as pd

DF1 = pd.DataFrame(columns=['Gene','Transcript','DNA'])
n = 0
m = 1
CombStr = '' 
Gene = ''
Transcript = ''

with open ('CDNA_Human.csv', 'r') as csv_file:
    csvread = csv.reader(csv_file, delimiter='|')
    for row in csvread:
        if len(row) != 0 and n >= 1 and row[0][0] == '>':
            
            DNAseq = CombStr 
            DF1 = DF1.append({'Gene': Gene, 'Transcript' : Transcript,'DNA': DNAseq}, ignore_index=True)
            Gene = row[0][1:]
            Transcript = row[1]
            m+=1
            print(m)
            CombStr = ''
            
        elif n == 1 and len(row) != 0:
            CombStr = "".join([CombStr, "".join(row)])  #Iterative joining of AA strings as long as it is not a header row in fasta file
       
        elif n == 0:                                    #first cylce allele name is entered, with no prior CombStr entry
            Gene = row[0][1:]
            Transcript = row[1]   #store Header allele name in Allele adjust as needed
            n += 1
        
DF1.to_csv('Human_full_CDS_library.csv', encoding='utf-8')