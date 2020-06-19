#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:06:49 2019

@author: Shalabh
"""

import numpy as np
import csv
    
m = 1    
## DeNovo refrence table construction form FASTA file download

Allele = np.array([]) # Allele and AAseq are array elements which contain those respective information
AAseq = np.array([])
CombStr = '' # fasta files split AAseq data by 60AA a line. Use .join to combine split strings
n = 0 # to initiate first loop iteration, can use row number instead but this works
 
with open ('CDNA_Human.csv', 'r') as csv_file:             # Loop taking forever, current time taken is 267 seconds for full (V2 w/ ISOonly is 118 seconds)
    csvread = csv.reader(csv_file, delimiter='|')
    
    for row in csvread:
        #row[0][0] == '>' denotes row as header row containing header information such as allele name in file
        
        if n >= 1 and row[0][0] == '>':
            
                                      #CombStr needs at least one AA before it is entered in. On first run this condition is not satisfied
                AAseq = np.append(AAseq,CombStr)    #store CombStr in AAseq corrosponding to header allele
                
                Allele = np.append(Allele,(row[0])
                m+=1
                
            CombStr = ''                        #clear CombStr for next header allele
            
        elif n == 1:
            CombStr = "".join([CombStr, "".join(row)])  #Iterative joining of AA strings as long as it is not a header row in fasta file
       
        elif n == 0:                                    #first cylce allele name is entered, with no prior CombStr entry
            Allele = np.append(Allele,(row[10].replace(' description','')))    #store Header allele name in Allele adjust as needed
            n += 1
            
AAseq = np.append(AAseq,CombStr) # last CombStr entry as loop has finished
CombStr = ''

RefTable = np.c_[Allele,AAseq] #everything can be stored directly into RefTable in loop, so no need for Allele,AAseq in final implementaion

CSVname = input("Save as .csv filename (i.e hla_prot_DB2.csv): ")
np.savetxt(CSVname, RefTable, fmt='%s')