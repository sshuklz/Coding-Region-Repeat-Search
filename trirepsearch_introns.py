from Bio.Seq import Seq
from Bio import BiopythonWarning
import itertools
import json
import numpy as np
import os
import pandas as pd
import pickle
import urllib.request
import time
from tqdm import tqdm
import math
from PIL import Image
import warnings

CDSlibrary = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/Human_full_CDS_libraryV2.csv'
RefTable = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/Exon2SeqDictionary.csv'
CancerStudy = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/cancer_study_id.csv'
GeneName2EnsembleID = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/GeneName2EnsembleID.csv'
GeneSyn2EnsembleID = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/GeneName2EnsembleID_V2.csv'
TissueProteinAtlas = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/TissueProteinAtlas.csv'
CellProteinAtlas = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/CellProteinAtlas.csv'
OrganelleProteinAtlas = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/OrganelleProteinAtlas.csv'
CancerProteinAtlas = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/CancerProteinAtlas.csv'
StringDB_GO = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/GoTerms.csv'
Mean_GO_String = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/GO_stringDB.csv'
UTRlib = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/UTR.csv'
ImCell = 'https://raw.githubusercontent.com/sshuklz/Coding-Region-Repeat-Search/master/Custom%20cell/'
Intron_library = 'https://media.githubusercontent.com/media/sshuklz/Coding-Region-Repeat-Search/master/Datasets/Intron_library.csv'

cell_organelles_png = ['a0_full.png',
                       'a1_cytosol.png',
                       'a2_intermediate_filaments.png',
                       'a3_actin_filaments.png',
                       'a4_focal_adhesion_sites.png',
                       'a5_centriolar_satellite.png',
                       'a6_centrosome.png',
                       'a7_microtubules.png',
                       'a8_microtubule_ends.png',
                       'a9_secreted_proteins.png',
                       'b0_lipid_droplets.png',
                       'b1_lysosomes.png',
                       'b2_peroxisomes.png',
                       'b3_endosomes.png',
                       'b4_endoplasmic_reticulum.png',
                       'b5_golgi_apparatus.png',
                       'b6_nucleoplasm.png',
                       'b7_nuclear_membrane.png',
                       'b8_nuclear_bodies.png',
                       'b9_nuclear_speckles.png',
                       'c0_nucleoli.png',
                       'c1_nucleoli_rim.png',
                       'c2_nucleoli_fibrillar_center.png',
                       'c3_rods_and_rings.png',
                       'c4_mitochondria.png',
                       'c5_cell_membrane.png']

folder_names_pickle_jar = ['gene_dict', 'ICGC_data_mut', 'ICGC_data_rep',
                           'Protein_Atlas', 'StringDB','mut_type',
                           'heat_map', 'blosum_score']

datasets = [CDSlibrary, RefTable, CancerStudy, GeneName2EnsembleID,
            GeneSyn2EnsembleID, TissueProteinAtlas, OrganelleProteinAtlas,
            CellProteinAtlas, CancerProteinAtlas, StringDB_GO,
            Mean_GO_String, UTRlib, Intron_library]

dataset_names = ['CDSlibrary', 'RefTable', 'CancerStudy',
                 'GeneName2EnsembleID', 'GeneSyn2EnsembleID',
                 'TissueProteinAtlas', 'OrganelleProteinAtlas',
                 'CellProteinAtlas', 'CancerProteinAtlas', 'StringDB_GO',
                 'Mean_GO_String', 'UTRlib','Intron_library'] ; i = 0

im_path = '../data/Cell_images/'
out_path = '../data/Output_CSVs/'
pickle_path = '../data/Pickle_jar/'
data_path = '../data/Datasets/'

try:

    if not os.path.isfile(out_path):

        os.makedirs(out_path)
        print('Creating output_CSV folder')

except:

    print('Output_CSV folder present')
    pass

for folder_name in folder_names_pickle_jar:

    try:

        if not os.path.isfile(pickle_path + folder_name + '/'):

            os.makedirs(pickle_path + folder_name + '/')
            print('Creating ' + folder_name + ' folder')

    except:

        print(folder_name + ' folder present')
        pass

try:

    if not os.path.isfile(im_path):

        os.makedirs(im_path + 'Custom_cell')
        print('\nCreating image folder')
        print('Downloading cell component images\n')

        urllib.request.urlretrieve(ImCell + 'Cell_normal.png',
            im_path + "Cell_normal.png")

        for cell_component in cell_organelles_png:

            urllib.request.urlretrieve(ImCell + cell_component,
                im_path + "Custom_cell/" + cell_component)

except:

    print('Image folder present')
    pass

try:

    if not os.path.isfile(data_path):

        os.makedirs(data_path)
        print('\nCreating datasets folder')

        for dataset in datasets:

            print('Downloading CSV for ' + dataset_names[i] + ' dataset')
            dataset_csv = pd.read_csv(dataset)

            pd.to_pickle(dataset_csv,os.path.join(data_path,
                dataset_names[i] + '.pickle'))

            dataset_csv.to_csv(data_path + dataset_names[i] + '.csv',
                index = False) ; i += 1

except:

    print('Datasets folder present')

warnings.filterwarnings('ignore', category = RuntimeWarning) 
warnings.filterwarnings('ignore', category = BiopythonWarning) 



def rep_search(bioc_mol_type,
               reps,
               min_size,
               max_size,
               gap_allowance,
               gene_flank_size,
               region):



    Codons = {
               "C": ["TGT", "TGC"],
               "D": ["GAT", "GAC"],
               "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
               "Q": ["CAA", "CAG"],
               "M": ["ATG"],
               "N": ["AAC", "AAT"],
               "P": ["CCT", "CCG", "CCA", "CCC"],
               "K": ["AAG", "AAA"],
               "*": ["TAG", "TGA", "TAA"],
               "T": ["ACC", "ACA", "ACG", "ACT"],
               "F": ["TTT", "TTC"],
               "A": ["GCA", "GCC", "GCG", "GCT"],
               "G": ["GGT", "GGG", "GGA", "GGC"],
               "I": ["ATC", "ATA", "ATT"],
               "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
               "H": ["CAT", "CAC"],
               "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
               "W": ["TGG"],
               "V": ["GTA", "GTC", "GTG", "GTT"],
               "E": ["GAG", "GAA"],
               "Y": ["TAT", "TAC"]
              }

    t_full = time.time()

    search_param = ('_Repeats-'
    + str(reps)
    + '_Min_Size-'
    + str(min_size)
    + '_Max_Size-'
    + str(max_size)
    + '_Gap-'
    + str(gap_allowance)
    + '_Flank-'
    + str(gene_flank_size)
    + '_Region-'
    + (region))

    search_param_short = ('Repeats-'
    + str(reps)
    + '_Min_Size-'
    + str(min_size)
    + '_Max_Size-'
    + str(max_size)
    + '_Region-'
    + (region))



    def AA2NUC(repeat):

        DNA_ambiguous = []

        for AA in range(0, len(repeat)):

            DNA_ambiguous.append(Codons[repeat[AA]])

        DNA_ambiguous = list(itertools.product(*DNA_ambiguous))

        for i in range(len(DNA_ambiguous)):

            DNA_ambiguous[i] = ''.join(DNA_ambiguous[i])

        return(DNA_ambiguous)



    def UTR_length(tscript_id):

        UTR_5_len = 0; UTR_3_len = 0

        for UTR_ref in UTR_library:

            if UTR_ref['Transcript stable ID version'] == tscript_id:

                if int(UTR_ref["5' UTR start"]) != 0:

                    UTR_5_len += (int(UTR_ref["5' UTR end"])
                                - int(UTR_ref["5' UTR start"]) + 1)

                if int(UTR_ref["3' UTR start"]) != 0:

                    UTR_3_len += (int(UTR_ref["3' UTR end"])
                                - int(UTR_ref["3' UTR start"]) + 1)

        return [UTR_5_len, UTR_3_len]



    def repeat_length(gene_name):

        full_repeat_length = 0 ; Repeat_Start = [0]

        for row in gene_dict.to_dict('records'):

            if row['Gene'] == gene_name:

                if int(row['Repeat_Start']) not in Repeat_Start:

                    Repeat_Start += [int(row['Repeat_Start'])]
                    full_repeat_length += int(row['Repeat_Length'])

        return full_repeat_length



    def full_coords(tscript_id, sense):

        coding_coords = [] ; intron_coords = [] ; n_count = 0
        introns_w_rep_coords = [] ; rep_coords = [] ; intron_number = 0
        UTR_5_coords = [] ; UTR_3_coords = [] ; ALL_coords = []
        Flank_U_coords = [] ; Flank_D_coords = []
        Gene_start_UPs = 0 ; Gene_end_DOWNs = 0
        Gene_Flank_UPs = 0 ; Gene_Flank_DOWNs = 0

        gene_index = gene_dict.index[
            gene_dict['Transcript'] == tscript_id].tolist()

        n_count = (gene_dict.at[gene_index[0], 'Coding_DNAseq']).count('N')

        for row in gene_dict.to_dict('records'):

            if row['Transcript'] == tscript_id:

                if int(row['Intron']) != intron_number:

                    intron_number = int(row['Intron'])

                    introns_w_rep_coords += list(range(row['Intron_Start'],
                                                     row['Intron_End'] + 1))

                    rep_coords += list(range(row['Repeat_Start'],
                                             row['Repeat_End'] + 1))

        for exon_ref in Exon_library:

            if exon_ref['Transcript_ID'] == tscript_id:

                coding_coords += list(range(int(exon_ref['Exon_Start']),
                                            int(exon_ref['Exon_End']) + 1))

        for UTR_ref in UTR_library:

            if UTR_ref['Transcript stable ID version'] == tscript_id:

                if int(UTR_ref["5' UTR start"]) != 0:

                    UTR_5_coords += list(range(int(UTR_ref["5' UTR start"]),
                                               int(UTR_ref["5' UTR end"]) + 1))

                if int(UTR_ref["3' UTR start"]) != 0:

                    UTR_3_coords += list(range(int(UTR_ref["3' UTR start"]),
                                               int(UTR_ref["3' UTR end"]) + 1))

        if sense == 1:

            if n_count != 0:

                for n in range(n_count):

                    coding_coords += [min(coding_coords) - 1]

            coding_coords.sort()

            if UTR_5_coords == []:

                UTR_5_coords = [min(coding_coords)]

            if UTR_3_coords == []:

                UTR_3_coords = [max(coding_coords)]

            Gene_start_UPs = min(UTR_5_coords)
            Gene_end_DOWNs = max(UTR_3_coords)
            Gene_Flank_UPs = Gene_start_UPs - gene_flank_size
            Gene_Flank_DOWNs = Gene_end_DOWNs + gene_flank_size
            ALL_coords = list(range(Gene_start_UPs, Gene_end_DOWNs + 1))
            Flank_U_coords = list(range(Gene_Flank_UPs, Gene_start_UPs + 1))
            Flank_D_coords = list(range(Gene_end_DOWNs, Gene_Flank_DOWNs + 1))

        else:

            if n_count != 0:

                for n in range(n_count):

                    coding_coords += [max(coding_coords) + 1]

            coding_coords.sort(reverse = True)

            if UTR_5_coords == []:

                UTR_5_coords = [max(coding_coords)]

            if UTR_3_coords == []:

                UTR_3_coords = [min(coding_coords)]

            Gene_end_DOWNs = min(UTR_3_coords)
            Gene_start_UPs = max(UTR_5_coords)
            Gene_Flank_UPs = Gene_start_UPs + gene_flank_size
            Gene_Flank_DOWNs = Gene_end_DOWNs - gene_flank_size
            ALL_coords = list(range(Gene_end_DOWNs, Gene_start_UPs + 1))
            Flank_U_coords = list(range(Gene_start_UPs, Gene_Flank_UPs + 1))
            Flank_D_coords = list(range(Gene_Flank_DOWNs, Gene_end_DOWNs + 1))

        Full_Gene_Length = abs(Gene_start_UPs - Gene_end_DOWNs) + 1

        for i in gene_index:

            Full_Intron_Length  = (Full_Gene_Length
                                 - gene_dict.at[i,'Coding_Length']
                                 - gene_dict.at[i,'Full_UTR_Length'])

            Non_Coding_Length = (gene_dict.at[i,'Full_UTR_Length']
                               + Full_Intron_Length
                               + (gene_flank_size * 2))

            gene_dict.at[i,'Flank_Up_Stream'] = Gene_Flank_UPs
            gene_dict.at[i,'Flank_Down_Stream'] = Gene_Flank_DOWNs
            gene_dict.at[i,'Gene_Start_Up_Stream'] = Gene_start_UPs
            gene_dict.at[i,'Gene_End_Down_Stream'] = Gene_end_DOWNs
            gene_dict.at[i,'Full_Gene_Length'] = Full_Gene_Length
            gene_dict.at[i,'Full_Intron_Length'] = Full_Intron_Length
            gene_dict.at[i,'Non_Coding_Length'] = Non_Coding_Length

        for coord in ALL_coords:

            if coord not in coding_coords:

                if coord not in UTR_5_coords:

                    if coord not in UTR_3_coords:

                        intron_coords += [coord]

        return [intron_coords, introns_w_rep_coords, rep_coords,
                coding_coords, UTR_5_coords, UTR_3_coords,
                Flank_U_coords, Flank_D_coords]



    def intron_stitch(intron_ref, gene_coord, gene_dict, utr_len):

        rep_st = 0 ; rep_ed = 0
        int_st  = int(intron_ref['Intron_Start'])
        int_ed  = int(intron_ref['Intron_End'])

        if int(intron_ref['Sense']) == 1:

            rep_st = (int_st
                   + (gene_coord[0]
                   - 1))

            rep_ed = (int_st
                    + (gene_coord[1]
                    - 2))

        if int(intron_ref['Sense']) == -1:

            rep_ed = (int_st
                   - (gene_coord[0]
                   + 1))

            rep_st = (int_st
                     - gene_coord[1])
            
            int_st  = int(intron_ref['Intron_End'])
            int_ed  = int(intron_ref['Intron_Start'])

        if rep_st > rep_ed:

            rep_st_new = rep_ed ; rep_ed = rep_st
            rep_st = rep_st_new

        DNAseq = intron_ref['Sequence']
        DNAseq_rep = intron_ref['Sequence'][gene_coord[0] : gene_coord[1]]
        Full_UTR_Length = utr_len[0] + utr_len[1]
        
        Coding_DNAseq = CDS_library[CDS_library['Transcript_ID'].isin([
            intron_ref['Transcript_ID']])].iloc[0]['Sequence']
        
        Coding_Length = len(Coding_DNAseq)

        gene_dict = gene_dict.append({
            'Gene': intron_ref['Gene'],
            'Transcript': intron_ref['Transcript_ID'],
            'CHR': intron_ref['CHR'],
            'Sense': int(intron_ref['Sense']),
            'Flank_Up_Stream' : 0,
            'Flank_Down_Stream' : 0,
            'Gene_Start_Up_Stream' : 0,
            'Gene_End_Down_Stream' : 0,
            'Intron': int(intron_ref['Intron_Num']),
            'Intron_Start': int_st,
            'Intron_End': int_ed,
            'Repeat_Length' : len(list(range(rep_st,rep_ed))) + 1,
            'Repeat_Start': rep_st,
            'Repeat_End': rep_ed,
            'Repeat_DNAseq': DNAseq_rep,
            'Repeat_DNAseq_Count': len(DNAseq_rep),
            'Intron_W_Rep_DNAseq' : DNAseq,
            'Coding_DNAseq' : Coding_DNAseq,
            'Full_Gene_Length' : 0,
            'Coding_Length': Coding_Length,
            'Full_Repeat_Length': 0,
            'UTR5_Length' : utr_len[0],
            'UTR3_Length' : utr_len[1],
            'Intron_W_Rep_Length': len(DNAseq),
            'Full_Intron_Length' : 0,
            'Full_UTR_Length' : Full_UTR_Length,
            'Transcript_Length' : Full_UTR_Length + Coding_Length,
            'Flank_Length' : gene_flank_size,
            'Non_Coding_Length' : 0
            }, ignore_index = True)

        return gene_dict



    def domain_match(gene_name, coord_start, coord_end):

        CodingMut = 'No' ; IntronWRepMut = 'No' ; RepMut = 'No'
        IntronMut = 'No' ; UTR5Mut = 'No' ; UTR3Mut = 'No'
        UpFlankMut = 'No' ; DownFlankMut = 'No' ; IntergeneMut ='No'
        NCodingMut = 'No' ; TranscriptMut = 'No' ; UTRMut = 'No'

        if any(x in coord_dict.at[gene_name,'coding_coords'] \
               for x in [coord_start, coord_end]):

            CodingMut = 'Yes' ; TranscriptMut = 'Yes'

        else:

            if any(x in coord_dict.at[gene_name,'UTR_3_coords'] \
                   for x in [coord_start, coord_end]):

                UTR3Mut = 'Yes' ; UTRMut = 'Yes'
                TranscriptMut = 'Yes' ; NCodingMut = 'Yes'

            else:

                if any(x in coord_dict.at[gene_name,'UTR_5_coords'] \
                       for x in [coord_start, coord_end]):

                    UTR5Mut = 'Yes' ; UTRMut = 'Yes'
                    TranscriptMut = 'Yes' ; NCodingMut = 'Yes'

                else:

                    if any(x in coord_dict.at[gene_name,'intron_coords'] \
                           for x in [coord_start, coord_end]):

                        IntronMut = 'Yes' ; NCodingMut = 'Yes'
                        
                        if any(x in coord_dict.at[gene_name,'introns_w_rep_coords'] \
                               for x in [coord_start, coord_end]):

                            IntronWRepMut = 'Yes'

                            if any(x in coord_dict.at[gene_name,'rep_coords'] \
                                  for x in [coord_start, coord_end]):

                                RepMut = 'Yes'

                    else:

                        if any(x in coord_dict.at[gene_name,'Flank_D_coords'] \
                               for x in [coord_start, coord_end]):

                            DownFlankMut = 'Yes' ; IntergeneMut = 'Yes'
                            NCodingMut = 'Yes'

                        else:

                            if any(x in coord_dict.at[gene_name,'Flank_U_coords'] \
                                   for x in [coord_start, coord_end]):

                                UpFlankMut = 'Yes' ; IntergeneMut = 'Yes'
                                NCodingMut = 'Yes'

        return [CodingMut, UTR5Mut, UTR3Mut, IntronMut, IntronWRepMut, RepMut,
                UTRMut, TranscriptMut, UpFlankMut, DownFlankMut,
                IntergeneMut, NCodingMut]

    def mutation_effect(gene_name, allele_w, allele_m, coord):

        AA_M ='' ; AA_W =''
        CDseq_dict = gene_dict.loc[(gene_dict['Gene'] == gene_name)].head(1)

        if CDseq_dict['Sense'].item() == -1:

            allele_w = str(Seq(allele_w).reverse_complement())
            allele_m = str(Seq(allele_m).reverse_complement())

        CDseq_W = CDseq_dict['Coding_DNAseq'].item() ; last_pos = len(CDseq_W)
        CDseq_coords = coord_dict.at[gene_name, 'coding_coords']
        Mut_pos = CDseq_coords.index(coord)
        CDseq_M = CDseq_W[: Mut_pos] + allele_m + CDseq_W[Mut_pos + 1 :]

        if CDseq_W[Mut_pos] != allele_w:

            print(coord, allele_w, CDseq_W[Mut_pos - 2 : Mut_pos + 3])
            raise Exception('mutation location incorrect')

        if Mut_pos % 3 == 0:

            Codon_M = CDseq_M[Mut_pos : Mut_pos + 3]
            Codon_W = CDseq_W[Mut_pos : Mut_pos + 3]

        elif Mut_pos % 3 == 1:

            Codon_M = CDseq_M[Mut_pos - 1 : Mut_pos + 2]
            Codon_W = CDseq_W[Mut_pos - 1 : Mut_pos + 2]

        elif Mut_pos % 3 == 2:

            Codon_M = CDseq_M[Mut_pos - 2 : Mut_pos + 1]
            Codon_W = CDseq_W[Mut_pos - 2 : Mut_pos + 1]

        for key, value in Codons.items():

            if Codon_M in value:

                AA_M = key

            if Codon_W in value:

                AA_W = key

        if AA_M == AA_W:

            effect = ' (silent)'

        else:

            effect = ' (missense)'

            if AA_M == '*':

                effect = ' (stop gained)'

            elif Mut_pos in [last_pos - 2, last_pos - 1, last_pos]:

                effect = ' (stop lost)'

            elif Mut_pos in [0, 1, 2]:

                effect = ' (start lost)'

        return effect

    """
    Genes with repeat search
    """

    t = time.time()

    if bioc_mol_type == 'Amino acid':

        reps_nuc = []

        for reps_AA in reps:

            reps_nuc += AA2NUC(reps_AA)

        reps = reps_nuc ; rep_nucs = len(reps[0])

    if os.path.isfile(os.path.join(pickle_path + 'gene_dict/',
        'gene_dict' + search_param + '.pickle')):

        full_dict = pickle.load(open(os.path.join(pickle_path + 'gene_dict/',
            'gene_dict' + search_param + '.pickle'), "rb"))

        gene_dict = full_dict[0] ; coord_dict = full_dict[1]
        gene_list = gene_dict.Gene.unique()

        print('\nUnpickled genes with repeats that match search for: '
             + search_param_short
             + '\n \n'
             + str(gene_list))

        print('\nGenerating gene dictionary and coordinate dictionary\n')

    else:

        rep_len = 0 # length of repeat in search
        res_num = 0 # amino acid residue count in search
        domain_len = 0 # full domain sequence count in search
        gap_num = 0 # non repeat in continous seq in search
        gene_coord = [] # DNA sequence gene_coordinates
        rep_nucs = len(reps[0]) # nucleotides in repeat input
        gene_ref_full_list = [] # list of all genes and transcripts
        gene_ref_accepted_list = [] # unique gene and transcript combinations
        gene_list = [] # unique genes after filter
        tscript_list = [] # unique transcripts after filter
        sense_list = [] # corrosponding sense of transcripts

        gene_dict = pd.DataFrame(columns = [
            'Gene',
            'Transcript',
            'CHR',
            'Sense',
            'Flank_Up_Stream',
            'Flank_Down_Stream',
            'Gene_Start_Up_Stream',
            'Gene_End_Down_Stream',
            'Intron',
            'Intron_Start',
            'Intron_End',
            'Repeat_Length',
            'Repeat_Start',
            'Repeat_End',
            'Repeat_DNAseq',
            'Repeat_DNAseq_Count',
            'Intron_DNAseq',
            'Coding_DNAseq',
            'Full_Gene_Length',
            'Coding_Length',
            'Full_Repeat_Length',
            'UTR5_Length',
            'UTR3_Length',
            'Intron_W_Rep_Length',
            'Full_Intron_Length',
            'Full_UTR_Length',
            'Transcript_Length',
            'Flank_Length',
            'Non_Coding_Length'
            ])

        coord_dict = pd.DataFrame(columns = [
            'intron_coords',
            'introns_w_rep_coords',
            'rep_coords',
            'coding_coords',
            'UTR_5_coords',
            'UTR_3_coords',
            'Flank_U_coords',
            'Flank_D_coords'
            ])

        Exon_library = (pickle.load(open(os.path.join(data_path,
            'RefTable.pickle'), "rb")).drop_duplicates()).to_dict('records')

        Intron_library = (pickle.load(open(os.path.join(data_path,
            'Intron_library.pickle'), "rb")).drop_duplicates())

        UTR_library = pickle.load(open(os.path.join(data_path,
            'UTRlib.pickle'), "rb")).to_dict('records')
        
        CDS_library = pickle.load(open(os.path.join(data_path,
            'CDSlibrary.pickle'), "rb"))

        print('\nSearch for genes containing '
              + search_param_short
              + ':\n', flush = True)

        pbar = tqdm(total = Intron_library.shape[0])

        for gene in Intron_library.to_dict('records'):

            pbar.update(1)
            domain_len = 0 ; res_num = 0 ; gap_num = 0 ; frame = -1
    
            try:
    
                for seq in [gene['Sequence'], \
                            gene['Sequence'][1:], \
                            gene['Sequence'][2:]]:
        
                    domain_len = 0 ; res_num = 0 ; gap_num = 0 ; frame += 1
        
                    for codon_pos in range(round(len(seq) / 3)):
        
                        domain_len += rep_nucs ; rep_len += rep_nucs
        
                        if seq[(domain_len - 3):domain_len] in reps:
        
                            res_num += 1 ; gap_num = 0
        
                        else:
        
                            gap_num += 1
        
                            if gap_num > gap_allowance:
        
                                if res_num >= min_size and res_num <= max_size:
                                    
                                    gene_coord = [
                                        ((codon_pos * rep_nucs)
                                            - rep_len
                                            + rep_nucs) + frame,
                                        ((codon_pos * rep_nucs)
                                            + ((gap_num - 1) * rep_nucs)) + frame
                                        ]
        
                                    gene_dict = intron_stitch(gene,
                                                gene_coord,
                                                gene_dict,
                                                UTR_length(gene['Transcript_ID']))
        
                                    gene_ref = (gene['Gene']
                                               + ' ('
                                               + gene['Transcript_ID']
                                               + ')')
        
                                    gene_ref_full_list += [gene_ref]
        
                                    if gene['Gene'] not in gene_ref_accepted_list or \
                                       gene_ref in gene_ref_full_list:
        
                                        if gene['Gene'] not in gene_ref_accepted_list:
        
                                            gene_ref_accepted_list += [gene['Gene']]
        
                                res_num = 0 ; rep_len = 0
            except:

                pass

        pbar.close()

        print('\ntotal search time: '
             + str(float((time.time() - t) / 60))[:5]
             + 'min\n')

        t = time.time()

        print('Sorting and pruning '
              + str(gene_dict.shape[0])
              + ' transcripts: ', flush = True)

        gene_dict = gene_dict.sort_values(by = 'Transcript_Length',
                                          ascending = False)

        gene_dict = gene_dict.drop_duplicates('Repeat_Start')
        gene_dict = gene_dict.drop_duplicates('Repeat_End')
        gene_dict = gene_dict.reset_index(drop = True)
        pbar = tqdm(total = gene_dict.shape[0]) ; i = 0

        for gene in gene_dict.to_dict('records'):

            gene_dict.at[i, 'Full_Repeat_Length'] = repeat_length(gene['Gene'])

            if gene_dict.at[i, 'Gene'] in gene_list:

                if gene_dict.at[i, 'Transcript'] not in tscript_list:

                    gene_dict = gene_dict.drop(i)

            else:

                gene_list += [gene_dict.at[i, 'Gene']]
                tscript_list += [gene_dict.at[i, 'Transcript']]
                sense_list += [gene_dict.at[i, 'Sense']]

            pbar.update(1) ; i += 1

        pbar.close()

        print('\nsort time: '
             + str(float((time.time() - t) / 60))[:5]
             + 'min\n', flush = True)

        t = time.time()

        print('Building gene_dict and coord_dict for: '
              + str(len(gene_list))
              + ' genes', flush = True)

        pbar = tqdm(total = len(gene_list)) ; i = 0

        for gene in gene_list:

            coord_dict.loc[gene] = full_coords(tscript_list[i], sense_list[i])
            pbar.update(1) ; i += 1

        pbar.close()

        print('build time: '
             + str(float((time.time() - t) / 60))[:5]
             + 'min\n', flush = True)
        
        pd.to_pickle([gene_dict,coord_dict],
            os.path.join(pickle_path + 'gene_dict/',
                'gene_dict' + search_param + '.pickle'))

    """
    ICGC MUT data
    """

    genes_missed = []
    gene_ID = [] # Enseble gene ID for ICGC query
    unpickled_gene_list_param = [' '] # list of pickled genes
    ICGC_full_mutations_list = []
    print('Searching ICGC database for mutation data: \n')
    t = time.time() ; t_master = t

    ICGC_df = pd.DataFrame(columns = [
        'id'
        'geneName',
        'type',
        'mutation',
        'referenceGenomeAllele',
        'start',
        'end',
        'testedDonorCount',
        'affectedDonorCountTotal',
        ])

    Gene_ID_Library = pickle.load(open(os.path.join(data_path,
        'GeneName2EnsembleID.pickle'), "rb")).to_dict('records')

    Gene_Syn_Library = pickle.load(open(os.path.join(data_path,
        'GeneSyn2EnsembleID.pickle'), "rb")).to_dict('records')

    for gene_name in Gene_ID_Library:

        if gene_name['Gene'] in gene_list:

            gene_ID += [[gene_name['Gene'], gene_name['ENSEMBL_ID']]]

    genes_missed = list(set(list(set(
        [i[0] for i in gene_ID]))).symmetric_difference(set(gene_list)))

    if genes_missed != []:

        print('Genes not found initially\n') ; print(genes_missed)

        for gene_name in Gene_Syn_Library:

            if any(x in gene_name['Gene'] for x in genes_missed):

                gene_ID += [[[x for x in gene_name['Gene'] if \
                              x in genes_missed][0], gene_name['ENSEMBL_ID']]]

        genes_missed = list(set(list(set(
            [i[0] for i in gene_ID]))).symmetric_difference(set(gene_list)))

    if genes_missed != []:

        print('Genes not found') ; print(genes_missed)
        raise Exception('rename genes or gene IDs!')

    gene_ID.sort()

    for ID in range(len(gene_ID)):

        no_muts = False

        if os.path.isfile(os.path.join(pickle_path + 'ICGC_data_rep/',
            gene_ID[ID][0] + '_' + gene_ID[ID][1] + search_param + '.pickle')):

            if gene_ID[ID][0] not in unpickled_gene_list_param:

                print('Unpickling '
                     + gene_ID[ID][0]
                     + ' mutation data for '
                     + search_param_short)

                ICGC_df_params = pickle.load(
                    open(os.path.join(pickle_path + 'ICGC_data_rep/',
                                     (gene_ID[ID][0]
                                     + '_'
                                     + gene_ID[ID][1]
                                     + search_param
                                     + '.pickle')), "rb"))

                unpickled_gene_list_param += [gene_ID[ID][0]]
                ICGC_full_mutations_list += [ICGC_df_params]

        elif gene_ID[ID][0] not in unpickled_gene_list_param:

            if os.path.isfile(os.path.join(pickle_path + 'ICGC_data_mut/',
                                          (gene_ID[ID][0]
                                          + '_'
                                          + gene_ID[ID][1]
                                          + '_mut_data.pickle'))):

                print('Unpickling ' + gene_ID[ID][0] + ' ICGC mutation data')

                ICGC_df = pickle.load(
                    open(os.path.join(pickle_path + 'ICGC_data_mut/',
                                      (gene_ID[ID][0]
                                      + '_'
                                      + gene_ID[ID][1]
                                      + '_mut_data.pickle')), "rb"))

            else:

                ICGC_url_custom_pages = ("https://dcc.icgc.org"
                    + "/api/v1/genes/"
                    + gene_ID[ID][1]
                    + "/mutations?field=id&field=type&"
                    + "field=chromosome&field="
                    + "start&field=end&field=mutation&"
                    + "field=referenceGenomeAllele"
                    + "&field=affectedDonorCountTotal&"
                    + "field=testedDonorCount&from=1"
                    + "&size=100")

                with urllib.request.urlopen(ICGC_url_custom_pages) as url:

                    json_data_ICGC = json.loads(url.read().decode())

                if json_data_ICGC["hits"] == []:

                    print('No ICGC mutations for '
                          + str(gene_ID[ID][0])
                          + '_'
                          + str(gene_ID[ID][1])
                          + '\n', flush = True)

                    no_muts = True

                else:

                    no_muts = False
                    mut_total = json_data_ICGC["pagination"]["total"]

                    print("pickling "
                         + str(mut_total)
                         + " ICGC mutations for: "
                         + str(gene_ID[ID][0])
                         + ' '
                         + str(gene_ID[ID][1]), flush = True)

                    data_pages = json_data_ICGC["pagination"]["pages"]
                    data_muts_list_gene_ICGC = []
                    pbar = tqdm(total = mut_total)

                    for page in (range(0, data_pages * 100, 100)):

                        pbar.update(100)
                        data_from_url_page_number = str(page + 1)

                        ICGC_url_mut_data = ("https://dcc.icgc.org/"
                            + "api/v1/genes/"
                            + gene_ID[ID][1]
                            + "/mutations?field=id&field=type&"
                            + "field=chromosome&field="
                            + "start&field=end&field=mutation&"
                            + "field=referenceGenomeAllele"
                            + "&field=affectedDonorCountTotal&"
                            + "field=testedDonorCount&from="
                            + data_from_url_page_number
                            + "&size=100")

                        with urllib.request.urlopen(ICGC_url_mut_data) as url:

                            json_data_ICGC = json.loads(url.read().decode())

                        data_muts_list_gene_ICGC += json_data_ICGC["hits"]

                    pbar.close()
                    ICGC_df = pd.DataFrame(data_muts_list_gene_ICGC)
                    ICGC_df['geneName'] = gene_ID[ID][0]

                    ICGC_df = ICGC_df.reindex(
                                columns = [
                                          'id',
                                          'geneName',
                                          'type',
                                          'mutation',
                                          'referenceGenomeAllele',
                                          'start',
                                          'end',
                                          'testedDonorCount',
                                          'affectedDonorCountTotal'
                                          ])

                    print(' mutations download time from ICGC: '
                         + str(float((time.time() - t)))[:5]
                         + 's', flush = True)

                    t = time.time()

                    pd.to_pickle(ICGC_df,
                        os.path.join(pickle_path + 'ICGC_data_mut/',
                            (gene_ID[ID][0]
                            + '_'
                            + gene_ID[ID][1]
                            + '_mut_data.pickle')))

            if no_muts == False:

                print(("Pickling "
                      + gene_ID[ID][0]
                      + '_'
                      + gene_ID[ID][1]
                      + " mutation data for search params: "
                      + search_param_short), flush = True)

                pbar = tqdm(total = ICGC_df.shape[0]) 
                Mut_df_data = [] ; type_new_list = []

                for row in ICGC_df.to_dict('records'):

                    mut_location = domain_match(row['geneName'],
                                                row['start'],
                                                row['end'])
                    
                    Mut_df_data += [mut_location]
                    type_new = row['type']

                    if ' of <=200bp' in type_new:

                        type_new = type_new[:-11]

                    elif ' (>=2bp and <=200bp)' in row['type']:

                        ref_size = len(row['referenceGenomeAllele'])

                        mut_size = len(row['mutation']) - \
                                   len(row['referenceGenomeAllele']) - 1

                        if ref_size == mut_size:

                            type_new = 'multiple base substitution'

                        elif ref_size > mut_size:

                            type_new = \
                                'multiple base substitution and deletion'

                        elif ref_size < mut_size:

                            type_new = \
                                'multiple base substitution and insertion'

                    if mut_location[0] == 'Yes':

                        if 'single base substitution' in type_new:
                            
                            try:
                                
                                mut_type = mutation_effect(gene_ID[ID][0],
                                                           row['mutation'][0],
                                                           row['mutation'][-1],
                                                           row['start'])
    
                                type_new += mut_type
                                
                            except:
                                
                                pass

                        elif 'deletion' in type_new:

                            if len(row['mutation'][:-2]) % 3 != 0:

                                type_new = 'deletion (frame shift)'

                            else:

                                type_new = 'deletion (in frame)'

                        elif 'insertion' in type_new:

                            if len(row['mutation'][2:]) % 3 != 0:

                                type_new = 'insertion (frame shift)'

                            else:

                                type_new = 'insertion (in frame)'

                        elif 'multiple base substitution' in type_new:

                            if abs(ref_size - mut_size) % 3 != 0:

                                type_new += ' (frame shift)'

                            else:

                                type_new += ' (in frame)'

                    else:

                        type_new += ' (non coding)'

                    type_new_list += [type_new] ; pbar.update(1)

                pbar.close()
                
                ICGC_df['type'] = type_new_list
                
                Mut_df = pd.DataFrame(Mut_df_data, columns = [
                                                            'coding_mut',
                                                            'UTR5_mut',
                                                            'UTR3_mut',
                                                            'intron_mut',
                                                            'intron_w_rep_mut',
                                                            'rep_mut',
                                                            'UTR_mut',
                                                            'transcript_mut',
                                                            'up_flank_mut',
                                                            'down_flank_mut',
                                                            'intergene_mut',
                                                            'non_coding_mut'
                                                            ])

                ICGC_df_params = pd.merge(ICGC_df, Mut_df,
                    left_index = True, right_index = True).set_index('id')

                ICGC_df_params = ICGC_df_params[
                    ~ICGC_df_params.index.duplicated(keep = 'first')]

                unpickled_gene_list_param += [gene_ID[ID][0]]

                print(' mutation domain search/calculation time: '
                     + str(float((time.time() - t)))[:5]
                     + 's\n', flush = True)

                t = time.time()

                pd.to_pickle(ICGC_df_params,
                    os.path.join(pickle_path + 'ICGC_data_rep/',
                                (gene_ID[ID][0]
                                + '_'
                                + gene_ID[ID][1]
                                + search_param
                                + '.pickle')))

                ICGC_full_mutations_list += [ICGC_df_params]

    mut_dict = pd.concat(ICGC_full_mutations_list)

    mut_dict_ref = mut_dict.loc[~((mut_dict['coding_mut'] == 'No') & \
                                  (mut_dict['non_coding_mut'] == 'No')), :]

    print('mutation dictionary build time for ('
          + str(len(gene_list))
          + ' genes): '
          + str(float((time.time() - t_master) / 60))[:5]
          + 'min\n')

    # genes_missed = [x for x in gene_list if x \
    #                 not in mut_dict_ref["geneName"].tolist()]

    # if genes_missed != []:

    #     print('\nGenes not found') ; print(genes_missed)
    #     raise Exception('rename genes or gene IDs!')

    """
    StringDB results for genes
    """

    if os.path.isfile(os.path.join(pickle_path + 'StringDB/',
        'GO' + search_param + '.pickle')):

        GO_df = pickle.load(open(os.path.join(pickle_path + 'StringDB/',
                            'GO' + search_param + '.pickle'), "rb"))

        print('Unpickling Gene Ontology data')

    else:

        print('Pickling StringDB Gene Ontology data')

        GO_DB = pickle.load(open(os.path.join(data_path,
            'StringDB_GO.pickle'), "rb")).dropna()

        GO_df = pd.DataFrame(0, index = GO_DB.Gene_name.unique(),
            columns = GO_DB.GO_term_name.unique())

        for GO_ID in GO_DB.to_dict('records'):

            GO_df.at[GO_ID['Gene_name'], GO_ID['GO_term_name']] = 1

        genes_SDB_search = list(set(list(gene_list)) &
                               set(list(GO_df.index.values)))

        GO_df = GO_df.loc[genes_SDB_search, :]
        GO_df_total = GO_df.iloc[0:-2].sum()
        GO_df_total.name = 'Total'
        GO_df = GO_df.append(GO_df_total.transpose())

        Mean_GO_all = pickle.load(open(os.path.join(data_path,
            'Mean_GO_String.pickle'), "rb")).set_index('Gene_name')

        GO_df = pd.concat([GO_df, Mean_GO_all])
        GO_df = GO_df.loc[:, (GO_df == 0).sum(axis = 0) < len(GO_df.index) - 1]
        GO_df = GO_df.loc[:, (GO_df == 1).sum(axis = 0) > 1]
        mean_GO_repeats = pd.DataFrame((GO_df[:-2]).mean(axis = 0)).T
        mean_GO_repeats['Gene_name'] = 'Mean_GO_ID_genes_with_repeat'
        mean_GO_repeats = mean_GO_repeats.set_index('Gene_name')
        GO_df = pd.concat([GO_df, mean_GO_repeats], sort = False)
        GO_IDs = GO_df.columns.tolist()

        for GO_ID in GO_IDs:

            if (GO_df.at['Mean_GO_ID_genes_with_repeat', GO_ID] / 2) < \
                GO_df.at['Mean_GO_ID_every_gene', GO_ID]:

                GO_df = GO_df.drop(columns = [GO_ID])

            elif GO_df.at['Total', GO_ID] < math.ceil(len(gene_list) / 10):

                GO_df = GO_df.drop(columns = [GO_ID])

        pd.to_pickle(GO_df, os.path.join(pickle_path + 'StringDB/',
            'GO' + search_param + '.pickle'))

    """
    Protein_Atlas tissue location
    """

    if os.path.isfile(os.path.join(pickle_path + 'Protein_Atlas/',
        'Tissue' + search_param + '.pickle')):

        tis_loc_df = pickle.load(
            open(os.path.join(pickle_path + 'Protein_Atlas/',
                'Tissue' + search_param + '.pickle'), "rb"))

        print('Unpickling tissue data')

    else:

        print('Pickling Tissue Protein Atlas data')

        tis_PA = pickle.load(open(os.path.join(data_path,
            'TissueProteinAtlas.pickle'), "rb")).fillna(0)

        tis_loc_df = pd.concat([tis_PA.loc[tis_PA['Gene'].isin(gene_list)],
            tis_PA.loc[tis_PA['Gene'].isin(['Mean_expression_every_gene'])]])

        tis_loc_df = tis_loc_df.loc[:,
            (tis_loc_df == 0).sum(axis = 0) < len(tis_loc_df.index) - 1]

        mean_tis = pd.DataFrame((tis_loc_df[:-1]).mean(axis = 0)).T
        mean_tis['Gene'] = 'Mean_expression_genes_with_repeat'

        tis_loc_df = pd.concat([tis_loc_df,mean_tis],
            sort = False).set_index('Gene')

        pd.to_pickle(tis_loc_df, os.path.join(pickle_path + 'Protein_Atlas/',
            'Tissue' + search_param + '.pickle'))
        
    """
    Protein_Atlas cancer expression
    """

    if os.path.isfile(os.path.join(pickle_path + 'Protein_Atlas/',
        'Cancer' + search_param + '.pickle')):

        can_loc_df = pickle.load(
            open(os.path.join(pickle_path + 'Protein_Atlas/',
                'Cancer' + search_param + '.pickle'), "rb"))

        print('Unpickling Cancer data')

    else:

        print('Pickling Cancer Protein Atlas data')

        can_PA = pickle.load(open(os.path.join(data_path,
            'CancerProteinAtlas.pickle'), "rb")).fillna('')

        can_loc_df = pd.concat([can_PA.loc[can_PA['Gene'].isin(gene_list)],
            can_PA.loc[can_PA['Gene'].isin(['Mean_expression_every_gene'])]])

        can_loc_df = can_loc_df.loc[:,
            (can_loc_df == 0).sum(axis = 0) < len(can_loc_df.index) - 1]

        mean_can = pd.DataFrame((can_loc_df[:-1]).mean(axis = 0)).T
        mean_can['Gene'] = 'Mean_expression_genes_with_repeat'

        can_loc_df = pd.concat([can_loc_df, mean_can],
            sort = False).set_index('Gene')

        pd.to_pickle(can_loc_df, os.path.join(pickle_path + 'Protein_Atlas/',
            'Cancer' + search_param + '.pickle'))

    """
    Protein_Atlas cell location
    """

    if os.path.isfile(os.path.join(pickle_path + 'Protein_Atlas/',
        'Cell' + search_param + '.pickle')):

        cell_loc_df = pickle.load(
            open(os.path.join(pickle_path + 'Protein_Atlas/',
                'Cell' + search_param + '.pickle'), "rb"))

        print('Unpickling Cell data')

    else:

        print('Pickling Cell Protein Atlas data')

        cell_PA = pickle.load(open(os.path.join(data_path,
            'CellProteinAtlas.pickle'), "rb")).fillna('')

        cell_loc_df = pd.concat([cell_PA.loc[cell_PA['Gene'].isin(gene_list)],
            cell_PA.loc[cell_PA['ID'].isin([row[1] for row in gene_ID])],
            cell_PA.iloc[-1:]])

        cell_loc_df = cell_loc_df.drop_duplicates()
        cell_loc_df = cell_loc_df.drop('ID', 1)
        mean_cell = pd.DataFrame((cell_loc_df[:-1]).mean(axis = 0)).T
        mean_cell['Gene'] = 'mean_genes_with_repeat'

        cell_loc_df = pd.concat([cell_loc_df, mean_cell],
                        sort = False).set_index('Gene')

        pd.to_pickle(cell_loc_df, os.path.join(pickle_path + 'Protein_Atlas/',
            'Cell' + search_param + '.pickle'))

    """
    Protein_Atlas organelle location
    """

    if os.path.isfile(os.path.join(pickle_path + 'Protein_Atlas/',
        'Organelle' + search_param + '.pickle')):

        org_loc_df = pickle.load(
            open(os.path.join(pickle_path + 'Protein_Atlas/',
                'Organelle' + search_param + '.pickle'), "rb"))

        print('Unpickling Organelle data')

    else:

        print('Pickling Organelle Protein Atlas data')

        org_PA = pickle.load(open(os.path.join(data_path,
            'OrganelleProteinAtlas.pickle'), "rb")).fillna('')

        org_loc_df = pd.concat([org_PA.loc[org_PA['Gene'].isin(gene_list)],
            org_PA.loc[org_PA['ID'].isin([row[1] for row in gene_ID])],
            org_PA.iloc[-1:]])

        org_loc_df = org_loc_df.drop_duplicates()
        org_loc_df = org_loc_df.drop('ID', 1)
        mean_org = pd.DataFrame((org_loc_df[:-1]).mean(axis = 0)).T
        mean_org['Gene'] = 'mean_genes_with_repeat'

        org_loc_df = pd.concat([org_loc_df, mean_org],
                        sort = False).set_index('Gene')
        
        pd.to_pickle(org_loc_df, os.path.join(pickle_path + 'Protein_Atlas/',
            'Organelle' + search_param + '.pickle'))

    """
    Protein_Atlas cell image creation
    """

    if os.path.isfile(os.path.join(im_path,
        'Cell_image' + search_param + '.png')):

        print('Opening cell image')
        avg = Image.open(im_path + 'Cell_image' + search_param + '.png')
        norm = Image.open(im_path + 'Cell_normal.png')
        
    else:
        
        try:
            
            print('Creating cell image')
            norm_vals = (org_loc_df.iloc[-1] / org_loc_df.iloc[-2]).values.tolist()
    
            norm_vals = (norm_vals - np.min(norm_vals)) / \
                        (np.max(norm_vals) - np.min(norm_vals))
    
            nucleoplasm_new = max([norm_vals[4],norm_vals[5]]) # kinetochore, mitotic chromosome
    
            if nucleoplasm_new > norm_vals[8]:
    
                norm_vals[8] = nucleoplasm_new
    
            microtubules_new = max([norm_vals[15],norm_vals[18], # cytokinetic bridge, mitotic spindle
                                    norm_vals[19],norm_vals[20]]) # midbody, midbody ring
    
            if microtubules_new > norm_vals[17]:
    
                norm_vals[17] = microtubules_new
    
            cytosol_new = max([norm_vals[21],norm_vals[22]]) # aggresome, cytoplasmic bodies
    
            if cytosol_new > norm_vals[23]:
    
                norm_vals[23] = cytosol_new
    
            if norm_vals[10] > norm_vals[9]: # cleavage furrow
    
                norm_vals[9] = norm_vals[10]
    
            if norm_vals[-8] > norm_vals[-7]: # cell junction
    
                norm_vals[-7] = norm_vals[-8]
    
            if norm_vals[-2] > norm_vals[-3]: # vesicles
    
                norm_vals[-3] = norm_vals[-2]
    
            if norm_vals[-2] > norm_vals[-4]:
    
                norm_vals[-4] = norm_vals[-2]
    
            if norm_vals[-2] > norm_vals[-5]:
    
                norm_vals[-5] = norm_vals[-2]
    
            if norm_vals[-2] > norm_vals[-6]:
    
                norm_vals[-6] = norm_vals[-2]
    
            org_weights = [norm_vals[23], norm_vals[12], norm_vals[9],
                           norm_vals[11], norm_vals[13], norm_vals[14],
                           norm_vals[17], norm_vals[16], norm_vals[-1],
                           norm_vals[-5], norm_vals[-4], norm_vals[-3],
                           norm_vals[-6], norm_vals[-10], norm_vals[-9],
                           norm_vals[8], norm_vals[0], norm_vals[6],
                           norm_vals[7], norm_vals[1], norm_vals[3],
                           norm_vals[2], norm_vals[24], norm_vals[25],
                           norm_vals[-7]]
    
            allfiles = os.listdir(im_path + 'Custom_cell')
            imlist = [file for file in allfiles if file[-4:] in ['.png']]
            imlist.sort() ; N = len(imlist)
            avg = Image.open(im_path + 'Custom_cell/' + imlist[0])
    
            for org in range(1, N):
    
                img = Image.open(im_path + 'Custom_cell/' + imlist[org])
                bands = list(img.split())
    
                if len(bands) == 4:
    
                    bands[3] = bands[3].point(lambda x: x* org_weights[org - 1])
    
                img = Image.merge(img.mode, bands) ; avg.paste(img, (0, 0), img)
    
            norm = Image.open('../data/Cell_images/Cell_normal.png')
            avg.save('../data/Cell_images/Cell_image' + search_param + '.png')
        
        except:
            
            pass

    """
    Calculating heat map
    """

    if os.path.isfile(os.path.join(pickle_path + 'mut_type/',
            'Mut_Type_all' + search_param + '.pickle')) and os.path.isfile(
        os.path.join(pickle_path + 'heat_map/',
            'Heat_Map' + search_param + '.pickle')) and os.path.isfile(
        os.path.join(pickle_path + 'mut_type/',
            'Mut_Type_general_all' + search_param + '.pickle')):

        print('Unpickling mutation type data (all)')

        mut_type_df = pickle.load(open(os.path.join(pickle_path + 'mut_type/',
            'Mut_Type_all' + search_param + '.pickle'), "rb"))

        print('Unpickling mutation type data (general)')

        general_df = pickle.load(open(os.path.join(pickle_path + 'mut_type/',
            'Mut_Type_general_all' + search_param + '.pickle'), "rb"))

        print('Unpickling heat map data\n')

        heat_map_df = pickle.load(open(os.path.join(pickle_path + 'heat_map/',
            'Heat_Map' + search_param + '.pickle'), "rb"))

    else:       

        domain_names = ['Coding', "5'UTR", "3'UTR", 'All Introns', 'Intron w/ repeat', 
                        'Repeat domain', 'All UTRs', 'Transcript',
                        'Up flank', 'Down flank', 'Intergene', 'Non coding']

        mut_names = ['coding_mut', 'UTR5_mut', 'UTR3_mut','intron_mut', 
                     'intron_w_rep_mut', 'rep_mut', 'UTR_mut', 'transcript_mut',
                     'up_flank_mut', 'down_flank_mut',
                     'intergene_mut', 'non_coding_mut']

        mut_general = ['single base substitution', 'insertion', 'deletion',
                       'multiple base substitution', 'in frame', 
                       'frame shifft']
        
        domain_num = len(domain_names) ; DomainLen_full = []

        heat_map_df = pd.DataFrame(index = gene_list,
            columns = (pd.MultiIndex.from_product(
                [['Mean tumor allele frequency',
                  'Mutation sites per BP',
                  'Impact score'],
                   domain_names])))

        mut_type_df = pd.DataFrame(0, index = gene_list,
                        columns = (pd.MultiIndex.from_product(
                            [mut_dict_ref.type.unique(), domain_names])))

        general_df = pd.DataFrame(0, index = gene_list,
                                columns = (pd.MultiIndex.from_product(
                                    [mut_general, domain_names])))
        

        print('\nGenerating data tables for '
              + str(len(gene_list))
              + ' genes:\n', flush = True) ; last = ''

        pbar = tqdm(total = len(gene_list))

        for gene in gene_dict.to_dict('records'):

            if gene['Gene'] != last:

                pbar.update(1)
                last = gene['Gene'] ; Mut_Pos_list = [] ; New_mut_pos = False
                all_df_list = [] ; mut_type_all_df = []
                Freq = np.zeros(domain_num)
                MutTotal = np.zeros(domain_num)
                MutSites = np.zeros(domain_num)
                SBSTotal = np.zeros(domain_num)
                INSTotal = np.zeros(domain_num)
                DELTotal = np.zeros(domain_num)
                MBSTotal = np.zeros(domain_num)

                DomainLen_gene = [gene["Coding_Length"],
                                  gene['UTR5_Length'],
                                  gene['UTR3_Length'],
                                  gene['Full_Intron_Length'],
                                  gene['Intron_W_Rep_Length'],
                                  gene['Full_Repeat_Length'],
                                  gene["Full_UTR_Length"],
                                  gene["Transcript_Length"],
                                  gene["Flank_Length"],
                                  gene["Flank_Length"],
                                  (gene["Flank_Length"] * 2),
                                  gene["Non_Coding_Length"]]

                DomainLen_full += [DomainLen_gene]

                mut_dict_ref_gene = mut_dict_ref.loc[
                    (mut_dict_ref['geneName'] == gene['Gene'])]

                for domain in range(domain_num):

                    mut_type_all_df = pd.DataFrame((mut_dict_ref_gene.loc[
                        (mut_dict_ref_gene[mut_names[domain]] == 'Yes')]
                        )['type'].value_counts())

                    mut_type_all_df = mut_type_all_df.rename(
                        columns = {'type' : domain_names[domain]})

                    all_df_list += [mut_type_all_df]

                mut_type_all = pd.concat(all_df_list, axis = 1, sort = False)

                for mut in range(len(mut_type_all.index)):

                    for domain in range(domain_num):

                        mut_type_df.loc[gene['Gene'],
                            (mut_type_all.index[mut], domain_names[domain])] \
                                = mut_type_all.loc[mut_type_all.index[mut],
                                                   domain_names[domain]]

            for mut in mut_dict_ref_gene.values:

                if mut[4] not in Mut_Pos_list or mut[5] not in Mut_Pos_list:

                    Mut_Pos_list += ([mut[4]] + [mut[5]])
                    Mut_Pos_list = list(set(Mut_Pos_list)) ; New_mut_pos = True

                else:

                    New_mut_pos = False

                for domain in range(domain_num):

                    if mut[8 + domain] == 'Yes':

                        if New_mut_pos is True:

                            MutSites[domain] += 1

                        MutTotal[domain] += 1
                        Freq[domain] += (int(mut[7]) / int(mut[6]))

                        if 'single base substitution' in mut[1]:

                            SBSTotal[domain] += 1

                        elif 'insertion' in mut[1]:

                            INSTotal[domain] += 1

                        elif 'deletion' in mut[1]:

                            DELTotal[domain] += 1

                        elif 'multiple base substitution' in mut[1]:

                            MBSTotal[domain] += 1

            gene_values = np.zeros(domain_num * 3)
            mut_values = np.zeros(domain_num * 6)

            for domain in range(domain_num):

                if not MutTotal[domain] == 0:

                    gene_values[domain] = Freq[domain] / MutTotal[domain]

                    gene_values[domain_num + domain] = \
                        MutSites[domain] / DomainLen_gene[domain]

                    gene_values[(domain_num * 2) + domain] = \
                        gene_values[domain] * gene_values[domain_num + domain]

                    mut_values[domain] = SBSTotal[domain]
                    mut_values[domain_num + domain] = INSTotal[domain]
                    mut_values[(domain_num * 2) + domain] = DELTotal[domain]
                    mut_values[(domain_num * 3) + domain] = MBSTotal[domain]

            heat_map_df.loc[gene['Gene']] = gene_values
            general_df.loc[gene['Gene']] = mut_values
            
        pbar.close()

        heat_map_df = heat_map_df.sort_values(
            [('Impact score', 'Repeat domain')],
                ascending = False).fillna(0).round(6)

        mean_len = [float(sum(col)) / len(col) for col in zip(*DomainLen_full)]
        heat_map_df.loc['mean'] = heat_map_df.mean()
        mut_type_df.loc['mean mut type count'] = mut_type_df.mean()

        mut_type_df.loc['mean domain length (BP)']  = mean_len * \
            int(len(mut_type_df.columns) / domain_num)

        mut_type_df.loc['mean mut type count per BP'] = \
            mut_type_df.loc['mean mut type count'] / \
            mut_type_df.loc['mean domain length (BP)']

        general_df.loc['mean mut type count'] = general_df.mean()

        general_df.loc['mean domain length (BP)']  = mean_len * \
            int(len(general_df.columns) / domain_num)

        general_df.loc['mean mut type count per BP'] = \
            general_df.loc['mean mut type count'] / \
            general_df.loc['mean domain length (BP)']

        print('Pickling mutation type data')

        pd.to_pickle(mut_type_df, os.path.join(pickle_path + 'mut_type/',
            'Mut_Type_all' + search_param + '.pickle'))

        print('Pickling general mutation type data')

        pd.to_pickle(general_df, os.path.join(pickle_path + 'mut_type/',
            'Mut_Type_general_all' + search_param + '.pickle'))

        print('Pickling heat map data\n')

        pd.to_pickle(heat_map_df, os.path.join(pickle_path + 'heat_map/',
            'Heat_Map' + search_param + '.pickle'))

    frequency_df = heat_map_df['Mean tumor allele frequency']
    mutations_per_BP_df = heat_map_df['Mutation sites per BP']
    impact_score_df = heat_map_df['Impact score']
    mut_type_df_colnames = list(mut_type_df.columns.levels[0])
    general_df_colnames = list(general_df.columns.levels[0])
    non_coding_colnames = [] ; coding_colnames = []

    for colname in mut_type_df_colnames:

        if 'non coding' in colname:

            non_coding_colnames += [colname]

        else:

            coding_colnames += [colname]

    mut_type_non_coding_df = mut_type_df[non_coding_colnames]
    mut_type_intron_df = mut_type_df[coding_colnames]
    mut_type_general_non_coding_df = general_df[general_df_colnames]
    mut_type_general_intron_df = general_df[general_df_colnames]

    mut_type_non_coding_df = mut_type_non_coding_df.loc[:,pd.IndexSlice[:,[
        "5'UTR", "3'UTR", 'All Introns', 'All UTRs', 'Up flank', 'Down flank',
        'Intergene', 'Non coding']]]
    
    mut_type_general_non_coding_df = mut_type_general_non_coding_df.loc \
        [:,pd.IndexSlice[:,[
        "5'UTR", "3'UTR", 'All Introns', 'All UTRs', 'Up flank', 'Down flank',
        'Intergene', 'Non coding']]]

    mut_type_intron_df = mut_type_intron_df.loc[:,pd.IndexSlice[:,[
        'All Introns', 'Intron w/ repeat', 'Repeat domain']]]
    
    mut_type_general_intron_df = mut_type_general_intron_df.loc \
        [:,pd.IndexSlice[:,['All Introns', 'Intron w/ repeat', 'Repeat domain']]]
        

    """
    Generating final output datatable CSVs
    """

    final_dfs_list = [gene_dict, mut_dict_ref, heat_map_df, frequency_df,
        mutations_per_BP_df, impact_score_df, mut_type_df, tis_loc_df,
        cell_loc_df, org_loc_df, can_loc_df, GO_df, mut_type_intron_df,
        mut_type_non_coding_df, general_df, mut_type_general_intron_df,
        mut_type_general_non_coding_df]

    final_dfs_names = ['gene_dict', 'mut_dict', 'mut_heat_map', 'frequency',
        'mutations_per_BP', 'impact_score', 'mut_type_all', 'tissue_location',
        'cell_location', 'organelle_location','cancer_expression',
        'GO_IDs_enriched', 'mut_type_coding', 'mut_type_non_coding',
        'mut_type_general', 'mut_type_general_intron_df', 
        'mut_type_general_non_coding_df'] ; i = 0
    
    for df in final_dfs_list:

        print('Creating CSV for ' + final_dfs_names[i] + ' dataset')

        df.to_csv(out_path + final_dfs_names[i] + search_param + '.csv',
            index = True) ; i += 1

    #Image.fromarray(np.hstack((np.array(norm), np.array(avg)))).show()
    print('\nmean: ' + str(impact_score_df.iloc[:-1].mean()[3]))
    print('gene_count: ' + str(len(gene_list)))
    print('stdev: ' + str(impact_score_df.iloc[:-1].std()[3]))
    
    print('\nFull execution time: '
         + str(float((time.time() - t_full) / 60))[:5]
         + 'min\n')

    return final_dfs_list

    return [mut_type_general_intron_df.tail(1).values.tolist(),
            [impact_score_df.iloc[:-1].mean()[5],
             len(gene_list),
             impact_score_df.iloc[:-1].std()[5]],
            cell_loc_df.tail(1).values.tolist()]

 

x = rep_search(bioc_mol_type = "Nucleotide",
                        reps = ['GAT'],
                        min_size = 21,
                        max_size = 1000,
                        gap_allowance = 0,
                        gene_flank_size = 500,
                        region = 'intron')



# codons = [['AAA'],['AAT'],['AAG'],['AAC'],
#  ['ATA'],['ATT'],['ATG'],['ATC'],
#  ['AGA'],['AGT'],['AGG'],['AGC'],
#  ['ACA'],['ACT'],['ACG'],['ACC'],
#  ['GGA'],['GGT'],['GGG'],['GGC'],
#  ['GTA'],['GTT'],['GTG'],['GTC'],
#  ['GAA'],['GAT'],['GAG'],['GAC'],
#  ['GCA'],['GCT'],['GCG'],['GCC'],
#  ['CCA'],['CCT'],['CCG'],['CCC'],
#  ['CTA'],['CTT'],['CTG'],['CTC'],
#  ['CGA'],['CGT'],['CGG'],['CGC'],
#  ['CAA'],['CAT'],['CAG'],['CAC'],
#  ['TTA'],['TTT'],['TTG'],['TTC'],
#  ['TAA'],['TAT'],['TAG'],['TAC'],
#  ['TGA'],['TGT'],['TGG'],['TGC'],
#  ['TCA'],['TCT'],['TCG'],['TCC']]

# data = []

# for codon in codons:
    
#     max_len = 1 ; max_achieved = False

#     while max_achieved is False:

#         try:
                
#             x = rep_search(bioc_mol_type = "Nucleotide",
#                             reps = codon,
#                             min_size = max_len,
#                             max_size = 10000,
#                             gap_allowance = 0,
#                             gene_flank_size = 500,
#                             region = 'intron')
            
#             max_len += 1
        
#         except:
            
#             print (codon,max_len-1)
#             data += [[codon,max_len-1]]
#             max_achieved = True
            
# #             pass

# min_len = 4 ; max_len = 5 ; y = [] ; y1 = [] ; y2 = []

# for n in [int(item) for item in np.linspace(min_len, max_len, max_len - min_len + 1)]:

#     try:
        
#         x = rep_search(bioc_mol_type = "Nucleotide",
#                         reps = ['CAG'],
#                         min_size = n,
#                         max_size = 1000,
#                         gap_allowance = 0,
#                         gene_flank_size = 500,
#                         region = 'intron')
        
#         y += x[0]
#         y1 += [x[1]]
#         y2 += [x[2]]
    
#     except:
        
#         y += [np.zeros(18)]
#         y1 += [np.zeros(3)]
#         y2 += [np.zeros(4)]
        
#         print('\nNo genes found')
#         pass

# z = pd.DataFrame(y)
# z1 = pd.DataFrame(y1)
# z2 = pd.DataFrame(y2)