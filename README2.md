# MHC_slot_dataset Instructions

mhc_slot_dataset is a binding slot data parser that looks to extract the meaningful region of a given protein that is suspected to participate in bindinging activity. 
There is need for an automated function to do this task for large protein datasets that can have thousands of allelic variants. These alleles are not neccasirly the same size
and may not even contatin the binding region at all if they are a spliced variant. Additionaly the region that participates in binding might be of variable size depending on the
substrate of interest. This rules out using amino acid coordinates to loacte the rgion of interest. Instead we employ a sequence match function against a provided refrence 
sequnce that is known to participate in binding (must be entered as a hyper parameter). In addition to specifying hypermater values to optimize search for the binding region,
there are also 3 user inputs when calling the mhc_slot_dataset funtion: allele, Binding_model, and Return_type. 

Below are some examples of how to correctly use the MHC_slot_function:

```
mhc_slot_dataset(allele = ['HLA-A*01:01'], binding_model = 'linear', return_type = 'strings')
mhc_slot_dataset(allele = ['HLA-DRA*01:01:01:01','HLA-DRB1*01:01:01:01'], binding_model = 'cyclic', return_type = 'smiles')
mhc_slot_dataset(allele = ['QJE37815_surface_glycoprotein'], binding_model = 'binding_slot', return_type = 'smiles')
```

## Hyper parameters
Currently there are 3 hyper parameters: BINDING_DOMAIN_REF, THRESHOLD_COVERAGE, G_LINKER_LENGTH 

* BINDING_DOMAIN_REF is a known binding sequence for a protein of interest. Elements of this binding refrence sequence should be shared across alleles. Right now this is a manual 
task on the end of the user to enter this refrence sequnce by looking at literature or looking at RCSB as a hyperparameter. Future iterations will add a directory of known binding
sequences (zinc finger domains, bromo domains etc.) that can be entered in as part of the function input instead of as a hyper parameter. 

* THRESHOLD_COVERAGE is a value between 0 and 1 that denotes the minimum percentage match between the allele of interest and the BINDING_DOMAIN_REF needed to be considered a 
functional binding domain. For higly variable domains a lower value might be recomended and for highly conserved domains a higher value might be recomended. Adjust accordingly.

* G_LINKER_LENGTH is a global variable only required for binding_models (see below). A G_LINKER_LENGTH = 4, is equivalent to -G-G-G-G-.

<img src="markdownmonstericon.png"
     alt="Markdown Monster icon"
     style="float: left; margin-right: 10px;" />

```
BINDING_DOMAIN_REF = [
  "GRFASFEAQGALANIAVDKANLEIMTKRSNYTPITNVPPEVTVLTNSPVELR", # Refrence Sequnece helix 1 for HLA-DRA1 
  "TELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVG" # Refrence Sequnece helix 2 for HLA-DRB1
  ] 

THRESHOLD_COVERAGE = 0.30 # percent alignment for a succesful match (0-1 scale)

G_LINKER_LENGTH = 4 # varies depending on size of binding slot of the protein
```



## Binding_model
Currently there are 4 models you can chose from that conatin varying degrees of information as it pertains to the binding slot. The linear and cyclic models are the most
simplified interpretation of the binding slot imagining the binding slot as binding regions linked together with G linkers.
*  full_complex - full protein complex 
*  binding_slot - isolated binding slot, starts from first aligned position and continues to last aligned position. Contains binding regions and everything else 
*  linear - multiple disparate binding regionslinked in a linear manner
*  cyclic - multiple disparate binding regions linked in a circular manner

## MHC allele naming convention
MHC allele is a list of strings (species-id*0n:0n, or species-id*0n:0n:0n:0n's', i.e HLA-A*02:06 or HLA-A*01:01:38L) that follow a specific naming convention: 
  * species prefix - (HLA, BoLA, H2 etc)
  * gene id - (A,B,DRB1 etc)
  * allele group / broad protein variant - (04, 11, 100 etc)
  * allele subtype / SNPs - (04, 11, 100 etc) 
  * silent mutations - (04, 11, 100 etc) can omit from entry
  * non-coding mutations - (04, 11, 100 etc) can omit from entry
  * expression status - (Q, N, L) can omit from entry
   
# note that some MHC alleles given in csv don't follow this naming convention (rat alleles), however all HLA will adhere to nomenclature

