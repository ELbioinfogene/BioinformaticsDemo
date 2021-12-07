#Eric Larsen
#12/3/2021
#Continued investigation of SARSCOV2 mutations

#Load functions from BioPython
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq

#FUNCTIONS:
def PARSE_SARS_2_GENOME(SEQ_RECORD):
    OUTPUT_DICT = {}
    #get accession Number
    ACC_ID = SEQ_RECORD.id
    ACC_ENTRY = {'Database Reference':ACC_ID}
    #isolate nucleotides
    NT_SEQ = SEQ_RECORD.seq
    #get total basepair length
    BP_TOTAL = len(NT_SEQ)
    SIZE_ENTRY = {'Genome Size':BP_TOTAL}
    #isolate each gene sequence
    #source: https://www.ncbi.nlm.nih.gov/gene/?term=NC_045512
    ORF1 = NT_SEQ[200:22000]
    ORF1_ENTRY = {'ORF1':ORF1}
    
    S_gene = NT_SEQ[21000:26000]
    S_ENTRY = {'Spike GlycoProtein':S_gene}
    
    E_Gene = NT_SEQ[26000:26500]
    E_ENTRY = {'Envelope':E_Gene}
    
    ORF3 = NT_SEQ[25300:26500]
    ORF3_ENTRY = {'ORF3':ORF3}
    
    N_gene = NT_SEQ[28000:29600]
    N_ENTRY = {'Nucleocapsid':N_gene}
    
    M_gene = NT_SEQ[26000:27200]
    M_ENTRY = {'Membrane':M_gene}
    
    ORF6 = NT_SEQ[27000:27500]
    ORF6_ENTRY = {'ORF6':ORF6}
    
    ORF7 = NT_SEQ[27200:28000]
    ORF7_ENTRY = {'ORF7':ORF7}
    
    ORF10 = NT_SEQ[29000::]
    ORF10_ENTRY = {'ORF10':ORF10}

    #Build OUTPUT_DICT
    OUTPUT_DICT.update(ACC_ENTRY)
    #OUTPUT_DICT.update(SIZE_ENTRY)
    OUTPUT_DICT.update(ORF1_ENTRY)
    OUTPUT_DICT.update(S_ENTRY)
    OUTPUT_DICT.update(E_ENTRY)
    OUTPUT_DICT.update(ORF3_ENTRY)
    OUTPUT_DICT.update(N_ENTRY)
    OUTPUT_DICT.update(M_ENTRY)
    OUTPUT_DICT.update(ORF6_ENTRY)
    OUTPUT_DICT.update(ORF7_ENTRY)
    OUTPUT_DICT.update(ORF10_ENTRY)
    return OUTPUT_DICT

def OUTPUT_GENE_FASTAS(SARS_COV2_DICT,NAME):
    ID_STRING=SARS_COV2_DICT['Database Reference']
    GENE_LIST=['Spike GlycoProtein','Envelope','Nucleocapsid','Membrane']
    for key in SARS_COV2_DICT:
        if key in GENE_LIST:
            GENE_NT = SARS_COV2_DICT[key]
            RECORD = SeqRecord.SeqRecord(GENE_NT)
            RECORD.id = ID_STRING
            RECORD.name = NAME
            RECORD.description = key
            SeqIO.write(RECORD,f'{NAME}-{key}.fasta','fasta')
    
#END OF FUNCTIONS

#FUNCTION_CALLS
#READ GENOME FILES and PARSE INDIVIDUAL GENES

#Original SARS-COV2
'''MN908947.3
    Collected: December 2019
    Location: China [Shanghai/Wuhan]'''
ALPHA_COVID_record = SeqIO.read('biopythonSARSCOV2.fasta', 'fasta')
'''NOTES:
    This is the SUBJECT of BLASTing while variants are the QUERY'''
ALPHA_GENES=PARSE_SARS_2_GENOME(ALPHA_COVID_record)
#FASTA acquired during inital alarm about strains

#'DELTA' variant
'''MZ577519.1
    Collected: July 12 2021
    Location:Louisiana USA'''
DELTA_COVID_record = SeqIO.read('biopythonSARSCOV2DELTAV.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: all single point mutations- looks natural
    ~~~
    Envelope: Identical to ALPHA
    Membrane: single point mutation at 634
    Nucleocapsid: Some point mutations, one Triplicate change is shared with OMICRON at 870
        SMALL INDEL at 270'''
DELTA_GENES=PARSE_SARS_2_GENOME(DELTA_COVID_record)
#FASTA acquired during inital alarm about strains

#'MU' variant
'''OK023693.1
    Collected: July 29 2021
    Location:Florida USA'''
MU_COVID_record = SeqIO.read('biopythonSARSCOV2MaybeMuV.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: point mutations and a 6 BP deletion in the high 900s
        An NGG is upstream of the 6 BP deletion, but that is nothing suspect
    ~~~
    Envelope: Identical to ALPHA
    Membrane: single point mutation at 731
    Nucleocapsid: a handful of point mutations
    '''
MU_GENES=PARSE_SARS_2_GENOME(MU_COVID_record)
#FASTA acquired during inital alarm about strains

#'OMICRON' variant
'''OL672836.1
    Collected: November 24 2021
    Location: Belgium'''
OMICRON_COVID_record = SeqIO.read('biopythonSARSCOV2OMICRON.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: MANY point mutations - including a Double Nucleotide Change
        2 deletion regions - one in the 600s and another around 900
        !An INDEL region around 1100!
        !Both deletions and the INDEL have CCN motifs upsteam and NGG downstream!
    ~~~
    Envelope: SINGLE POINT MUTATION AT 195 - !btwn a CGG and an AGG!
    Membrane: MULTIPLE POINT MUTATIONS
        !ALL OF THEM ARE Cas9 PAM ADJACENT!
    Nucleocapsid: alt codon? at 797 shared with DELTA
        !9 BP DELETION in PAM ENRICHED REGION!
        '''
OMICRON_GENES=PARSE_SARS_2_GENOME(OMICRON_COVID_record)
#Conclusion - THIS IS CAUSE FOR CONCERN

#MORE GENOME RECORDS FOR ADDITIONAL REFERENCE:
'''OV018878.1
    Collected: November 10 2021
    Location: UK'''
OTHER1_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_1.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: most point mutations in common with MU - including the 6 BP deletion
        handful of unique point mutants including one at 846 that is shared by OMICRON
    Envelope: Identical to ALPHA
    Membrane: Point mutant in common with MU
    Nucleocapsid: !UNIQUE DELETION at 240!
        Point Mutations shared with MU except for one at 903 shared with OL627097.1 [NJ<->UK TRAVEL]'''
OTHER1_GENES=PARSE_SARS_2_GENOME(OTHER1_COVID_record)
#evidence of active recombination

'''OL651917.1
    Collected:November 8 2021
    Location:New Mexico USA'''
OTHER2_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_2.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: Identical to ALPHA - except for some NNNs near the 3' (noise)
    ~~~
    Envelope: Indefinite NNNN regions - but identical to ALPHA
    Membrane: Indefinite NNNN regions - but identical to ALPHA
    Nucleocapsid: Indefinite NNNN regions - but identical to ALPHA'''
OTHER2_GENES=PARSE_SARS_2_GENOME(OTHER2_COVID_record)
#Alpha variant still in the wild - interesting

'''OL627097.1
    Collected: November 8 2021
    Location:New Jersey USA'''
OTHER3_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_3.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: similar to 'MU' - some point mutants
    ~~~
    Envelope: Identical to ALPHA
    Membrane: Point mutant in common with MU
    Nucleocapsid: Point mutant in common with MU  - one point is shared only with OV018878.1 [NJ<->UK TRAVEL]'''
OTHER3_GENES=PARSE_SARS_2_GENOME(OTHER3_COVID_record)
#evidence of active recombination

'''OV054768.1
    Collected: November 2021
    Location:Germany [Greifswald]'''
OTHER4_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_4.fasta', 'fasta')
'''NOTES:
    Spike GlycoProtein: Single point mutations - most in common with DELTA
        2 deletion sites: one is 6 BP around 756, the other is 3 BP around 980
        !BOTH ARE IN CCN and NGG enriched areas!
    ~~~
    Envelope: Identical to ALPHA
    Membrane: Identical to ALPHA
    Nucleocapsid: Point mutations and an INDEL very similar to DELTA'''
OTHER4_GENES=PARSE_SARS_2_GENOME(OTHER4_COVID_record)
#evidence of active recombination

'''
#Output all genes to .FASTA files for offline alignment in ApE
print('Writing Genes to files...')
OUTPUT_GENE_FASTAS(ALPHA_GENES,'ALPHA')
OUTPUT_GENE_FASTAS(DELTA_GENES,'DELTA')
OUTPUT_GENE_FASTAS(MU_GENES,'MU')
OUTPUT_GENE_FASTAS(OMICRON_GENES,'OMICRON')
OUTPUT_GENE_FASTAS(OTHER1_GENES,'ADDTL1')
OUTPUT_GENE_FASTAS(OTHER2_GENES,'ADDTL2')
OUTPUT_GENE_FASTAS(OTHER3_GENES,'ADDTL3')
OUTPUT_GENE_FASTAS(OTHER4_GENES,'ADDTL4')
print('..Done!')'''
