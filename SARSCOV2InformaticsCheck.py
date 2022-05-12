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
    ~~~
    Envelope: SINGLE POINT MUTATION AT 195
    Membrane: MULTIPLE POINT MUTATIONS
    Nucleocapsid: alt codon? at 797 shared with DELTA
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
    ~~~
    Envelope: Identical to ALPHA
    Membrane: Identical to ALPHA
    Nucleocapsid: Point mutations and an INDEL very similar to DELTA'''
OTHER4_GENES=PARSE_SARS_2_GENOME(OTHER4_COVID_record)
#evidence of active recombination

'''MZ372258
    Collected: May 2021
    Location:Illinois USA'''
OTHER5_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_5.fasta', 'fasta')
'''Earlier Sample of Omicron PangoLineage'''
OTHER5_GENES=PARSE_SARS_2_GENOME(OTHER5_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER5_GENES,'ADDTL5')

#5/12/2022 ADDITIONAL IMPACT - masking returns
#look at newer samples
'''ON464308.1
    Collected: April 2022
    Location:Colorado USA'''
OTHER6_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_6.fasta', 'fasta')
OTHER6_GENES=PARSE_SARS_2_GENOME(OTHER6_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER6_GENES,'ADDTL6')
#5/12/2022 - similar to ON480325

'''ON495673
    Collected: January 2022
    Location: Pennsylvania USA'''
OTHER7_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_7.fasta', 'fasta')
OTHER7_GENES=PARSE_SARS_2_GENOME(OTHER7_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER7_GENES,'ADDTL7')
#5/12/2022 - spike is very similar to omicron (drops a large insert?)
#other genes match omicron (not surprising)

'''ON480325
    Collected: February 2022
    Location: New Jersey USA'''
OTHER8_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_8.fasta', 'fasta')
OTHER8_GENES=PARSE_SARS_2_GENOME(OTHER8_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER8_GENES,'ADDTL8')
#5/12/2022 - spike has many features in common with omicron but some alpha recomb?
#membrane shows similar recomb, nucleocapsid is omicron + a SNP
#envelope is omicron
#this seems to be the 'new strain' that is spreading now

'''ON481333
    Collected: April 2022
    Location: Florida USA'''
OTHER9_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_9.fasta', 'fasta')
OTHER9_GENES=PARSE_SARS_2_GENOME(OTHER9_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER9_GENES,'ADDTL9')
#5/12/2022 - similar to ON480325

'''ON481762
    Collected: May 2022
    Location: California USA'''
OTHER10_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_10.fasta', 'fasta')
OTHER10_GENES=PARSE_SARS_2_GENOME(OTHER10_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER10_GENES,'ADDTL10')
#5/12/2022 - similar to ON480325

'''ON481750
    Collected: May 2022
    Location: New Hampshire USA'''
OTHER11_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_11.fasta', 'fasta')
OTHER11_GENES=PARSE_SARS_2_GENOME(OTHER11_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER11_GENES,'ADDTL11')
#5/12/2022 - similar to ON480325

'''ON481757
    Collected: May 2022
    Location: Michigan USA'''
OTHER12_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_12.fasta', 'fasta')
OTHER12_GENES=PARSE_SARS_2_GENOME(OTHER12_COVID_record)
#OUTPUT_GENE_FASTAS(OTHER12_GENES,'ADDTL12')
#5/12/2022 - similar to ON480325

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
print('..Done!')
'''
