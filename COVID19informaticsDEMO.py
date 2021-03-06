#Eric Larsen 2021
#Side project elaborating on udemy course lessons
#https://www.udemy.com/course/bioinformatics-with-python

#Load functions from BioPython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio import pairwise2

#load genomes from FASTA files
#load canonical COVID genome
CN_COVID_record = SeqIO.read('biopythonSARSCOV2.fasta', 'fasta')
CN_COVID_NT = CN_COVID_record.seq
#isolate spike gene
CN_SPIKE_GENE = CN_COVID_NT[21000:26000]

#load 'DELTA'
LA_COVID_record = SeqIO.read('biopythonSARSCOV2DELTAV.fasta', 'fasta')
LA_COVID_NT = LA_COVID_record.seq
#isolate spike gene
LA_SPIKE_GENE = LA_COVID_NT[21000:26000]

#load 'Mu' - Florida sample
FL_COVID_record = SeqIO.read('biopythonSARSCOV2MaybeMuV.fasta', 'fasta')
FL_COVID_NT = FL_COVID_record.seq
#isolate spike gene
FL_SPIKE_GENE = FL_COVID_NT[21000:26000]

#load 'OMICRON' - Belgium sample
OM_COVID_record = SeqIO.read('biopythonSARSCOV2OMICRON.fasta', 'fasta')
OM_COVID_NT = OM_COVID_record.seq
#isolate spike gene - wider window!
OM_SPIKE_GENE = OM_COVID_NT[21000:26000]

#load 'OMICRONb' - New Jersey sample
OMb_COVID_record = SeqIO.read('biopythonSARSCOV2SUPPLEMENTAL_8.fasta', 'fasta')
OMb_COVID_NT = OMb_COVID_record.seq
#isolate spike gene - wider window!
OMb_SPIKE_GENE = OMb_COVID_NT[21000:26000]
#end File Loading

#FUNCTIONS:
#returns the length of a sequence (needed to sort list by seq length)
def SEQLENGTH(seq):
    return len(seq)

#Finds all Start Codons and creates a list of possible proteins
def Protein_Estimator(DNA_QUERY):
    #make input ALLCAPS string
    DNA_STRING = str(DNA_QUERY).upper()
    PROTEIN_LIST = []
    #find DNA start codon(s) 'ATG'
    START_CODON_CANDIDATES = nt_search(DNA_STRING,'ATG')
    OUTPUT_SIZE = len(START_CODON_CANDIDATES)
    #remove [0] entry which is 'ATG'
    START_CODON_COORDINATES = START_CODON_CANDIDATES[1:OUTPUT_SIZE]
    #FROM EACH DETECTED START CODON - BEGIN IN-FRAME READING
    for n, S in enumerate(START_CODON_COORDINATES):
        #n is the candidate #
        #S is the address in DNA_STRING of the start codon
        PROTEIN_LIST.append(CODON_READER(DNA_QUERY,S))
    #test list for 'None' entries and remove them
    while None in PROTEIN_LIST:
        PROTEIN_LIST.remove(None)
    return PROTEIN_LIST
#7/31 works well

#Create a protein sequence from a single start codon
def CODON_READER(DNA_SEQ,START_ADDRESS):
    TRUNCATE_SEQUENCE = DNA_SEQ[START_ADDRESS::]
    SEQ_LENGTH = len(TRUNCATE_SEQUENCE)
    PROTEIN_BUILD = ''
    codon_count = 0
    #7/31 codon framing looks good
    for N in range(0,SEQ_LENGTH,3):
        codon = TRUNCATE_SEQUENCE[N:N+3]
        amino_acid = codon.translate(stop_symbol='#')
        #TO DO - address psuedo codon warning
        if str(amino_acid) != '#':
            PROTEIN_BUILD = PROTEIN_BUILD + amino_acid
            codon_count = codon_count + 1
        else:
            #DEBUG:
            #print('Protein sequnce made. Length:',str(codon_count),' amino acids')
            return PROTEIN_BUILD
#7/31 works well
#TO DO - address psuedo codon warning

#Takes 2 lists of posible proteins and creates a report on sequence identity
def PROTEIN_COMPARER(P_LIST1,P_LIST2,CULL_N):
    #REMOVE ALL ENTRIES LESS THAN CULL_N AA long
    CULLED_LIST1 = []
    CULLED_LIST2 = []
    PROT_OUTPUT_LIST = []
    for N, SEQ in enumerate(P_LIST1):
        if len(SEQ) > CULL_N:
            CULLED_LIST1.append(SEQ)
    for N, SEQ in enumerate(P_LIST2):
        if len(SEQ) > CULL_N:
            CULLED_LIST2.append(SEQ)
    #7/31 culling works
    #test using set()
    CULLED_SET1 = set(CULLED_LIST1)
    CULLED_SET2 = set(CULLED_LIST2)

    #ONLY COMPARE SEQUENCES THAT ARE DIFFERENT!
    SET1_UNIQUES = list(CULLED_SET1.difference(CULLED_SET2))
    SET2_UNIQUES = list(CULLED_SET2.difference(CULLED_SET1))

    #SORT UNIQUE SEQUENCES BY LENGTH
    SET1_UNIQUES.sort(key=SEQLENGTH)
    SET2_UNIQUES.sort(key=SEQLENGTH)

    #8/11/21
    #Compare lengths of sets so the shorter one is used for enumerating the loop
    SET_SIZE_DIFFERENCE = len(SET1_UNIQUES) - len(SET2_UNIQUES)
    #SET1_UNIQUES is larger or the same size as SET2_UNIQUES
    if SET_SIZE_DIFFERENCE >= 0:
        UNIQUE_INDEX = enumerate(SET2_UNIQUES)
    #SET1_UNIQUES is smaller than SET2_UNIQUES
    if SET_SIZE_DIFFERENCE < 0:
        UNIQUE_INDEX = enumerate(SET1_UNIQUES)
    
    #go through all SET1 uniques and do a pairwise identity score with SET2 uniques
    for N, SEQ in UNIQUE_INDEX:
        SET1_SEQ = SET1_UNIQUES[N]
        SET2_SEQ = SET2_UNIQUES[N]
        #TO DO - WHAT IF SET2_Uniques is too small?
        #Fixed - 8/11/21
        AA1_count = len(SET1_SEQ)
        AA2_count = len(SET2_SEQ)
        SCORE = pairwise2.align.globalxx(SET1_SEQ,SET2_SEQ,one_alignment_only=True,score_only=True)
        AAdiff = AA1_count-AA2_count
        if AAdiff > 0:
            IDENTITY = (SCORE/AA2_count)*100
        if AAdiff <= 0:
            IDENTITY = (SCORE/AA1_count)*100
        SCORE_DICT = {'SEQA':SET1_SEQ,
                      'SEQB':SET2_SEQ,
                      'A_length':AA1_count,
                      'B_length':AA2_count,
                      'SCORE':SCORE,
                      'IDENTITY':IDENTITY}
        PROT_OUTPUT_LIST.append(SCORE_DICT)
    return PROT_OUTPUT_LIST
#7/31 seems to work!!

#8/11/21 additional function for displaying output of PROTEIN_COMPARER()
def DISPLAY_PROTEIN_COMPARISON(PROT_LIST,DISPLAY_INDEX,NAME_SET,ALIGNMENT_TYPE):
    #inputs are a list of compared amino acid sequences
    #(PROT_OUTPUT_LIST produced by function PROTEIN_COMPARER),
    #specific index integer, a list of 2 strings to name the proteins
    #and a string for the type of alignment (global or local)
    #8/11/21 TO DO - add global/local switch for PROTEIN_COMPARER
    #get names
    PROTEIN_A_NAME = NAME_SET[0]
    PROTEIN_B_NAME = NAME_SET[1]
    #get lengths as strings
    PROTEIN_A_LENGTH = str(PROT_LIST[DISPLAY_INDEX]['A_length'])
    PROTEIN_B_LENGTH = str(PROT_LIST[DISPLAY_INDEX]['B_length'])
    #get identity score as a string from a 4 digit floating
    SEQUENCE_IDENTITY_SCORE = str(round(PROT_LIST[DISPLAY_INDEX]['IDENTITY'],4))
    print(f'Protein {PROTEIN_A_NAME} is {PROTEIN_A_LENGTH} amino acids long.\n')
    print(f'Protein {PROTEIN_B_NAME} is {PROTEIN_B_LENGTH} amino acids long.\n')
    print(f'{PROTEIN_A_NAME} and {PROTEIN_B_NAME} have {SEQUENCE_IDENTITY_SCORE}% sequence identity in a {ALIGNMENT_TYPE} alignment\n')

def VAX_DRIFT_REPORT(VAX_SPIKE,QUERY_SPIKE,QUERY_NAME,SIG_FIGS):
    #which is longer?
    V_SIZE = len(VAX_SPIKE)
    Q_SIZE = len(QUERY_SPIKE)
    if V_SIZE>Q_SIZE:
        SIZE_NORM = V_SIZE
    if Q_SIZE>V_SIZE:
        SIZE_NORM = Q_SIZE
    if Q_SIZE==V_SIZE:
        SIZE_NORM = V_SIZE
    PAIRWISE_SCORE = pairwise2.align.globalxx(VAX_SPIKE,QUERY_SPIKE,one_alignment_only=True,score_only=True)
    IDENTITY_PERCENT = (PAIRWISE_SCORE/SIZE_NORM)*100
    ROUND_OFF_PERCENT = round(IDENTITY_PERCENT,SIG_FIGS)
    print(f'The Pfizer mRNA spike protein is {ROUND_OFF_PERCENT}% identical to {QUERY_NAME} variant spike protein')
    print('***\n')

def PROT_DRIFT_REPORT(PROTEIN_A,PROTEIN_B,NAME_LIST,SIG_FIGS):
    #which is longer?
    A_SIZE = len(PROTEIN_A)
    B_SIZE = len(PROTEIN_B)
    if A_SIZE>B_SIZE:
        SIZE_NORM = A_SIZE
    if B_SIZE>A_SIZE:
        SIZE_NORM = B_SIZE
    if A_SIZE==B_SIZE:
        SIZE_NORM = A_SIZE
    PAIRWISE_SCORE = pairwise2.align.globalxx(PROTEIN_A,PROTEIN_B,one_alignment_only=True,score_only=True)
    IDENTITY_PERCENT = (PAIRWISE_SCORE/SIZE_NORM)*100
    ROUND_OFF_PERCENT = round(IDENTITY_PERCENT,SIG_FIGS)
    A_NAME = NAME_LIST[0]
    B_NAME = NAME_LIST[1]
    PROT_NAME = NAME_LIST[2]
    print(f'The {PROT_NAME} AA sequence of {B_NAME} variant is {ROUND_OFF_PERCENT}% identical to the {A_NAME} variant.')
    print('===\n')

#END OF FUNCTIONS

#Function Execution
#Make Lists of Potential Proteins
CN_PROTEINS = Protein_Estimator(CN_SPIKE_GENE)
LA_PROTEINS = Protein_Estimator(LA_SPIKE_GENE)
FL_PROTEINS = Protein_Estimator(FL_SPIKE_GENE)
OM_PROTEINS = Protein_Estimator(OM_SPIKE_GENE)
OMb_PROTEINS = Protein_Estimator(OMb_SPIKE_GENE)
#Compare unique proteins of similar length - allows user to select a cull threshold
DELTA_LIST = PROTEIN_COMPARER(CN_PROTEINS,LA_PROTEINS,100)
MU_LIST = PROTEIN_COMPARER(CN_PROTEINS,FL_PROTEINS,100)
OM_LIST = PROTEIN_COMPARER(CN_PROTEINS,OM_PROTEINS,100)
OMb_LIST = PROTEIN_COMPARER(CN_PROTEINS,OMb_PROTEINS,100)

#stop using this function for now in favor of PROT_DRIFT_REPORT()
#DISPLAY_PROTEIN_COMPARISON(DELTA_LIST,10,['COVID19 ALPHA SPIKE','COVID19 DELTA SPIKE'],'global')

#Analyze VAX target mutations
#Pfizer Sequence from WHO
#Source: https://berthub.eu/articles/posts/reverse-engineering-source-code-of-the-biontech-pfizer-vaccine/
#used microsoft word to isolate the bolded reading frame, remove numbers and spaces, converted greek special base to Uracil for ease of use
PFIZER_WHO_VAX=Seq('AUGUUCGUGUUCCUGGUGCUGCUGCCUCUGGUGUCCAGCCAGUGUGUGAACCUGACCACCAGAACACAGCUGCCUCCAGCCUACACCAACAGCUUUACCAGAGGCGUGUACUACCCCGACAAGGUGUUCAGAUCCAGCGUGCUGCACUCUACCCAGGACCUGUUCCUGCCUUUCUUCAGCAACGUGACCUGGUUCCACGCCAUCCACGUGUCCGGCACCAAUGGCACCAAGAGAUUCGACAACCCCGUGCUGCCCUUCAACGACGGGGUGUACUUUGCCAGCACCGAGAAGUCCAACAUCAUCAGAGGCUGGAUCUUCGGCACCACACUGGACAGCAAGACCCAGAGCCUGCUGAUCGUGAACAACGCCACCAACGUGGUCAUCAAAGUGUGCGAGUUCCAGUUCUGCAACGACCCCUUCCUGGGCGUCUACUACCACAAGAACAACAAGAGCUGGAUGGAAAGCGAGUUCCGGGUGUACAGCAGCGCCAACAACUGCACCUUCGAGUACGUGUCCCAGCCUUUCCUGAUGGACCUGGAAGGCAAGCAGGGCAACUUCAAGAACCUGCGCGAGUUCGUGUUUAAGAACAUCGACGGCUACUUCAAGAUCUACAGCAAGCACACCCCUAUCAACCUCGUGCGGGAUCUGCCUCAGGGCUUCUCUGCUCUGGAACCCCUGGUGGAUCUGCCCAUCGGCAUCAACAUCACCCGGUUUCAGACACUGCUGGCCCUGCACAGAAGCUACCUGACACCUGGCGAUAGCAGCAGCGGAUGGACAGCUGGUGCCGCCGCUUACUAUGUGGGCUACCUGCAGCCUAGAACCUUCCUGCUGAAGUACAACGAGAACGGCACCAUCACCGACGCCGUGGAUUGUGCUCUGGAUCCUCUGAGCGAGACAAAGUGCACCCUGAAGUCCUUCACCGUGGAAAAGGGCAUCUACCAGACCAGCAACUUCCGGGUGCAGCCCACCGAAUCCAUCGUGCGGUUCCCCAAUAUCACCAAUCUGUGCCCCUUCGGCGAGGUGUUCAAUGCCACCAGAUUCGCCUCUGUGUACGCCUGGAACCGGAAGCGGAUCAGCAAUUGCGUGGCCGACUACUCCGUGCUGUACAACUCCGCCAGCUUCAGCACCUUCAAGUGCUACGGCGUGUCCCCUACCAAGCUGAACGACCUGUGCUUCACAAACGUGUACGCCGACAGCUUCGUGAUCCGGGGAGAUGAAGUGCGGCAGAUUGCCCCUGGACAGACAGGCAAGAUCGCCGACUACAACUACAAGCUGCCCGACGACUUCACCGGCUGUGUGAUUGCCUGGAACAGCAACAACCUGGACUCCAAAGUCGGCGGCAACUACAAUUACCUGUACCGGCUGUUCCGGAAGUCCAAUCUGAAGCCCUUCGAGCGGGACAUCUCCACCGAGAUCUAUCAGGCCGGCAGCACCCCUUGUAACGGCGUGGAAGGCUUCAACUGCUACUUCCCACUGCAGUCCUACGGCUUUCAGCCCACAAAUGGCGUGGGCUAUCAGCCCUACAGAGUGGUGGUGCUGAGCUUCGAACUGCUGCAUGCCCCUGCCACAGUGUGCGGCCCUAAGAAAAGCACCAAUCUCGUGAAGAACAAAUGCGUGAACUUCAACUUCAACGGCCUGACCGGCACCGGCGUGCUGACAGAGAGCAACAAGAAGUUCCUGCCAUUCCAGCAGUUUGGCCGGGAUAUCGCCGAUACCACAGACGCCGUUAGAGAUCCCCAGACACUGGAAAUCCUGGACAUCACCCCUUGCAGCUUCGGCGGAGUGUCUGUGAUCACCCCUGGCACCAACACCAGCAAUCAGGUGGCAGUGCUGUACCAGGACGUGAACUGUACCGAAGUGCCCGUGGCCAUUCACGCCGAUCAGCUGACACCUACAUGGCGGGUGUACUCCACCGGCAGCAAUGUGUUUCAGACCAGAGCCGGCUGUCUGAUCGGAGCCGAGCACGUGAACAAUAGCUACGAGUGCGACAUCCCCAUCGGCGCUGGAAUCUGCGCCAGCUACCAGACACAGACAAACAGCCCUCGGAGAGCCAGAAGCGUGGCCAGCCAGAGCAUCAUUGCCUACACAAUGUCUCUGGGCGCCGAGAACAGCGUGGCCUACUCCAACAACUCUAUCGCUAUCCCCACCAACUUCACCAUCAGCGUGACCACAGAGAUCCUGCCUGUGUCCAUGACCAAGACCAGCGUGGACUGCACCAUGUACAUCUGCGGCGAUUCCACCGAGUGCUCCAACCUGCUGCUGCAGUACGGCAGCUUCUGCACCCAGCUGAAUAGAGCCCUGACAGGGAUCGCCGUGGAACAGGACAAGAACACCCAAGAGGUGUUCGCCCAAGUGAAGCAGAUCUACAAGACCCCUCCUAUCAAGGACUUCGGCGGCUUCAAUUUCAGCCAGAUUCUGCCCGAUCCUAGCAAGCCCAGCAAGCGGAGCUUCAUCGAGGACCUGCUGUUCAACAAAGUGACACUGGCCGACGCCGGCUUCAUCAAGCAGUAUGGCGAUUGUCUGGGCGACAUUGCCGCCAGGGAUCUGAUUUGCGCCCAGAAGUUUAACGGACUGACAGUGCUGCCUCCUCUGCUGACCGAUGAGAUGAUCGCCCAGUACACAUCUGCCCUGCUGGCCGGCACAAUCACAAGCGGCUGGACAUUUGGAGCAGGCGCCGCUCUGCAGAUCCCCUUUGCUAUGCAGAUGGCCUACCGGUUCAACGGCAUCGGAGUGACCCAGAAUGUGCUGUACGAGAACCAGAAGCUGAUCGCCAACCAGUUCAACAGCGCCAUCGGCAAGAUCCAGGACAGCCUGAGCAGCACAGCAAGCGCCCUGGGAAAGCUGCAGGACGUGGUCAACCAGAAUGCCCAGGCACUGAACACCCUGGUCAAGCAGCUGUCCUCCAACUUCGGCGCCAUCAGCUCUGUGCUGAACGAUAUCCUGAGCAGACUGGACCCUCCUGAGGCCGAGGUGCAGAUCGACAGACUGAUCACAGGCAGACUGCAGAGCCUCCAGACAUACGUGACCCAGCAGCUGAUCAGAGCCGCCGAGAUUAGAGCCUCUGCCAAUCUGGCCGCCACCAAGAUGUCUGAGUGUGUGCUGGGCCAGAGCAAGAGAGUGGACUUUUGCGGCAAGGGCUACCACCUGAUGAGCUUCCCUCAGUCUGCCCCUCACGGCGUGGUGUUUCUGCACGUGACAUAUGUGCCCGCUCAAGAGAAGAAUUUCACCACCGCUCCAGCCAUCUGCCACGACGGCAAAGCCCACUUUCCUAGAGAAGGCGUGUUCGUGUCCAACGGCACCCAUUGGUUCGUGACACAGCGGAACUUCUACGAGCCCCAGAUCAUCACCACCGACAACACCUUCGUGUCUGGCAACUGCGACGUCGUGAUCGGCAUUGUGAACAAUACCGUGUACGACCCUCUGCAGCCCGAGCUGGACAGCUUCAAAGAGGAACUGGACAAGUACUUUAAGAACCACACAAGCCCCGACGUGGACCUGGGCGAUAUCAGCGGAAUCAAUGCCAGCGUCGUGAACAUCCAGAAAGAGAUCGACCGGCUGAACGAGGUGGCCAAGAAUCUGAACGAGAGCCUGAUCGACCUGCAAGAACUGGGGAAGUACGAGCAGUACAUCAAGUGGCCCUGGUACAUCUGGCUGGGCUUUAUCGCCGGACUGAUUGCCAUCGUGAUGGUCACAAUCAUGCUGUGUUGCAUGACCAGCUGCUGUAGCUGCCUGAAGGGCUGUUGUAGCUGUGGCAGCUGCUGCAAGUUCGACGAGGACGAUUCUGAGCCCGUGCUGAAGGGCGUGAAACUGCACUACACAUGAUGA')
#'Bespoke' solution:
PFIZER_SPIKE = PFIZER_WHO_VAX.translate()
VAX_AA_SIZE = len(PFIZER_SPIKE)
#ISOLATE 'CANON' spike proteins
CN_SPIKE_CANON = DELTA_LIST[10]['SEQA']
CN_AA_SIZE = len(CN_SPIKE_CANON)
LA_SPIKE_CANON = DELTA_LIST[10]['SEQB']
LA_AA_SIZE = len(LA_SPIKE_CANON)
#12/2/21 this matches reported Mu spike protein QZW69948.1
MU_SPIKE_CANON = MU_LIST[8]['SEQB']
MU_AA_SIZE = len(MU_SPIKE_CANON)
#12/2/21 - this matches reporteded Omicron spike protein UFO69279.1
OM_SPIKE_CANON = OM_LIST[8]['SEQB']
OM_AA_SIZE = len(OM_SPIKE_CANON)
#5/12/2022 - this matches reported Omicron 'b' spike protein UQI51632.1
OMb_SPIKE_CANON = OMb_LIST[10]['SEQB']
OMb_AA_SIZE = len(OMb_SPIKE_CANON)

#report Spike protein changes
PROT_DRIFT_REPORT(CN_SPIKE_CANON,LA_SPIKE_CANON,['Alpha','Delta','Spike'],3)
PROT_DRIFT_REPORT(CN_SPIKE_CANON,OM_SPIKE_CANON,['Alpha','Omicron','Spike'],3)
PROT_DRIFT_REPORT(CN_SPIKE_CANON,OMb_SPIKE_CANON,['Alpha','Omicron b','Spike'],3)

VAX_DRIFT_REPORT(PFIZER_SPIKE,CN_SPIKE_CANON,'Alpha',3)
VAX_DRIFT_REPORT(PFIZER_SPIKE,LA_SPIKE_CANON,'Delta',3)
VAX_DRIFT_REPORT(PFIZER_SPIKE,MU_SPIKE_CANON,'Mu',3)
VAX_DRIFT_REPORT(PFIZER_SPIKE,OM_SPIKE_CANON,'Omicron',3)
VAX_DRIFT_REPORT(PFIZER_SPIKE,OMb_SPIKE_CANON,'Omicron b',3)

#Look at genes besides the spike
#(not targeted by Vaccine, significance even less understood)
#source: https://www.ncbi.nlm.nih.gov/gene/?term=NC_045512
print('Analyzing other SARS-COV2 proteins...\n~~~~')

#Isolate gene sequences:
#ORF1ab - this is the largest reading frame
#always run Protein_Estimator() on this LAST as it will be SLOW
ALPH_ORF1_GENE = CN_COVID_NT[200:22000]
DELT_ORF1_GENE = LA_COVID_NT[200:22000]
OMIC_ORF1_GENE = OM_COVID_NT[200:22000]

#ORF3a
ALPH_ORF3_GENE = CN_COVID_NT[25300:26500]
DELT_ORF3_GENE = LA_COVID_NT[25300:26500]
OMIC_ORF3_GENE = OM_COVID_NT[25300:26500]

#Envelope
ALPH_ENVELOPE_GENE = CN_COVID_NT[26000:26500]
DELT_ENVELOPE_GENE = LA_COVID_NT[26000:26500]
OMIC_ENVELOPE_GENE = OM_COVID_NT[26000:26500]

#Membrane
ALPH_MEMBRANE_GENE = CN_COVID_NT[26000:27200]
DELT_MEMBRANE_GENE = LA_COVID_NT[26000:27200]
OMIC_MEMBRANE_GENE = OM_COVID_NT[26000:27200]

#Nucleocapsid
ALPH_NUCLEO_GENE = CN_COVID_NT[28000:29600]
DELT_NUCLEO_GENE = LA_COVID_NT[28000:29600]
OMIC_NUCLEO_GENE = OM_COVID_NT[28000:29600]

#ORF6
ALPH_ORF6_GENE = CN_COVID_NT[27000:27500]
DELT_ORF6_GENE = LA_COVID_NT[27000:27500]
OMIC_ORF6_GENE = OM_COVID_NT[27000:27500]

#ORF7
ALPH_ORF7_GENE = CN_COVID_NT[27200:28000]
DELT_ORF7_GENE = LA_COVID_NT[27200:28000]
OMIC_ORF7_GENE = OM_COVID_NT[27200:28000]

#ORF10 - this is the 'tail' of the genome
ALPH_ORF10_GENE = CN_COVID_NT[29000::]
DELT_ORF10_GENE = LA_COVID_NT[29000::]
OMIC_ORF10_GENE = OM_COVID_NT[29000::]
#End of gene collection

#Produce protein candidates and look at mutations in the 'true' AA seqs:
#membrane protein
ALPH_MEMBRANE_AAs = Protein_Estimator(ALPH_MEMBRANE_GENE)
DELT_MEMBRANE_AAs = Protein_Estimator(DELT_MEMBRANE_GENE)
OMIC_MEMBRANE_AAs = Protein_Estimator(OMIC_MEMBRANE_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_MEMBRANE_PROTEIN = ALPH_MEMBRANE_AAs[8]
DELT_MEMBRANE_PROTEIN = DELT_MEMBRANE_AAs[8]
OMIC_MEMBRANE_PROTEIN = OMIC_MEMBRANE_AAs[7]
#report results
PROT_DRIFT_REPORT(ALPH_MEMBRANE_PROTEIN,DELT_MEMBRANE_PROTEIN,['Alpha','Delta','Membrane'],3)
PROT_DRIFT_REPORT(ALPH_MEMBRANE_PROTEIN,OMIC_MEMBRANE_PROTEIN,['Alpha','Omicron','Membrane'],3)

#NucleoCapsid protein
ALPH_NUCLEO_AAs = Protein_Estimator(ALPH_NUCLEO_GENE)
DELT_NUCLEO_AAs = Protein_Estimator(DELT_NUCLEO_GENE)
OMIC_NUCLEO_AAs = Protein_Estimator(OMIC_NUCLEO_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_NUCLEO_PROTEIN = ALPH_NUCLEO_AAs[4]
DELT_NUCLEO_PROTEIN = DELT_NUCLEO_AAs[4]
OMIC_NUCLEO_PROTEIN = OMIC_NUCLEO_AAs[3]
#report results
PROT_DRIFT_REPORT(ALPH_NUCLEO_PROTEIN,DELT_NUCLEO_PROTEIN,['Alpha','Delta','Nucleocapsid'],3)
PROT_DRIFT_REPORT(ALPH_NUCLEO_PROTEIN,OMIC_NUCLEO_PROTEIN,['Alpha','Omicron','Nucleocapsid'],3)

#Envelope protein
ALPH_ENV_AAs = Protein_Estimator(ALPH_ENVELOPE_GENE)
DELT_ENV_AAs = Protein_Estimator(DELT_ENVELOPE_GENE)
OMIC_ENV_AAs = Protein_Estimator(OMIC_ENVELOPE_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_ENV_PROTEIN = ALPH_ENV_AAs[-1]
DELT_ENV_PROTEIN = DELT_ENV_AAs[-1]
OMIC_ENV_PROTEIN = OMIC_ENV_AAs[-1]
#report results
PROT_DRIFT_REPORT(ALPH_ENV_PROTEIN,DELT_ENV_PROTEIN,['Alpha','Delta','Envelope'],3)
PROT_DRIFT_REPORT(ALPH_ENV_PROTEIN,OMIC_ENV_PROTEIN,['Alpha','Omicron','Envelope'],3)

#ORF3
ALPH_ORF3_AAs = Protein_Estimator(ALPH_ORF3_GENE)
DELT_ORF3_AAs = Protein_Estimator(DELT_ORF3_GENE)
OMIC_ORF3_AAs = Protein_Estimator(OMIC_ORF3_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_ORF3_PROTEIN = ALPH_ORF3_AAs[1]
DELT_ORF3_PROTEIN = DELT_ORF3_AAs[1]
OMIC_ORF3_PROTEIN = OMIC_ORF3_AAs[0]
#report results
PROT_DRIFT_REPORT(ALPH_ORF3_PROTEIN,DELT_ORF3_PROTEIN,['Alpha','Delta','ORF3'],3)
PROT_DRIFT_REPORT(ALPH_ORF3_PROTEIN,OMIC_ORF3_PROTEIN,['Alpha','Omicron','ORF3'],3)

#ORF6
ALPH_ORF6_AAs = Protein_Estimator(ALPH_ORF6_GENE)
DELT_ORF6_AAs = Protein_Estimator(DELT_ORF6_GENE)
OMIC_ORF6_AAs = Protein_Estimator(OMIC_ORF6_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_ORF6_PROTEIN = ALPH_ORF6_AAs[0]
DELT_ORF6_PROTEIN = DELT_ORF6_AAs[0]
OMIC_ORF6_PROTEIN = OMIC_ORF6_AAs[0]
#report results
PROT_DRIFT_REPORT(ALPH_ORF6_PROTEIN,DELT_ORF6_PROTEIN,['Alpha','Delta','ORF6'],3)
PROT_DRIFT_REPORT(ALPH_ORF6_PROTEIN,OMIC_ORF6_PROTEIN,['Alpha','Omicron','ORF6'],3)

#ORF7(a & b)
ALPH_ORF7_AAs = Protein_Estimator(ALPH_ORF7_GENE)
DELT_ORF7_AAs = Protein_Estimator(DELT_ORF7_GENE)
OMIC_ORF7_AAs = Protein_Estimator(OMIC_ORF7_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_ORF7a_PROTEIN = ALPH_ORF7_AAs[4]
DELT_ORF7a_PROTEIN = DELT_ORF7_AAs[3]
OMIC_ORF7a_PROTEIN = OMIC_ORF7_AAs[2]
ALPH_ORF7b_PROTEIN = ALPH_ORF7_AAs[5]
DELT_ORF7b_PROTEIN = DELT_ORF7_AAs[4]
OMIC_ORF7b_PROTEIN = OMIC_ORF7_AAs[3]
#report results
PROT_DRIFT_REPORT(ALPH_ORF7a_PROTEIN,DELT_ORF7a_PROTEIN,['Alpha','Delta','ORF7a'],3)
PROT_DRIFT_REPORT(ALPH_ORF7a_PROTEIN,OMIC_ORF7a_PROTEIN,['Alpha','Omicron','ORF7a'],3)
PROT_DRIFT_REPORT(ALPH_ORF7b_PROTEIN,DELT_ORF7b_PROTEIN,['Alpha','Delta','ORF7b'],3)
PROT_DRIFT_REPORT(ALPH_ORF7b_PROTEIN,OMIC_ORF7b_PROTEIN,['Alpha','Omicron','ORF7b'],3)

#ORF10
ALPH_ORF10_AAs = Protein_Estimator(ALPH_ORF10_GENE)
DELT_ORF10_AAs = Protein_Estimator(DELT_ORF10_GENE)
OMIC_ORF10_AAs = Protein_Estimator(OMIC_ORF10_GENE)
#Isolate canon versions (manually checked NCBI)
ALPH_ORF10_PROTEIN = ALPH_ORF10_AAs[8]
DELT_ORF10_PROTEIN = DELT_ORF10_AAs[8]
OMIC_ORF10_PROTEIN = OMIC_ORF10_AAs[7]
#report results
PROT_DRIFT_REPORT(ALPH_ORF10_PROTEIN,DELT_ORF10_PROTEIN,['Alpha','Delta','ORF10'],3)
PROT_DRIFT_REPORT(ALPH_ORF10_PROTEIN,OMIC_ORF10_PROTEIN,['Alpha','Omicron','ORF10'],3)

#ORF1 - this can be SLOW
ALPH_ORF1_AAs = Protein_Estimator(ALPH_ORF1_GENE)
DELT_ORF1_AAs = Protein_Estimator(DELT_ORF1_GENE)
OMIC_ORF1_AAs = Protein_Estimator(OMIC_ORF1_GENE)
#Isolate canon versions (manually checked NCBI)
#Note - this is N terminus part, some splicing
#   or 'ribosomal slippage' occurs during translation
#   Documented Alpha is 'MESL...VNN'
#   while documented Delta & Omicron are 'MESL...FAV'
#   did a truncation occur?
ALPH_ORF1_PROTEIN = ALPH_ORF1_AAs[0]
DELT_ORF1_PROTEIN = DELT_ORF1_AAs[0]
OMIC_ORF1_PROTEIN = OMIC_ORF1_AAs[0]
#report results
PROT_DRIFT_REPORT(ALPH_ORF1_PROTEIN,DELT_ORF1_PROTEIN,['Alpha','Delta','ORF1 N terminus'],3)
PROT_DRIFT_REPORT(ALPH_ORF1_PROTEIN,OMIC_ORF1_PROTEIN,['Alpha','Omicron','ORF1 N terminus'],3)
#End of Other Protein Analysis
