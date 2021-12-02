# BioinformaticsDemo
Independent study from a uDemy bioinformatics course, this script uses Biopython functions to analyze SARSCOV2 genomes

With Biopython dependences some 'general use' functions are defined here: Protein_Estimator() and CODON_READER() - which can create potential amino acids sequences from any given nucleotide sequence.

Running this script creates potential protein candidates from the spike gene segment of the SARSCOV2 genomes as well as the Pfizer mRNA sequence and displays results for DELTA_LIST[10] which corresponds to NCBI YP_009724390.1, MU_LIST[8] which corresponds to NCBI QZW69948.1, and OM_LIST[8] which corresponds to NCBI UFO69279.1.
