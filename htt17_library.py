###   htt17-polyQ97 library, taking into account Hamming Distance when evaluating the mutations   ###
###   Marina Vallejo Vallés   ### 

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

DNA_Codons = {
    "A":["GCT","GCC","GCA","GCG"],                #alanine
    "C":["TGT","TGC"],                            #cysteine
    "D":["GAT","GAC"],                            #asparagine
    "E":["GAA","GAG"],                            #glutamate
    "F":["TTT","TTC"],                            #phenylalanine
    "G":["GGT","GGC","GGA","GGG"],                #glycine
    "H":["CAT","CAC"],                            #histidine
    "I":["ATA","ATT","ATC"],                      #isoleucine
    "K":["AAA","AAG"],                            #lysine
    "L":["TTA","TTG","CTT","CTC","CTA","CTG"],    #leucine
    "M":["ATG"],                                  #methionine
    "N":["AAT","AAC"],                            #asparagine
    "P":["CCT","CCC","CCA","CCG"],                #proline
    "Q":["CAA","CAG"],                            #glutamine
    "R":["CGT","CGC","CGA","CGG","AGA","AGG"],    #arginine
    "S":["TCT","TCC","TCA","TCG","AGT","AGC"],    #serine
    "T":["ACT","ACC","ACA","ACG"],                #threonine
    "V":["GTT","GTC","GTA","GTG"],                #valine
    "W":["TGG"],                                  #tryptophan
    "Y":["TAT","TAC"]                             #tyrosine
}


def HammingDistance(mutation, wt_codon): 
    ''' Calculate nucleotide changes between the evaluated mutation and the wt codon '''
    aa = 0
    i = 0
    while(aa < len(mutation)):
        if(mutation[aa] != wt_codon[aa]):
            i += 1
        aa += 1
    return i

def Obtain_WT_sequence():
    ''' Obtain htt WT sequence from FASTA file '''
    for record in SeqIO.parse("htt17.fa", "fasta"):
        return (str(record.seq))

split_strings = []
n  = 3 # 3 nt = 1 codon

for index in range(0, len(Obtain_WT_sequence()), n):   # WT sequence to chunks of codons
    split_strings.append(Obtain_WT_sequence()[index : index + n])

dna_seq = []
aa_seq = []
index = 0   

while index < (len(Obtain_WT_sequence())/3):  
    for key in DNA_Codons:
        codons =[]
        hamming_distance = {}
        change_wt = split_strings[index] 
        for x in DNA_Codons[key]:
            str1 = x
            codons.append(str1)
            hd = HammingDistance(str1, change_wt)
            hamming_distance[x]=hd
        print(hamming_distance)
        max_key = max(hamming_distance, key=hamming_distance.get)
        print(max_key)
        split_strings[index] = max_key #get the codon with the highest hamming distance, in case of a tie it selects the first one
        mut = ''.join(split_strings)
        dna_seq.append(mut) #contains all the mutated dna sequences 
        mut_to_prot = Seq(mut)
        mut_to_prot = mut_to_prot.translate() #translation of the sequence to protein
        aa_seq.append(mut_to_prot) #contains all the mutated protein sequences
        split_strings = []
        n  = 3
        for a in range(0, len(Obtain_WT_sequence()), n):   
            split_strings.append(Obtain_WT_sequence()[a : a + n])
        codons =[]
        hamming_distance = {}
    index = index+1


dna_seq.append(Obtain_WT_sequence())
dna_seq = pd.DataFrame(dna_seq)
dna_seq = dna_seq.drop_duplicates()


WT_protein = Seq(Obtain_WT_sequence())
WT_protein = WT_protein.translate()
aa_seq.append(WT_protein)
aa_seq = pd.DataFrame(aa_seq)
aa_seq = aa_seq.drop_duplicates()

dna_seq.to_csv('htt17_sequences_DNA.txt', header=None, index=None, sep=' ', mode='a')
aa_seq.to_csv('htt17_sequences_proteins.txt', header=None, index=None, sep=' ', mode='a')
