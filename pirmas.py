from Bio import SeqIO
from Bio.SeqUtils import CodonUsage as CU
import Bio.Alphabet
import os
import numpy
import pandas as pd
from scipy.spatial import distance_matrix

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    answer.sort()
    return answer

def get_codon_frequency(orf_list):
    alphabet_list = list(Bio.Alphabet.IUPAC.IUPACProtein.letters)
    codon_alphabet = dict.fromkeys(alphabet_list, 0)

    total = 0
    for start, end, strand, pro in orf_list:
        total = total + len(pro)
        for codon in list(pro):
            codon_alphabet[codon] += 1
    for key in codon_alphabet:
        codon_alphabet[key] = codon_alphabet[key] / total
    return codon_alphabet

def get_dicodon_frequency(orf_list):
    alphabet_list = list(Bio.Alphabet.IUPAC.IUPACProtein.letters)
    dicodon_alphabet = dict.fromkeys([a + b for a in alphabet_list for b in alphabet_list], 0)

    total = 0
    for start, end, strand, pro in orf_list:
        if (len(pro) % 2 != 0):
            pro = pro[:-1]
        total = total + (len(pro) / 2)
        dicodon_list = [pro[i:i+2] for i in range(0, len(pro), 2)]
        for dicodon in dicodon_list:
            dicodon_alphabet[dicodon] += 1
    for key in dicodon_alphabet:
        dicodon_alphabet[key] = dicodon_alphabet[key] / total
    return dicodon_alphabet

def get_all_frequencies():   
    codon_frequencies = {}
    dicodon_frequencies = {}
    table = 11
    min_pro_len = 100
    for filename in os.listdir('data'):
        filepath = 'data/' + filename
        record = SeqIO.read(filepath, "fasta")
        orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
        codon_frequency = get_codon_frequency(orf_list)
        dicodon_frequency = get_dicodon_frequency(orf_list)
        codon_frequencies[record.name] = codon_frequency.values()
        dicodon_frequencies[record.name] = dicodon_frequency.values()
    return codon_frequencies, dicodon_frequencies

def print_distance_matrix(frequencies):
    names = []
    data = []
    for key in frequencies:
        names.append(key)
        data.append(list(frequencies[key]))
    df = pd.DataFrame(data, index=names)
    print(pd.DataFrame(distance_matrix(df.values, df.values), index=df.index))


codon_frequencies, dicodon_frequencies = get_all_frequencies()

print('\nCodon distances:')
print_distance_matrix(codon_frequencies)

print('\nDicodon distances:')
print_distance_matrix(dicodon_frequencies)