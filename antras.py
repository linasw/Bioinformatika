from Bio import SeqIO


def find_lowest_and_highest_ascii():
    lowest_quality = 999
    highest_quality = -999
    with open("data/reads_for_analysis.fastq", "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            quality_list = record.letter_annotations['phred_quality']
            lowest = min(quality_list)
            highest = max(quality_list)
            if (lowest_quality > lowest):
                lowest_quality = lowest
            if (highest_quality < highest):
                highest_quality = highest
                if (highest_quality == 40):
                    print("KURWA")
    return chr(lowest_quality+33), chr(highest_quality+33)
            
# lowest, highest = find_lowest_and_highest_ascii()

import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8})

def generate_percentage_dictionary():
    my_dict = {}
    for i in range(0, 101, 1):
        my_dict[i] = 0
    return my_dict

def get_c_g_percentage():
    percent_dict = generate_percentage_dictionary()
    with open("data/reads_for_analysis.fastq", "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = record.seq
            c_g_sum = seq.count("C") + seq.count("G")
            total_sum = len(seq)
            percentage = int((c_g_sum / total_sum) * 100)
            percent_dict[percentage] += 1
    return percent_dict

def print_bar_chart(percentage_dict):
    objects = percentage_dict.keys()
    y_pos = np.arange(len(objects))
    performance = percentage_dict.values()

    plt.bar(y_pos, performance, align='center', alpha=0.5, width=0.3)
    plt.xticks(y_pos, objects)
    plt.xlabel('C/G Nukleotidu dalis sekoje')
    plt.ylabel('Readu skaicius')
    plt.title('C/G Nukleotidu dalis sekoje')

    plt.show()

# percent_dict = get_c_g_percentage()
# print_bar_chart(percent_dict)

from Bio.Blast import NCBIWWW

def get_peaks_seqs(peaks_list):
    sequences_dict = {}
    with open("data/reads_for_analysis.fastq", "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = record.seq
            c_g_sum = seq.count("C") + seq.count("G")
            total_sum = len(seq)
            percentage = int((c_g_sum / total_sum) * 100)
            if (percentage in peaks_list):
                sequences_dict[record.id] = seq
                peaks_list.remove(percentage)
                if not peaks_list:
                    return sequences_dict

def do_blast_search(sequences):
    with open('results.xml', 'w') as save_file:
        for value in sequences.values():
            result_handle = NCBIWWW.qblast("blastn", "nt", value)
            blast_results = result_handle.read() 
            save_file.write(blast_results)

# peaks = [32, 33, 34, 35, 36, 50, 51, 52, 53, 54, 68, 69, 70, 71, 72]
# peaks_sequences = get_peaks_seqs(peaks)
# do_blast_search(peaks_sequences)
