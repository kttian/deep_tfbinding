#https://www.biostars.org/p/710/
from itertools import groupby
def fasta_iter(fasta_name):
    """
        given a fasta file, yield tuples of (header, sequence)
    """
    fh = open(fasta_name) # file handle
    # ditch the boolean (x[0]) and just keep the header or sequence since they alternate
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        header = header.next()[1:].strip() # drop the ">" from the header
        seq = "".join(s.strip() for s in fa_iter.next()) # join all sequence lines to one
        yield header, seq

def gen_allele(input_name):
    orig_seqs = []
    alleles  = [[], [], []]
    fh = {}
    for i in range(3):
        fh[i] = open(input_name + str(i+1), "w")

    fasta = fasta_iter(input_name)

    letters = {'A', 'C', 'G', 'T'}
    off = int((1000-1)/2)
    for header, seq in fasta:   
        u_seq = list(seq.upper())
        orig_seqs.append(u_seq)
        ref = u_seq[off]
        if ref not in letters:
            allele_list = ['N', 'N', 'N']
        else:
            allele_list = letters.difference({ref})

        for idx, allele in enumerate(allele_list):
            u_seq[off] = allele
            fh[idx].write(">" + header + "\n")
            fh[idx].write("".join(u_seq) + "\n")

gen_allele("interpret.fa")


