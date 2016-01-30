# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE
@author: EMILY YEH
"""

import random
from amino_acids import aa, codons, aa_table
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

def get_complement(nucleotide):
    """ Returns the complementary nucleotide."""
    if nucleotide == "A":
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "G":
        return "C"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    def complement(nucleotide):
        list_nucleotides = list(dna)
        res = []
        for item in list_nucleotides:
            res.append(get_complement(item))
        return res
    
    complement_dna = complement(dna)

    def backwards(complement_dna):
        res = []
        index = len(complement_dna) - 1
        while index < len(complement_dna) and index >= 0:
            letter = complement_dna[index]
            res.append(letter)
            index = index - 1
        return res

    backwards_dna = backwards(complement_dna)

    dna_string = ''.join(backwards_dna)
    return dna_string


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
    codon and returns the sequence up to but not including the
    first in frame stop codon.  If there is no in frame stop codon,
    returns the whole string.
    dna: a DNA sequence
    returns: the open reading frame represented as a string"""
    codons = ["TAG", "TAA", "TGA"]
    def find_substring(dna,codons):
        for item in codons:
            index = dna.find(item)
            if index != -1:
                return index
        return -1

    if find_substring(dna,codons) == -1:
        return dna
    else:
        open_reading_frame = dna[0:find_substring(dna,codons)]
        return open_reading_frame


def find_all_ORFs_oneframe(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    res = []
    split = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    codons = ["TAG", "TAA", "TGA"]
    start_codon = "ATG"

    def find_substring_end(split,codons):
        res_end = []
        for item in codons:
            if item in split:
                res_end.append(split.index(item))
        return res_end

    def find_substring_start(split,start_codon):
        res_start = []
        for i in range(0,len(split)):
            if split[i] == start_codon:
                res_start.append(i)
        return res_start

    for i in find_substring_start(split,start_codon):
        for n in find_substring_end(split,codons):
            if i < n:
                new_list = split[i:(n)]
                new_string = ''.join(new_list)
                res.append(new_string)
            else:
                new_list = split[i:]
                new_string = ''.join(new_list)
                res.append(new_string)
    return res


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    ORFs_list = []
    split1 = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    split2 = [dna[i + 1:i + 4] for i in range(0, len(dna), 3)]
    split3 = [dna[i + 2:i + 5] for i in range(0, len(dna), 3)]
    codons = ["TAG", "TAA", "TGA"]
    start_codon = "ATG"

    def find_substring_end(split,codons):
        res_end = []
        test = []
        for item in codons:
            if item in split:
                test.append(split.index(item))
            else:
                test.append(-1)
        maximum = max(test)
        res_end.append(maximum)
        return res_end

    def find_substring_start(split,start_codon):
        res_start = []
        for i in range(0,len(split)):
            if split[i] == start_codon:
                res_start.append(i)
        return res_start

    def find_ORFs(split,start_codon,codons):
        for i in find_substring_start(split,start_codon):
            for n in find_substring_end(split,codons):
                if len(find_substring_start(split, start_codon)) > 1:
                    minimum = min(find_substring_start(split, start_codon))
                    if minimum < n and (n-minimum) > 0:
                        new_list = split[minimum:n]
                        new_string = ''.join(new_list)
                        ORFs_list.append(new_string)
                else:
                    if i < n:
                        new_list = split[i:n]
                        new_string = ''.join(new_list)
                        ORFs_list.append(new_string)
                    else:
                        new_list = split[i:]
                        new_string = ''.join(new_list)
                        ORFs_list.append(new_string)

    ORFs1 = find_ORFs(split1,start_codon,codons)
    ORFs2 = find_ORFs(split2,start_codon,codons)
    ORFs3 = find_ORFs(split3,start_codon,codons)
    ORFs_list.append(ORFs1)
    ORFs_list.append(ORFs2)
    ORFs_list.append(ORFs3)

    final_list = [x for x in ORFs_list if x is not None]

    res = list(set(final_list))
    return res


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']"""
    res = []
    strand1 = find_all_ORFs(dna)
    res.append(strand1)

    new_string = get_reverse_complement(dna)

    strand2 = find_all_ORFs(new_string)
    res = strand1 + strand2
    return res

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    res = []
    ORFs = find_all_ORFs_both_strands(dna)
    maximum_length = max(len(x) for x in ORFs)
    for item in ORFs:
        if len(item) == maximum_length:
            res.append(item)
    dna_string = ''.join(res)
    return dna_string

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

# if __name__ == "__main__":
#     import doctest
#     doctest.testmod()

# print rest_of_ORF("ATGCATGAATGTAGATAGATGTGCCC")
# print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
# print find_all_ORFs("ATGCATGAATGTAG")
# print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
print longest_ORF("ATGCGAATGTAGCATCAAA")