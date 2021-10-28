#!/usr/bin/env python3
"""Exract portions of scaffolds from a fasta file and reverse comp if needed

Usage:
    <program> input_fasta input_info output_fasta

Where input_info has 5 tab-separated columns:
    - One random ID (Assemblytics_w_1)
    - Scaffold or sequence FULL name line
    - Start position
    - End position
    - Stand [-+]
"""

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return self.name + " " + self.sequence[:31]

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with myopen(input_file) as f:
        sequence = ""
        name = ""
        begun = False

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if begun:
                    yield Fasta(name, sequence)

                name = line[1:]
                sequence = ""
                begun = True

            else:
                sequence += line

        if name != "":
            yield Fasta(name, sequence)

def complement(seq):
    """Return the complement of a sequence, *NOT* it's reverse complement
    """

    if not seq.isalpha():
        print("The sequence contained non-alphabetic characters")
        print(seq)

    if not seq.isupper():
        seq = seq.upper()

    seq = seq.replace("A","1").replace("T","2").replace("C","3").replace("G","4")
    seq = seq.replace("a","5").replace("t","6").replace("c","7").replace("t","8")
    seq = seq.replace("1","T").replace("2","A").replace("3","G").replace("4","C")
    seq = seq.replace("5","t").replace("6","a").replace("7","g").replace("8","c")

    return seq

def reverse_complement(seq):

    return complement(seq)[::-1]

# Module
from collections import defaultdict
import gzip
import sys

# Parse user input
try:
    input_fasta = sys.argv[1]
    input_info = sys.argv[2]
    output_fasta = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Parse info file
wanted = defaultdict(list)

with myopen(input_info) as infile:
    for line in infile:
        id_, name, start, end, strand = line.strip().split()

        wanted[name].append(
                (int(start), int(end), id_, strand)
                )

# Read fasta file and extract needed regions
sequences = fasta_iterator(input_fasta)

with myopen(output_fasta, "wt") as outfile:
    for s in sequences:
        if s.name in wanted:

            for region in wanted[s.name]:
                start, end, id_, strand = region
                region_name = s.name + "_" + str(start) + "-" + str(end)
                region_seq = s.sequence[start: end + 1]

                if strand == "-":
                    region_seq = reverse_complement(region_seq)

                outfile.write("\t".join([id_, s.name, str(start), str(end), strand, region_seq]) + "\n")
