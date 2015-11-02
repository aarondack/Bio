__author__ = 'Aaron'

# Counting DNA Nucleotides
# A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.
#
# An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."
#
# Given: A DNA string s of length at most 1000 nt.
#
# Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.

def count_nucleotides(file):
    with open(file, 'r') as fileObject:
            a,c,g,t = 0,0,0,0
            for line in fileObject:
                for char in line:
                    if char == "A":
                        a+=1
                    elif char == "G":
                        g+=1
                    elif char == "T":
                        t+=1    
                    elif char == "C":
                        c+=1
    return a,c,g,t

# Transcribing DNA into RNA
# An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
#
# Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u.
#
# Given: A DNA string t having length at most 1000 nt.
#
# Return: The transcribed RNA string of t.

def transcribe_to_rna(dna,rna):
    with open(dna, 'r') as dnaObject:
        readIn = dnaObject.read()
        rnaData = readIn.replace("T","U")
        with open(rna, 'w') as rnaObject:
            rnaObject.write(rnaData)

# Complementing a Strand of DNA
# In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.
#
# The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").
#
# Given: A DNA string s of length at most 1000 bp.

# Return: The reverse complement sc of s.
def complement_dna(dna,cDNA):
        with open(dna,'r') as dnaObject:
            readIn = dnaObject.read()
            processed = []
            for line in readIn:
                  for nucleotide in line:
                      if nucleotide == "A":
                        processed.append("T")
                      elif nucleotide == "T":
                        processed.append("A")
                      elif nucleotide == "C":
                        processed.append("G")
                      elif nucleotide == "G":
                        processed.append("C")
            complement = ''.join(processed)[::-1]
            with open(cDNA, 'w') as cDNAObject:
                cDNAObject.write(complement)