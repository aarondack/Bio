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

# Rabbits and Recurrence Relations
# A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence (π,−2‾√,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.
#
# A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.
#
# When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.
#
# Given: Positive integers n≤40 and k≤5.
#
# Return: The total number of rabbit pairs that will be present after n months if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
def fibonacci(n,k):
        if n==0:
            return 1
        if n==1:
            return k
        elif n==2:
            return k
        else:
            return fibonacci(n-2,k)+fibonacci(n-1,k)


# Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), is the number of corresponding symbols that differ in s and t. See Figure 2.
#
# Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
#
# Return: The Hamming distance dH(s,t).
def hamming_distance(file):
    with open(file, 'r+') as file_object:
        for line in file_object:
            line2 = file_object.readline()
            data = [char_line for char_line,char_line2 in zip(line,line2) if char_line != char_line2]
            return(len(data))
#
# The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.
#
# The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
#
# Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
#
# Return: The protein string encoded by s.
def translating_rna_to_protein(rna):
    protein_string=""
    codon_dict = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L","UCU":"S","UCC":"S","UCA":"S", "UCG":"S",
    "UAU":"Y","UAC":"Y","UGU":"C","UGC":"C","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L","CCU":"P","CCC":"P","CCA":"P", "CCG":"P","CAU":"H",
    "CAC":"H","CAA":"Q","CAG":"Q","CGU":"R","CGC":"R","CGA":"R", "CGG":"R","AUU":"I","AUC":"I",
    "AUA":"I","AUG":"M","ACU":"T","ACC":"T","ACA":"T","ACG":"T","AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R","GUU":"V","GUC":"V","GUA":"V","GUG":"V","GCU":"A","GCC":"A",
    "GCA":"A","GCG":"A","GAU":"D","GAC":"D","GAA":"E","GAG":"E","GGU":"G","GGC":"G","GGA":"G","GGG":"G"
    }
    with open(rna, 'r') as rnaObject:
        for rna_string in rnaObject:
                for nucleotide in range(0,len(rna_string),3):
                   processing = rna_string[nucleotide:nucleotide+3]
                   if processing in codon_dict:
                            protein_string+=codon_dict[processing]
                   elif processing == "UAA" or "UAG" or "UGA":
                       break
                return protein_string