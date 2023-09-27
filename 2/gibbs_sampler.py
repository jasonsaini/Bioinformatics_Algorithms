import random

# Generate a profile matrix from a set of motifs
def generate_profile(motifs):
    k = len(motifs[0])

    # build profile matrix w/pseudocoutns
    profile = []
    
    # create a 4 (ACGT) x k (motif length) matrix 
    for _ in range(4):
        row = []
        for _ in range(k):
            # add pseudo count
            row.append(1)
        profile.append(row)
        
    for i in range(k):
        column = [motif[i] for motif in motifs]
        
        # populate profile with nucleotides
        for j in range(4):
            profile[j][i] = (column.count("ACGT"[j]) + 1) / (len(column) + 4)
    
    return profile

# Find the most probable k-mer in a text given a profile matrix
def profile_most_probable(text, k, profile):
    max_prob = -1
    most_probable = ""
    
    for i in range(len(text) - k + 1):
        # read in k-mer
        kmer = text[i:i+k]
        prob = 1
        
        # assign nucleotide to correct position in profile
        for j in range(k):
            if kmer[j] == 'A':
                prob *= profile[0][j]
            elif kmer[j] == 'C':
                prob *= profile[1][j]
            elif kmer[j] == 'G':
                prob *= profile[2][j]
            elif kmer[j] == 'T':
                prob *= profile[3][j]
        
        if prob > max_prob:
            max_prob = prob
            most_probable = kmer
    
    return most_probable

# Calculate the consensus string score of a set of motifs
def calculate_score(motifs):
    consensus = ""
    k = len(motifs[0])
    
    for j in range(k):
        column = [motif[j] for motif in motifs]
        max_count = max([column.count("A"), column.count("C"), column.count("G"), column.count("T")])
        
        if column.count("A") == max_count:
            consensus += "A"
        elif column.count("C") == max_count:
            consensus += "C"
        elif column.count("G") == max_count:
            consensus += "G"
        elif column.count("T") == max_count:
            consensus += "T"
            
    score = 0
    
    for motif in motifs:
        for i in range(k):
           # If the consensus base at this position is different from the motif base, increment score counter
            if consensus[i] != motif[i]:
                score += 1
    return score

# Perform one iteration of Gibbs Sampler
def gibbs_sampler_iter(dna, k, t, motifs):
    i = random.randint(0, t - 1)
    motifs_except_i = motifs[:i] + motifs[i+1:]
    profile = generate_profile(motifs_except_i)
    most_probable = profile_most_probable(dna[i], k, profile)
    motifs[i] = most_probable
    
    return motifs

def gibbs_sampler(dna, k, t, N):
    best_motifs = [seq[:k] for seq in dna]
    
    for _ in range(N):
        best_motifs = gibbs_sampler_iter(dna, k, t, best_motifs)
    
    return best_motifs

# Sample Data & Output
k, t, N = 8, 5, 100
Dna = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

# Run GibbsSampler 20 times with pseudocounts
best_motifs = None
best_score = float('inf')

for _ in range(20):
    motifs = gibbs_sampler(Dna, k, t, N)
    score = calculate_score(motifs)
    if score < best_score:
        best_score = score
        best_motifs = motifs

# Print the resulting motifs
for motif in best_motifs:
    print(motif)
















print("ACGTCCACCGGCGTC AAGCGCACCGGGGTG ACCCTTACCGGGGTG AAGTTCCTCGGGGTG AAGTTTTATGGGGTG AAGTTTACCGGGTGC AAGTTTCGAGGGGTG CTGTTTACCGGGGTA AAGTTGCTCGGGGTG AAACATACCGGGGTG AAGTTTAGGAGGGTG AAGGAAACCGGGGTG AAGTTTACACAGGTG TAGTTTACCGGGGAT CCTTTTACCGGGGTG AAGTGAGCCGGGGTG AAGTCGTCCGGGGTG AAGTTTACCGGACAG AAGTTTACCAATGTG AAGTTTACCGTCATG")