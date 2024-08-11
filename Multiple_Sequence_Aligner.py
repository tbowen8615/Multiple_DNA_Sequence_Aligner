"""
=======================================================================================================================
Title           : Multiple DNA Sequence Aligner
Description     : This program performs Multiple Sequence Alignment of multiple DNA sequences. First, all possible
                    pairs of input sequences are aligned and a distance matrix is constructed. Then, a normalized
                    distance matrix is constructed. Using the normalized distance matrix as input, a guide tree is
                    constructed according to the UPGMA method. The Progressive Alignment function then takes this guide
                    tree as input and performs either sequence to sequence, sequence to MSA, or MSA to MSA alignments
                    until the final alignment is complete. A score is calculated and assigned to the final alignment
                    using the sum of pairs method. The final alignment is then output to the output file specified by
                    the user.
Author          : Tyler Bowen
Date            : July 20, 2024
Version         : 1
Usage           : In the Command Prompt, navigate to the directory containing this python script, an input FASTA
                    (.fasta) file, and a scoring matrix, then enter following command:
                        python Multiple_Sequence_Aligner.py -i MSA_input.fasta -o MSA_complete.fasta
                        -s BLOSUM50.mtx
                    You will specify an output FASTA file that the final alignment will be written to in the same
                    directory. Example input and scoring matrix files are included in the repository.
                        Input file: MSA_input.fasta
                        Scoring Matrix: BLOSUM50.mtx or BLOSUM62.mtx or nucleotide.mtx
                        Output: You will specify an output file name e.g. MSA_complete.fasta
Notes           :
Python Version  : 3.11.7
=======================================================================================================================
"""
import sys
import numpy as np
import itertools
import math

# Read sequences from a FASTA file
def read_fasta(file_path):
    with open(file_path, 'r') as f:
        headers, sequences = [], []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                headers.append(line[1:])
                sequences.append("")
            else:
                sequences[-1] += line
    return headers, sequences


# Write sequences to a FASTA file
def write_fasta(file_path, headers, sequences, scores):
    with open(file_path, 'w') as f:
        for header, sequence, score in zip(headers, sequences, scores):
            f.write(f">{header}; score={score}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + '\n')


# Read a scoring matrix from a file
def read_scoring_matrix(file_path):
    with open(file_path, 'r') as f:
        lines = [line.strip().split() for line in f if line.strip()]
    matrix = {}
    for line in lines[1:]:
        for col, score in zip(lines[0], line[1:]):
            matrix[(line[0], col)] = int(score)  # Populate the scoring matrix dictionary
    return matrix


def print_scoring_matrix(scoring_matrix):
    # Extract unique elements
    elements = sorted(set([key[0] for key in scoring_matrix.keys()] + [key[1] for key in scoring_matrix.keys()]))
    # Print the header row
    print("Scoring Matrix:")
    print("    " + " ".join(f"{el:>4}" for el in elements))
    for el1 in elements:
        row = [f"{el1:>4}"]
        for el2 in elements:
            row.append(f"{scoring_matrix.get((el1, el2), 0):>4}")
        print(" ".join(row))
    print("\n")


# Implement the Needleman-Wunsch algorithm
def needleman_wunsch(seq1, seq2, scoring_matrix):
    m, n = len(seq1), len(seq2)

    # Initialize scoring and traceback matrices
    F = np.zeros((m + 1, n + 1), dtype=int)
    traceback = np.zeros((m + 1, n + 1), dtype=int)

    # Initialize the first row and column
    for i in range(1, m + 1):
        F[i, 0] = F[i - 1, 0] + scoring_matrix[(seq1[i - 1], '-')]  # Gap penalty for seq1
        traceback[i, 0] = 1

    for j in range(1, n + 1):
        F[0, j] = F[0, j - 1] + scoring_matrix[('-', seq2[j - 1])]  # Gap penalty for seq2
        traceback[0, j] = 2

    # Fill in the scoring matrix and traceback matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = F[i - 1, j - 1] + scoring_matrix[(seq1[i - 1], seq2[j - 1])]  # Match/Mismatch score
            delete = F[i - 1, j] + scoring_matrix[(seq1[i - 1], '-')]  # Deletion score
            insert = F[i, j - 1] + scoring_matrix[('-', seq2[j - 1])]  # Insertion score
            F[i, j], traceback[i, j] = max((match, 0), (delete, 1), (insert, 2))  # Max score and traceback

    # Trace back to find the optimal alignment
    align1 = ''
    align2 = ''
    i = m
    j = n

    while i > 0 or j > 0:
        if traceback[i, j] == 0:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif traceback[i, j] == 1:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    alignment_length = len(align1)
    return align1, align2, F[m, n], alignment_length


# Calculate the normalized score
def calculate_normalized_score(Sa, Smin, Smax):
    if Smax == Smin:
        return 0  # Avoid division by zero
    normalized_score = (-1) * math.log((Sa - Smin) / (Smax - Smin))
    return normalized_score


# Construct guide tree
def construct_guide_tree(distance_matrix):
    # Replace None with np.inf to represent infinite distances
    distance_matrix = np.array([[np.inf if x is None else x for x in row] for row in distance_matrix], dtype=float)

    cluster_sizes = [1] * len(distance_matrix)
    labels = [(i,) for i in range(len(distance_matrix))]  # Initialize labels as tuples containing single integers

    guide_tree = []

    while distance_matrix.shape[0] > 1:
        # Initialize variables to store the minimum distance and corresponding indices
        min_distance = np.inf
        closest_pair = None

        # Iterate through the upper triangular part of the matrix to find the minimum distance
        for i in range(len(distance_matrix)):
            for j in range(i + 1, len(distance_matrix)):
                if distance_matrix[i, j] < min_distance:
                    min_distance = distance_matrix[i, j]
                    closest_pair = (i, j)

        i, j = closest_pair

        # Record the merging order using integer pairs
        guide_tree.append((labels[i][0], labels[j][0]))

        # Calculate the new distances
        for k in range(len(distance_matrix)):
            if k != i and k != j:
                distance_matrix[min(i, k), max(i, k)] = (distance_matrix[min(i, k), max(i, k)] * cluster_sizes[i] +
                                                         distance_matrix[min(j, k), max(j, k)] * cluster_sizes[j]) / (
                                                                cluster_sizes[i] + cluster_sizes[j])

        # Merge clusters
        labels[i] = (labels[i][0], labels[j][0])
        cluster_sizes[i] += cluster_sizes[j]

        # Remove the merged node's row and column (node j)
        distance_matrix = np.delete(distance_matrix, j, axis=0)
        distance_matrix = np.delete(distance_matrix, j, axis=1)
        del labels[j]
        del cluster_sizes[j]

    return guide_tree


def align_sequence_with_msa(seq, msa, scoring_matrix):
    best_score = float('-inf')
    best_alignment = None
    best_msa = None

    for y in msa:
        align1, align2, score, _ = needleman_wunsch(seq, y, scoring_matrix)
        if score > best_score:
            best_score = score
            best_alignment = (align1, align2)
            best_msa = [align1 if y == seq else y for y in msa]

    # Copy gaps from the best alignment to all sequences in the MSA
    aligned_msa = []
    for y in msa:
        aligned_seq = ''
        seq_idx = 0
        for char in best_alignment[1]:
            if char == '-':
                aligned_seq += '-'
            else:
                aligned_seq += y[seq_idx]
                seq_idx += 1
        aligned_msa.append(aligned_seq)

    return best_alignment[0], aligned_msa

def align_msa_with_msa(msa1, msa2, scoring_matrix):
    best_score = float('-inf')
    best_pair = None
    best_msa1 = None
    best_msa2 = None

    for x in msa1:
        for y in msa2:
            align1, align2, score, _ = needleman_wunsch(x, y, scoring_matrix)
            if score > best_score:
                best_score = score
                best_pair = (align1, align2)
                best_msa1 = [align1 if x == seq else seq for seq in msa1]
                best_msa2 = [align2 if y == seq else seq for seq in msa2]

    # Copy gaps from best alignment to all sequences in both MSAs
    aligned_msa1 = []
    for x in msa1:
        aligned_seq = ''
        seq_idx = 0
        for char in best_pair[0]:
            if char == '-':
                aligned_seq += '-'
            else:
                aligned_seq += x[seq_idx]
                seq_idx += 1
        aligned_msa1.append(aligned_seq)

    aligned_msa2 = []
    for y in msa2:
        aligned_seq = ''
        seq_idx = 0
        for char in best_pair[1]:
            if char == '-':
                aligned_seq += '-'
            else:
                aligned_seq += y[seq_idx]
                seq_idx += 1
        aligned_msa2.append(aligned_seq)

    return aligned_msa1, aligned_msa2

def replace_xs_with_gaps(msa):
    return [seq.replace('X', '-') for seq in msa]

def progressive_alignment(sequences, guide_tree, scoring_matrix):
    # Initialize clusters with individual sequences
    clusters = [[seq] for seq in sequences]

    for (i, j) in guide_tree:
        # Identify clusters to be aligned
        cluster_i = clusters[i]
        cluster_j = clusters[j]

        if len(cluster_i) == 1 and len(cluster_j) == 1:
            # Sequence to sequence alignment
            print(f"Alignment Type: sequence " + str(cluster_i) + " to sequence " + str(cluster_j) + "\n")
            aligned_seq1, aligned_seq2, _, _ = needleman_wunsch(cluster_i[0], cluster_j[0], scoring_matrix)
            clusters[i] = [aligned_seq1]
            clusters[j] = [aligned_seq2]
        elif len(cluster_i) > 1 and len(cluster_j) == 1:
            # MSA to sequence alignment
            print(f"Alignment Type: MSA " + str(cluster_i) + " to sequence " + str(cluster_j) + "\n")
            aligned_seq, aligned_msa = align_sequence_with_msa(cluster_j[0], cluster_i, scoring_matrix)
            clusters[j] = [aligned_seq]
            clusters[i] = aligned_msa
        elif len(cluster_i) == 1 and len(cluster_j) > 1:
            # Sequence to MSA alignment
            print(f"Alignment Type: sequence " + str(cluster_i) + " to MSA " + str(cluster_j) + "\n")
            aligned_seq, aligned_msa = align_sequence_with_msa(cluster_i[0], cluster_j, scoring_matrix)
            clusters[i] = [aligned_seq]
            clusters[j] = aligned_msa
        else:
            # MSA to MSA alignment
            print(f"Alignment Type: MSA " + str(cluster_i) + " to MSA " + str(cluster_j) + "\n")
            aligned_msa1, aligned_msa2 = align_msa_with_msa(cluster_i, cluster_j, scoring_matrix)
            clusters[i] = aligned_msa1
            clusters[j] = aligned_msa2

        # Merge clusters
        clusters[i].extend(clusters[j])
        clusters[j] = []

    # Merge all clusters into a single MSA
    final_msa = [seq for cluster in clusters if cluster for seq in cluster]
    final_msa = replace_xs_with_gaps(final_msa)
    return final_msa

def sum_of_pairs(final_msa, scoring_matrix):
    score = 0
    m = len(final_msa[0])
    for i in range(m):
        column = [seq[i] for seq in final_msa]
        for k in range(len(column)):
            for l in range(k + 1, len(column)):
                score += scoring_matrix.get((column[k], column[l]), 0)
    return score


def main():
    print("Assignment 4 :: R11474743")
    try:
        input_fp = sys.argv[sys.argv.index('-i') + 1]
        output_fp = sys.argv[sys.argv.index('-o') + 1]
        scoring_matrix_fp = sys.argv[sys.argv.index('-s') + 1]
        headers, sequences = read_fasta(input_fp)
        scoring_matrix = read_scoring_matrix(scoring_matrix_fp)

        # Print the scoring matrix
        print_scoring_matrix(scoring_matrix)

        # Find the min and max values in the scoring matrix
        matrix_values = list(scoring_matrix.values())
        Smin = min(matrix_values)
        Smax = max(matrix_values)

        # Initialize matrices to store the alignment scores and normalized scores
        num_sequences = len(sequences)
        score_matrix = np.zeros((num_sequences, num_sequences), dtype=int)
        normalized_score_matrix = np.zeros((num_sequences, num_sequences), dtype=float)
        aligned_pairs = []

        # Perform alignment for every pair of sequences
        for i, j in itertools.combinations(range(num_sequences), 2):
            align1, align2, score, alignment_length = needleman_wunsch(sequences[i], sequences[j], scoring_matrix)
            score_matrix[i, j] = score
            Smin_scaled = Smin * alignment_length
            Smax_scaled = Smax * alignment_length
            normalized_score = calculate_normalized_score(score, Smin_scaled, Smax_scaled)
            normalized_score_matrix[i, j] = normalized_score
            normalized_score_matrix[j, i] = normalized_score  # Populate both upper and lower triangle
            aligned_pairs.append((i, j, align1, align2, score, normalized_score, alignment_length, Smin_scaled, Smax_scaled))

        # Construct guide tree
        guide_tree = construct_guide_tree(normalized_score_matrix)

        # Print the alignment score matrix
        print("Initial Distance Matrix:")
        for i in range(num_sequences):
            for j in range(num_sequences):
                if j > i:
                    print(f"{score_matrix[i, j]:5}", end=" ")
                else:
                    print("     ", end=" ")
            print()

        # Print the normalized score matrix
        print("\nNormalized Distance Matrix:")
        for i in range(num_sequences):
            for j in range(num_sequences):
                if j > i:
                    print(f"{normalized_score_matrix[i, j]:>10.4f}", end=" ")
                else:
                    print(f"{'':>10}", end=" ")
            print()

        # Print each aligned pair
        print("\nAligned sequence pairs:")
        for (
                i, j, align1, align2, score, normalized_score, alignment_length, Smin_scaled,
                Smax_scaled) in aligned_pairs:
            print(f"\nAlignment between {headers[i]} and {headers[j]}:")
            print(f"{headers[i]}: {align1}")
            print(f"{headers[j]}: {align2}")
            print(f"Score: {score}")
            print(f"Alignment Length: {alignment_length}")
            print(f"Smin: {Smin_scaled}")
            print(f"Smax: {Smax_scaled}")
            print(f"Normalized Score: {normalized_score}")

        # Print the guide tree
        print("\nGuide Tree:")
        for pair in guide_tree:
            print(f"{pair}")
        print("\n")

        final_msa = progressive_alignment(sequences, guide_tree, scoring_matrix)

        # Calculate sum of pairs and assign the same score to each sequence in the final_msa list
        scores = [sum_of_pairs(final_msa, scoring_matrix)] * len(final_msa)

        # Write the final MSA to the output file
        write_fasta(output_fp, headers, final_msa, scores)

        print("\nFinal Alignment:")
        for header, sequence, score in zip(headers, final_msa, scores):
            print(f">{header}; score={score}")
            for i in range(0, len(sequence), 80):
                print(sequence[i:i + 80])

        print("\nAlignments performed successfully. Output written to", output_fp)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
