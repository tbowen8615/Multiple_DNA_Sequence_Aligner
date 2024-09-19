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
    print("Scoring Matrix:\n")
    print("    " + " ".join(f"{el:>4}" for el in elements))
    for el1 in elements:
        row = [f"{el1:>4}"]
        for el2 in elements:
            row.append(f"{scoring_matrix.get((el1, el2), 0):>4}")
        print(" ".join(row))
    print("")


# Implement the Needleman-Wunsch optimal global alignment algorithm
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


# Calculate the normalized distance
def calculate_normalized_distance(Sa, Smin, Smax):
    if Smax == Smin:
        return 0  # Avoid division by zero
    normalized_distance = (-1) * math.log((Sa - Smin) / (Smax - Smin))
    return normalized_distance


# Construct guide tree
def construct_guide_tree(normalized_distance_matrix):
    print("\nStep 3) Use the Normalized Distance Matrix as input to construct a Guide Tree\n")
    guide_tree = []

    # Keep track of the size of each cluster
    cluster_sizes = [1] * normalized_distance_matrix.shape[0]

    # Labels for each cluster (initially, each point is its own cluster)
    labels = [(i,) for i in range(normalized_distance_matrix.shape[0])]

    # Repeat the following while the distance matrix size is greater than 1
    while normalized_distance_matrix.shape[0] > 1:
        min_distance = np.inf
        closest_pair = None

        # Finding the closest pair
        for row in range(len(normalized_distance_matrix)):
            for column in range(row + 1, len(normalized_distance_matrix)):
                if normalized_distance_matrix[row, column] < min_distance:
                    min_distance = normalized_distance_matrix[row, column]
                    closest_pair = (row, column)

        row, column = closest_pair

        # Update guide tree with the merge
        guide_tree.append((labels[row], labels[column]))

        # Calculate the new merged row/column using the weighted average formula
        for i in range(normalized_distance_matrix.shape[0]):
            if i != row and i != column:
                new_distance = (normalized_distance_matrix[min(i, row), max(i, row)] * cluster_sizes[row] +
                                normalized_distance_matrix[min(i, column), max(i, column)] * cluster_sizes[column]) / (
                                    cluster_sizes[row] + cluster_sizes[column])
                normalized_distance_matrix[min(i, row), max(i, row)] = new_distance

        # Update cluster size and labels
        cluster_sizes[row] += cluster_sizes[column]
        labels[row] = labels[row] + labels[column]

        # Remove the old row and column corresponding to the second part of the pair
        normalized_distance_matrix = np.delete(normalized_distance_matrix, column, axis=0)
        normalized_distance_matrix = np.delete(normalized_distance_matrix, column, axis=1)
        cluster_sizes.pop(column)
        labels.pop(column)

        # Print statements for debugging
        print(f"Closest Pair: {closest_pair}")
        print("Updated Matrix:\n", normalized_distance_matrix)
        print("Updated Cluster Sizes:", cluster_sizes)
        print("Guide Tree:", guide_tree, "\n")

    return guide_tree

def sum_of_pairs(final_msa, scoring_matrix):
    score = 0
    # print("\nscore: " + str(score) + "\n")
    num_of_columns = len(final_msa[0])
    for column_index in range(num_of_columns):
        column = [seq[column_index] for seq in final_msa]
        for pair_element_1 in range(len(column)):
            for pair_element_2 in range(pair_element_1 + 1, len(column)):
                score += scoring_matrix.get((column[pair_element_1], column[pair_element_2]), 0)
                # print("pair_element_1: " + str(column[pair_element_1]))
                # print("pair_element_2: " + str(column[pair_element_2]))
                # print("Change in Score: " + str(scoring_matrix.get((column[pair_element_1], column[pair_element_2]), 0)))
                # print("\nscore: " + str(score) + "\n")
    return score

def progressive_alignment(guide_tree, sequences, scoring_matrix):
    alignments = {i: [sequences[i]] for i in range(len(sequences))}

    for cluster1, cluster2 in guide_tree:
        if len(cluster1) == 1 and len(cluster2) == 1:
            # Sequence to sequence alignment
            seq1 = alignments[cluster1[0]][0]
            seq2 = alignments[cluster2[0]][0]
            align1, align2, _, _ = needleman_wunsch(seq1, seq2, scoring_matrix)
            alignments[cluster1[0]] = [align1]
            alignments[cluster2[0]] = [align2]
        elif len(cluster1) == 1:
            # Sequence to MSA alignment
            seq = alignments[cluster1[0]][0]
            msa = alignments[cluster2[0]]
            align1, align2 = align_sequence_to_msa(seq, msa, scoring_matrix)
            alignments[cluster1[0]] = [align1]
            alignments[cluster2[0]] = align2
        elif len(cluster2) == 1:
            # Sequence to MSA alignment
            seq = alignments[cluster2[0]][0]
            msa = alignments[cluster1[0]]
            align1, align2 = align_sequence_to_msa(seq, msa, scoring_matrix)
            alignments[cluster2[0]] = [align1]
            alignments[cluster1[0]] = align2
        else:
            # MSA to MSA alignment
            msa1 = alignments[cluster1[0]]
            msa2 = alignments[cluster2[0]]
            align1, align2 = align_msa_to_msa(msa1, msa2, scoring_matrix)
            alignments[cluster1[0]] = align1
            alignments[cluster2[0]] = align2

    # The final alignment will be in the last remaining cluster
    final_alignment = []
    for key in sorted(alignments.keys()):
        final_alignment.extend(alignments[key])
    return final_alignment

def align_sequence_to_msa(seq, msa, scoring_matrix):
    best_score = float('-inf')
    best_alignment = None

    for msa_seq in msa:
        align1, align2, score, _ = needleman_wunsch(seq, msa_seq, scoring_matrix)
        if score > best_score:
            best_score = score
            best_alignment = (align1, align2)

    align1, align2 = best_alignment

    # Copy gaps from the best alignment to all sequences in the MSA
    aligned_msa = []
    for msa_seq in msa:
        new_seq = ''
        msa_index = 0
        for char in align2:
            if char == '-':
                new_seq += '-'
            else:
                new_seq += msa_seq[msa_index]
                msa_index += 1
        aligned_msa.append(new_seq)

    # Propagate gaps to the aligned sequence
    new_align1 = ''
    msa_index = 0
    for char in align2:
        if char == '-':
            new_align1 += '-'
        else:
            new_align1 += align1[msa_index]
            msa_index += 1

    return new_align1, aligned_msa

def align_msa_to_msa(msa1, msa2, scoring_matrix):
    best_score = float('-inf')
    best_alignment = None

    # Find the best alignment between any pair of sequences from the two MSAs
    for seq1 in msa1:
        for seq2 in msa2:
            align1, align2, score, _ = needleman_wunsch(seq1, seq2, scoring_matrix)
            if score > best_score:
                best_score = score
                best_alignment = (align1, align2)

    align1, align2 = best_alignment

    # Copy gaps from the best alignment to all sequences in the respective MSAs
    aligned_msa1 = []
    aligned_msa2 = []

    for msa_seq in msa1:
        new_seq = ''
        msa_index = 0
        for char in align1:
            if char == '-':
                new_seq += '-'
            else:
                new_seq += msa_seq[msa_index]
                msa_index += 1
        aligned_msa1.append(new_seq)

    for msa_seq in msa2:
        new_seq = ''
        msa_index = 0
        for char in align2:
            if char == '-':
                new_seq += '-'
            else:
                new_seq += msa_seq[msa_index]
                msa_index += 1
        aligned_msa2.append(new_seq)

    return aligned_msa1, aligned_msa2

def main():
    try:
        input_fp = sys.argv[sys.argv.index('-i') + 1]
        output_fp = sys.argv[sys.argv.index('-o') + 1]
        scoring_matrix_fp = sys.argv[sys.argv.index('-s') + 1]
        headers, sequences = read_fasta(input_fp)
        scoring_matrix = read_scoring_matrix(scoring_matrix_fp)

        # Print the input sequences
        print("\nInput Sequences:")
        with open(input_fp,'r') as file:
            for line in file:
                print(line, end='')
        print("\n")

        # Print the scoring matrix
        print_scoring_matrix(scoring_matrix)

        # Find the min and max values in the scoring matrix
        matrix_values = list(scoring_matrix.values())
        Smin = min(matrix_values)
        Smax = max(matrix_values)

        # Initialize matrices to store the alignment scores and normalized scores
        num_sequences = len(sequences)
        score_matrix = np.zeros((num_sequences, num_sequences), dtype=int)
        normalized_distance_matrix = np.zeros((num_sequences, num_sequences), dtype=float)
        aligned_pairs = []

        # Perform alignment for every pair of sequences
        for i, j in itertools.combinations(range(num_sequences), 2):
            align1, align2, score, alignment_length = needleman_wunsch(sequences[i], sequences[j], scoring_matrix)
            score_matrix[i, j] = score
            Smin_scaled = Smin * alignment_length
            Smax_scaled = Smax * alignment_length
            normalized_distance = calculate_normalized_distance(score, Smin_scaled, Smax_scaled)
            normalized_distance_matrix[i, j] = normalized_distance
            normalized_distance_matrix[j, i] = normalized_distance  # Populate both upper and lower triangle
            aligned_pairs.append((i, j, align1, align2, score, normalized_distance, alignment_length, Smin_scaled, Smax_scaled))

        # Print each aligned pair
        print("\n\nStep 1) Align all possible sequence pairs:")
        for (
                i, j, align1, align2, score, normalized_distance, alignment_length, Smin_scaled,
                Smax_scaled) in aligned_pairs:
            print(f"\nAlignment Pair: {headers[i]} and {headers[j]}:")
            print(f"Aligned Sequence {headers[i]}: {align1}")
            print(f"Aligned Sequence {headers[j]}: {align2}")
            print(f"Alignment Score: {score}")
            print(f"Alignment Length: {alignment_length}")
            print(f"Min Possible Alignment Score: {Smin_scaled}")
            print(f"Max Possible Alignment Score: {Smax_scaled}")
            print(f"Normalized Distance: {normalized_distance}")

        # Print the alignment score matrix
        print("\n\n\nStep 2) Construct a Normalized Distance Matrix:")
        print("\nInitial Distance Matrix (uses each Alignment Pair's Alignment Score):\n")
        for i in range(num_sequences):
            for j in range(num_sequences):
                if j > i:
                        print(f"{score_matrix[i, j]:5}", end=" ")
                else:
                    print("     ", end=" ")
            print()

        # Replace lower triangular portion including diagonal of normalized_distance_matrix with inf
        for i in range(num_sequences):
            for j in range(i + 1):
                normalized_distance_matrix[i, j] = np.inf

        # Print the normalized distance matrix
        print("Normalized Distance Matrix (uses each Alignment Pair's Normalized Distance):\n")
        for i in range(num_sequences):
            for j in range(num_sequences):
                if j > i:
                    print(f"{normalized_distance_matrix[i, j]:>10.4f}", end=" ")
                else:
                    print(f"{'':>10}", end=" ")
            print()

        # Print Raw Normalized Score Matrix
        print("Normalized Distance Matrix (sent to construct_guide_tree):\n")
        print(normalized_distance_matrix)
        print("\n")

        guide_tree = construct_guide_tree(normalized_distance_matrix)

        # Print the guide tree
        print("Guide Tree:")
        for pair in guide_tree:
            print(f"{pair}")

        # Perform progressive alignment
        final_alignment = progressive_alignment(guide_tree, sequences, scoring_matrix)

        # Calculate the final score
        final_score = sum_of_pairs(final_alignment, scoring_matrix)

        # Print the final alignment with scores
        print("\nFinal Alignment:")
        for header, sequence in zip(headers, final_alignment):
            print(f">{header}; score={final_score}")
            for i in range(0, len(sequence), 80):
                print(sequence[i:i + 80])

        # Write the final alignment to the output file
        write_fasta(output_fp, headers, final_alignment, [final_score] * len(final_alignment))

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()