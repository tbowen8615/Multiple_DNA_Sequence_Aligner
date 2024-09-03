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

# Align Sequence to MSA
def align_sequence_with_msa(input_seq, msa, scoring_matrix):
    best_score = float('-inf')
    best_alignment = None
    best_msa = None

    # Perform sequence to MSA alignment
    for sequence in msa:
        print("      input_sequence: ", input_seq)
        print("   sequence_from_msa: ", sequence)
        aligned_input_seq, aligned_seq_from_msa, score, _ = needleman_wunsch(input_seq, sequence, scoring_matrix)
        print("   aligned_input_seq: ", aligned_input_seq)
        print("aligned_seq_from_msa: ", aligned_seq_from_msa,)
        print("               score: ", score)
        if score > best_score:
            best_score = score
            best_alignment = (aligned_input_seq, aligned_seq_from_msa)
            best_msa = [aligned_input_seq if sequence == input_seq else sequence for sequence in msa]
            print("\n                best_score: ", best_score)
            print("            best_alignment: ", best_alignment[0], "\n", "                           ", best_alignment[1], "\n")

    # Copy gaps from the best alignment to all sequences in the MSA
    print("\n\n                        STEP 4 SUBTASK D: Copy gaps from the best alignment to all sequences in the MSA\n")
    aligned_msa = []
    print("                        aligned_msa: ", aligned_msa, "\n")
    for sequence in msa:
        print("                        for sequence in msa: ")
        print("                                   sequence: ", sequence)
        print("                                        msa: ", msa, "\n")
        aligned_seq = ''
        seq_index = 0
        print("                                             aligned_seq: ", aligned_seq)
        print("                                               seq_index: ", seq_index, "\n")
        for char in best_alignment[1]:
            print("                                     char: ", char)
            print("                        best_alignment[1]: ", best_alignment[1], "\n")
            if char == '-':
                aligned_seq += '-'
                print("                                              aligned_seq: ", aligned_seq)
                print("                                                seq_index: ", seq_index, "\n")
            else:
                aligned_seq += sequence[seq_index]
                seq_index += 1
                print("                                              aligned_seq: ", aligned_seq)
                print("                                                seq_index: ", seq_index, "\n")
        aligned_msa.append(aligned_seq)
        print("                              aligned_msa: ", aligned_msa)
        print("                        best_alignment[0]: ", best_alignment[0], "\n")
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

# Progressive Alignment
def progressive_alignment(sequences, guide_tree, scoring_matrix):
    print("\n\n\nStep 4) Use the Input Sequences, Guide Tree, and Scoring Matrix as input to perform Multiple Sequence Alignment")
    # Print the input sequences for reference
    input_fp = sys.argv[sys.argv.index('-i') + 1]
    print("\nInput Sequences:")
    with open(input_fp, 'r') as file:
        for line in file:
            print(line, end='')

    print("\n\nGuide Tree:")
    for pair in guide_tree:
        print(f"{pair}")
    print()

    print_scoring_matrix(scoring_matrix)

    # Initialize clusters with individual sequences
    clusters = [[seq] for seq in sequences]
    cluster_map = {i: i for i in range(len(clusters))}

    for pair in guide_tree:
        print("\n\nSTEP 4 SUBTASK A: Clustering")
        print("clusters: " + str(clusters))
        print("cluster_map: " + str(cluster_map))
        print("Guide Tree pair: " + str(pair))
        m, n = pair

        # Find the indices in the current clusters
        i = cluster_map[m[0]]
        j = cluster_map[n[0]]
        print("i: " + str(i))
        print("j: " + str(j))

        cluster_i = clusters[i]
        cluster_j = clusters[j]
        print("cluster[" + str(i) + "]: " + str(cluster_i))
        print("cluster[" + str(j) + "]: " + str(cluster_j) + "\n")

        if len(cluster_i) == 1 and len(cluster_j) == 1:
            # Sequence to sequence alignment
            print("        STEP 4 SUBTASK B: Perform sequence to sequence alignment")
            print(f"        Alignment Type: sequence " + str(cluster_i) + " to sequence " + str(cluster_j) + "\n")
            aligned_seq1, aligned_seq2, _, _ = needleman_wunsch(cluster_i[0], cluster_j[0], scoring_matrix)
            print("        aligned_seq1: " + str(aligned_seq1))
            print("        aligned_seq2: " + str(aligned_seq2) + "\n")
            clusters[i] = [aligned_seq1]
            clusters[j] = [aligned_seq2]
            print("        cluster[" + str(i) + "] post alignment: " + str(clusters[i]))
            print("        cluster[" + str(j) + "] post alignment: " + str(clusters[j]))
            print("          clusters post alignment:\n        " + str(clusters))
        elif len(cluster_i) > 1 and len(cluster_j) == 1:
            # MSA to sequence alignment
            print("        STEP 4 SUBTASK B: Perform MSA to sequence alignment")
            print(f"        Alignment Type: MSA " + str(cluster_i) + " to sequence " + str(cluster_j) + "\n")
            aligned_seq, aligned_msa = align_sequence_with_msa(cluster_j[0], cluster_i, scoring_matrix)
            print("        aligned_seq: " + str(aligned_seq))
            print("        aligned_msa: " + str(aligned_msa) + "\n")
            clusters[j] = [aligned_seq]
            clusters[i] = aligned_msa
            print("        cluster[" + str(i) + "] post alignment: " + str(clusters[i]))
            print("        cluster[" + str(j) + "] post alignment: " + str(clusters[j]))
            print("          clusters post alignment:\n        " + str(clusters))
        elif len(cluster_i) == 1 and len(cluster_j) > 1:
            # Sequence to MSA alignment
            print("        STEP 4 SUBTASK B: Perform sequence to MSA alignment")
            print(f"        Alignment Type: sequence " + str(cluster_i) + " to MSA " + str(cluster_j) +"\n")
            aligned_seq, aligned_msa = align_sequence_with_msa(cluster_i[0], cluster_j, scoring_matrix)
            print("        aligned_seq: " + str(aligned_seq))
            print("        aligned_msa: " + str(aligned_msa) + "\n")
            clusters[i] = [aligned_seq]
            clusters[j] = aligned_msa
            print("        cluster[" + str(i) + "] post alignment: " + str(clusters[i]))
            print("        cluster[" + str(j) + "] post alignment: " + str(clusters[j]))
            print("          clusters post alignment:\n        " + str(clusters))
        else:
            # MSA to MSA alignment
            print("        STEP 4 SUBTASK B: Perform MSA to MSA alignment")
            print(f"        Alignment Type: MSA " + str(cluster_i) + " to MSA " + str(cluster_j) + "\n")
            aligned_msa1, aligned_msa2 = align_msa_with_msa(cluster_i, cluster_j, scoring_matrix)
            print("        aligned_msa1: " + str(aligned_msa1))
            print("        aligned_msa2: " + str(aligned_msa2) + "\n")
            clusters[i] = aligned_msa1
            clusters[j] = aligned_msa2
            print("        cluster[" + str(i) + "] post alignment: " + str(clusters[i]))
            print("        cluster[" + str(j) + "] post alignment: " + str(clusters[j]))
            print("          clusters post alignment:\n        " + str(clusters))

        # Merge clusters
        print("\n                STEP 4 SUBTASK C: Merge cluster[i] and cluster[j]")
        clusters[i].extend(clusters[j])
        clusters[j] = []
        print("                cluster[" + str(i) + "] after merge: " + str(clusters[i]))
        print("                cluster[" + str(j) + "] after merge: " + str(clusters[j]))
        print("                  clusters after merge:\n                " + str(clusters))

        # Update the cluster map to reflect the merged cluster
        for item in n:
            print("\n                Update cluster_map: " + str(cluster_map))
            print("                item: " + str(item))
            cluster_map[item] = i
            print("                cluster_map[item]: " + str(cluster_map[item]))
            print("                New cluster_map: " + str(cluster_map) + "\n")

    # Merge all clusters into a single MSA
    final_msa = [seq for cluster in clusters if cluster for seq in cluster]
    print("\nfinal_msa: " + str(final_msa))
    return final_msa

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
