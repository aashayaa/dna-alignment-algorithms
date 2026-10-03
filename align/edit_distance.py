"""
edit_distance.py

"""

import numpy as np


def edit_distance(seq1: str, seq2: str) -> int:
    #Compute the edit distance between two sequences.
    n, m = len(seq1), len(seq2)
    table = np.zeros((n + 1, m + 1), dtype=int)

    # Base cases
    for i in range(n + 1):
        table[i][0] = i
    for j in range(m + 1):
        table[0][j] = j

    # Fill the rest of the table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                cost = 0
            else:
                cost = 1

            table[i][j] = min(
                table[i - 1][j] + 1,       # deletion
                table[i][j - 1] + 1,       # insertion
                table[i - 1][j - 1] + cost # match or substitution
            )

    return table[n][m]
