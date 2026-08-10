"""
edit_distance.py
=================
Levenshtein edit distance: the minimum number of single-character
insertions, deletions, or substitutions needed to turn one string
into another. This is the simplest dynamic programming alignment
algorithm, and the DP table it builds is the same shape used by
Needleman-Wunsch and Smith-Waterman later.
"""

import numpy as np


def edit_distance(seq1: str, seq2: str) -> int:
    """Compute the edit distance between two sequences.

    Builds an (len(seq1)+1) x (len(seq2)+1) DP table where
    table[i][j] = edit distance between seq1[:i] and seq2[:j].
    """
    n, m = len(seq1), len(seq2)
    table = np.zeros((n + 1, m + 1), dtype=int)

    # Base cases: turning a string into "" costs len(string) deletions,
    # and turning "" into a string costs len(string) insertions.
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