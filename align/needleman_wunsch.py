"""
needleman_wunsch.py
"""

import numpy as np


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
    return_table: bool = False,
):
    """Compute the optimal global alignment of two sequences.
    Returns:
        (score, aligned_seq1, aligned_seq2) 
    """
    n, m = len(seq1), len(seq2)
    table = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(n + 1):
        table[i][0] = i * gap_penalty
    for j in range(m + 1):
        table[0][j] = j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag = table[i - 1][j - 1] + match_score
            else:
                diag = table[i - 1][j - 1] + mismatch_score

            up = table[i - 1][j] + gap_penalty
            left = table[i][j - 1] + gap_penalty

            table[i][j] = max(diag, up, left)

    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        current = table[i][j]

        if i > 0 and j > 0:
            match = seq1[i - 1] == seq2[j - 1]
            diag_score = table[i - 1][j - 1] + (match_score if match else mismatch_score)
        else:
            diag_score = None

        if i > 0 and j > 0 and current == diag_score:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and current == table[i - 1][j] + gap_penalty:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    if return_table:
        return table[n][m], "".join(aligned1), "".join(aligned2), table
    return table[n][m], "".join(aligned1), "".join(aligned2)
