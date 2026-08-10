"""
smith_waterman.py
==================
Smith-Waterman: local alignment of two sequences. Finds the
highest-scoring matching subsequence, rather than aligning the
sequences end-to-end like Needleman-Wunsch does.
"""

import numpy as np


def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
) -> tuple[int, str, str]:
    """Compute the optimal local alignment of two sequences.

    Returns:
        (score, aligned_seq1, aligned_seq2) -- the best local alignment
        score, and the aligned substrings (no leading/trailing gaps).
    """
    n, m = len(seq1), len(seq2)
    table = np.zeros((n + 1, m + 1), dtype=int)

    # Base cases: unlike Needleman-Wunsch, local alignment always allows
    # "starting fresh" at score 0, so the first row/column stay zero
    # instead of accumulating gap penalties.

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag = table[i - 1][j - 1] + match_score
            else:
                diag = table[i - 1][j - 1] + mismatch_score

            up = table[i - 1][j] + gap_penalty
            left = table[i][j - 1] + gap_penalty

            # The key difference from Needleman-Wunsch: clamp at 0.
            # A negative score means "start a new local alignment here"
            # rather than carrying forward a bad partial alignment.
            table[i][j] = max(0, diag, up, left)

            if table[i][j] > max_score:
                max_score = table[i][j]
                max_pos = (i, j)

    # Traceback: start from the highest-scoring cell (not the corner!),
    # and stop as soon as we hit a 0 -- that's where the local match begins.
    aligned1, aligned2 = [], []
    i, j = max_pos

    while i > 0 and j > 0 and table[i][j] != 0:
        current = table[i][j]

        match = seq1[i - 1] == seq2[j - 1]
        diag_score = table[i - 1][j - 1] + (match_score if match else mismatch_score)

        if current == diag_score:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current == table[i - 1][j] + gap_penalty:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return max_score, "".join(aligned1), "".join(aligned2)