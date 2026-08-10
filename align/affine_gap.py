"""
affine_gap.py
=============
Global alignment with affine gap penalties (Gotoh's algorithm).

To track this, we need THREE matrices instead of one:
- M[i][j]  = best score of an alignment of seq1[:i], seq2[:j] that ENDS
             in a match/mismatch (no gap at the very end)
- Ix[i][j] = best score that ends with a gap in seq2 (i.e. seq1[i-1] is
             aligned against a gap)
- Iy[i][j] = best score that ends with a gap in seq1 (i.e. seq2[j-1] is
             aligned against a gap)

The overall best score at (i, j) is whichever of the three is highest.
"""

import numpy as np

NEG_INF = float("-inf")


def affine_gap_alignment(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_open: int = -3,
    gap_extend: int = -1,
) -> tuple[int, str, str]:
    """Compute optimal global alignment with affine gap penalties.

    Returns:
        (score, aligned_seq1, aligned_seq2)
    """
    n, m = len(seq1), len(seq2)

    M = np.full((n + 1, m + 1), NEG_INF)
    Ix = np.full((n + 1, m + 1), NEG_INF)  # gap in seq2 (consumes seq1 char)
    Iy = np.full((n + 1, m + 1), NEG_INF)  # gap in seq1 (consumes seq2 char)

    M[0][0] = 0

    # First row/column: only reachable by opening then extending a single gap
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        Iy[0][j] = gap_open + (j - 1) * gap_extend

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score

            # M: previous cell must have ended in ANY state (match, or
            # either kind of gap) -- a match/mismatch can follow anything.
            M[i][j] = max(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1]) + s

            # Ix: either OPEN a new gap from an M cell, or EXTEND an
            # existing gap from an Ix cell. (Biologically: you don't
            # usually switch straight from a seq1-gap to a seq2-gap
            # without a match between them, so we don't allow Iy -> Ix here.)
            Ix[i][j] = max(
                M[i - 1][j] + gap_open,
                Ix[i - 1][j] + gap_extend,
            )

            Iy[i][j] = max(
                M[i][j - 1] + gap_open,
                Iy[i][j - 1] + gap_extend,
            )

    # Best final score is whichever matrix is highest at the bottom-right corner
    final_scores = {"M": M[n][m], "Ix": Ix[n][m], "Iy": Iy[n][m]}
    state = max(final_scores, key=final_scores.get)
    score = final_scores[state]

    # Traceback: walk backward, following whichever matrix we're "in"
    aligned1, aligned2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        if state == "M":
            s = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            # figure out which matrix we came from
            prev = {"M": M[i - 1][j - 1], "Ix": Ix[i - 1][j - 1], "Iy": Iy[i - 1][j - 1]}
            state = max(prev, key=prev.get)
            i -= 1
            j -= 1
        elif state == "Ix":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            # came from opening (M) or extending (Ix) a gap
            if Ix[i][j] == M[i - 1][j] + gap_open:
                state = "M"
            else:
                state = "Ix"
            i -= 1
        else:  # state == "Iy"
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            if Iy[i][j] == M[i][j - 1] + gap_open:
                state = "M"
            else:
                state = "Iy"
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return int(score), "".join(aligned1), "".join(aligned2)