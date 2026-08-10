"""
visualize.py
============
Renders a DP alignment matrix as a heatmap, with the traceback path
drawn on top -- turns the abstract table of numbers into a picture of
the alignment actually being built.
"""

import numpy as np
import matplotlib.pyplot as plt


def build_traceback_path(seq1: str, seq2: str, table: np.ndarray,
                          match_score: int, mismatch_score: int, gap_penalty: int,
                          start: tuple[int, int] = None, stop_at_zero: bool = False) -> list[tuple[int, int]]:
    """Re-walk a Needleman-Wunsch or Smith-Waterman traceback, but just
    return the list of (i, j) cell coordinates visited, instead of
    rebuilding aligned strings. Used only for plotting.
    """
    n, m = len(seq1), len(seq2)
    i, j = start if start else (n, m)
    path = [(i, j)]

    while i > 0 or j > 0:
        if stop_at_zero and table[i][j] == 0:
            break

        if i > 0 and j > 0:
            match = seq1[i - 1] == seq2[j - 1]
            diag_score = table[i - 1][j - 1] + (match_score if match else mismatch_score)
        else:
            diag_score = None

        if i > 0 and j > 0 and table[i][j] == diag_score:
            i, j = i - 1, j - 1
        elif i > 0 and table[i][j] == table[i - 1][j] + gap_penalty:
            i -= 1
        elif j > 0:
            j -= 1
        else:
            break

        path.append((i, j))

    return path


def plot_alignment_matrix(
    seq1: str,
    seq2: str,
    table: np.ndarray,
    path: list[tuple[int, int]] = None,
    title: str = "Alignment DP Matrix",
    save_path: str = None,
):
    """Plot a DP matrix as a heatmap with sequence letters as axis labels,
    and (optionally) a traceback path drawn as a line over the cells it visits.

    Args:
        seq1, seq2: the two sequences (used for axis labels).
        table: the (len(seq1)+1) x (len(seq2)+1) DP matrix.
        path: list of (i, j) coordinates from build_traceback_path,
            drawn as a connected line if given.
        title: plot title.
        save_path: if given, saves the figure to this path instead of
            just displaying it (useful for embedding in a README).
    """
    fig, ax = plt.subplots(figsize=(len(seq2) * 0.6 + 2, len(seq1) * 0.6 + 2))

    im = ax.imshow(table, cmap="viridis")
    fig.colorbar(im, ax=ax, label="score")

    ax.set_xticks(range(len(seq2) + 1))
    ax.set_xticklabels(["-"] + list(seq2))
    ax.set_yticks(range(len(seq1) + 1))
    ax.set_yticklabels(["-"] + list(seq1))

    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            ax.text(j, i, int(table[i][j]), ha="center", va="center",
                     color="white", fontsize=8)

    if path:
        ys = [p[0] for p in path]
        xs = [p[1] for p in path]
        ax.plot(xs, ys, color="red", linewidth=2, marker="o", markersize=4)

    ax.set_title(title)
    ax.set_xlabel("seq2")
    ax.set_ylabel("seq1")
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
    else:
        plt.show()