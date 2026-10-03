"""
run_examples.py
"""

import sys
import os


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from align.edit_distance import edit_distance
from align.needleman_wunsch import needleman_wunsch
from align.smith_waterman import smith_waterman
from align.affine_gap import affine_gap_alignment
from align.visualize import build_traceback_path, plot_alignment_matrix


def main():
    seq1, seq2 = "GATTACA", "GCATGCU"

    print(f"Comparing '{seq1}' and '{seq2}'\n")

    dist = edit_distance(seq1, seq2)
    print(f"Edit distance: {dist}")

    #Needleman-Wunsch
    score, a1, a2, table = needleman_wunsch(seq1, seq2, return_table=True)
    print(f"\nNeedleman-Wunsch score: {score}")
    print(a1)
    print(a2)

    #Smith-Waterman
    local_seq1, local_seq2 = "TGTTACGG", "GGTTGACTA"
    sw_score, sw_a1, sw_a2 = smith_waterman(local_seq1, local_seq2)
    print(f"\nSmith-Waterman score: {sw_score}")
    print(sw_a1)
    print(sw_a2)

    #Affine gap penalties 
    affine_seq1, affine_seq2 = "GATTACAGATTACA", "GATTAGATTACA"
    aff_score, aff_a1, aff_a2 = affine_gap_alignment(
        affine_seq1, affine_seq2, gap_open=-3, gap_extend=-1
    )
    print(f"\nAffine gap alignment score: {aff_score}")
    print(aff_a1)
    print(aff_a2)

    #traceback 
    os.makedirs("examples/output", exist_ok=True)
    path = build_traceback_path(
        seq1, seq2, table, match_score=1, mismatch_score=-1, gap_penalty=-2
    )
    plot_alignment_matrix(
        seq1, seq2, table, path=path,
        title=f"Needleman-Wunsch: {seq1} vs {seq2}",
        save_path="examples/output/nw_matrix.png",
    )
    print("\nSaved visualization to examples/output/nw_matrix.png")


if __name__ == "__main__":
    main()
