# dna-alignment-algorithms

These are python implementations of the classic DP algorithms for DNA sequence
alignment, built from scratch with NumPy. No existing alignment libraries were
used. I wanted to understand what tools like BWA or BLAST are
doing internally.

## What's in here

All in `align/`:

- **`needleman_wunsch.py`** - global alignment (aligns two full
  sequences end to end), single scoring matrix, linear gap penalty.
  Needleman-Wunsch and edit distance are almost identical.
  Edit distance is really just Needleman-Wunsch with a simpler scoring
  scheme (counting edits instead of comparing and maximizing a match score)
  and no traceback, since you only care about the final number, not the actual
  alignment.
- **`smith_waterman.py`** - local alignment (finds the best matching
  region between two sequences instead of forcing the whole thing to
  align). Scores can't go below 0, so the alignment can restart anywhere,
  and traceback starts from the highest scoring cell instead of the bottom
  right corner
- **`edit_distance.py`** - standard Levenshtein distance, minimum
  edits to turn one sequence into another
- **`affine_gap.py`** - global alignment but with affine gap penalties,
  meaning opening a gap costs more than extending one already open.
  More realistic biologically since a 5 base deletion is usually one
  event, not five. Uses three matrices (M, Ix, Iy) to track whether
  you're in a match or one of the two gap states. Affine gap is really
  a modification of Needleman-Wunsch, not local  alignment at all, it
  still aligns the full sequences end to end. Affine gap instead needs
  to know whether a gap was just opened or is already in progress

## How to use it

```python
from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman
from edit_distance import edit_distance
from affine_gap import affine_gap_alignment

score, a1, a2 = needleman_wunsch("GATTACA", "GCATGCU")
print(score, a1, a2)

score, a1, a2 = smith_waterman("ACACACTA", "AGCACACA")
print(score, a1, a2)

dist = edit_distance("kitten", "sitting")
print(dist)  # 3

score, a1, a2 = affine_gap_alignment("GATTACA", "GCATGCU", gap_open=-3, gap_extend=-1)
print(score, a1, a2)
```

Each one fills a DP table, then traces back through it to reconstruct
the actual alignment. Affine gap is the hardest part, needs 3 matrices
instead of 1 since extending a gap and opening one cost different
amounts, and traceback has to know which matrix each cell came from.
