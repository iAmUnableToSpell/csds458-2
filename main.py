import random
import string
import numpy as np
from sequence_alignment import *

# alphabet to be used in the case no inputs are supplied
alphabet = ['A', 'C', 'G', 'T']

with open("./input.data") as file_in:
    lines = []
    for line in file_in:
        lines.append(line.replace('\n', ''))
isglobal = lines[0] == 'g'
score_match, score_mismatch, score_gap = lines[1].split(' ')

# if sample sequences are provided, write them to s1 and s2
try:
    s1 = lines[2]
    s2 = lines[3]
# if not, randomize from alphabet and write to s1 and s2
except Exception as e:
    s1 = ''.join(random.choices(alphabet, k=10))
    s2 = ''.join(random.choices(alphabet, k=(10 if isglobal else 20)))

# create sequence_alignment object
seq = sequence_alignment(s1, s2, int(score_gap), int(score_match), int(score_mismatch), not isglobal)

# find the most optimal alignment from the end of both strings (this populates the table)
seq.recursive_table_search(len(s1), len(s2))


# get coordinates for max values in score matrix
alignments = np.where(seq.T == seq.T.max())
# count number of ideal alignments (alignments with maximum score)
num_alignments = np.count_nonzero(seq.T == seq.T.max())

if isglobal:
    # get result strings from backtracking from the end of both s1 and s2
    a, b = seq.backtracking(len(s1), len(s2))
    # print results
    print(f"Globally Optimal Sequence ({int(seq.T[len(s1), len(s2)])}) at: [{len(s1)}, {len(s2)}]")
    print("".join(a))
    print("".join(b))
else:
    for i in range(len(alignments[0])):
        # print results for all alignments
        print(f"Optimal Sequence ({int(seq.T.max())}) at: [{alignments[0][i]}, {alignments[1][i]}]")
        a, b = seq.backtracking(alignments[0][i], alignments[1][i])
        print("\t" + "".join(a))
        print("\t" + "".join(b))
    print(f"{num_alignments} Optimal Alignment{'' if num_alignments == 1 else 's'}")








