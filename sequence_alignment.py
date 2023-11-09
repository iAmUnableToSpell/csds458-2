import numpy as np


class sequence_alignment:
    def __init__(self, _s1: str, _s2: str, _gap_score: int, _match_score: int, _mismatch_score: int, _local: bool):
        self.s1 = _s1
        self.s2 = _s2
        self.gap_score = _gap_score
        self.match_score = _match_score
        self.mismatch_score = _mismatch_score
        self.local = _local
        self.T = self.init_matrix(len(self.s1), len(self.s2))
        self.num_paths = np.zeros(self.T.shape) - 1

    def init_matrix(self, n: int, m: int):
        # returns an n x m matrix with [:,0] and [0,:] initialized to reflect score_gap
        T = np.zeros((n + 1, m + 1))
        T[:, :] = np.NaN
        for i in range(n + 1):
            T[i, 0] = i * self.gap_score if not self.local else 0
        for i in range(m + 1):
            T[0, i] = i * self.gap_score if not self.local else 0
        return T

    def recursive_table_search(self, i: int, j: int):
        # recursively fill score table starting from i, j
        if not np.isnan(self.T[i, j]):
            # value has already been memoized, so return it
            return self.T[i, j]
        else:
            # value has not been memoized, so memoize it and write to T
            replace = self.recursive_table_search(i - 1, j - 1) + self.does_match(i, j)
            s1_miss = self.recursive_table_search(i - 1, j) + self.gap_score
            s2_miss = self.recursive_table_search(i, j - 1) + self.gap_score
            limit = 0 if self.local else float('-inf')
            self.T[i, j] = max(replace, s1_miss, s2_miss, limit)
            return self.T[i, j]

    def does_match(self, i: int, j: int):
        # helper function to return score of some (i, j) match attempt
        if self.s1[i - 1] == self.s2[j - 1]:
            return self.match_score
        else:
            return self.mismatch_score

    def get_all_alginments(self, i, j):
        alignments = []
        options = []
        if self.T[i - 1, j] + self.gap_score == self.T[i, j]:
            options.append(self.T[i - 1, j])
        if self.T[i - 1, j - 1] + self.match_score == self.T[i, j] or self.T[i - 1, j - 1] + self.mismatch_score == self.T[i, j]:
            options.append(self.T[i - 1, j - 1])
        if self.T[i, j - 1] + self.gap_score == self.T[i, j]:
            options.append(self.T[i, j - 1])
        if not options:
            return [([''], [''])]
        maxval = max(options)
        if j != 0 and self.T[i, j - 1] == maxval and self.T[i, j - 1] + self.gap_score == self.T[i, j]:
            # a gap in s1
            valid_seq = (self.get_all_alginments(i, j - 1))
            for path in valid_seq:
                path[0].append('-')
                path[1].append(self.s2[j - 1])
            alignments.extend(valid_seq)
        if i != 0 and self.T[i - 1, j] == maxval and self.T[i - 1, j] + self.gap_score == self.T[i, j]:
            # a gap in s1
            valid_seq = (self.get_all_alginments(i - 1, j))
            for path in valid_seq:
                path[0].append(self.s1[i - 1])
                path[1].append('-')
            alignments.extend(valid_seq)
        if (i != 0 and j != 0) and (self.T[i - 1, j - 1] == maxval) and (self.T[i - 1, j - 1] + self.match_score == self.T[i, j] or self.T[i - 1, j - 1] + self.mismatch_score == self.T[i, j]):
            # a match or missmatch valid
            valid_seq = self.get_all_alginments(i - 1, j - 1)
            for path in valid_seq:
                path[0].append(self.s1[i - 1])
                path[1].append(self.s2[j - 1])
            alignments.extend(valid_seq)
        return alignments
