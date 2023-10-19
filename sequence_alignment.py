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

    def init_matrix(self, n: int, m: int):
        # returns an n x m matrix with [:,0] and [0,:] initialized to reflect score_gap
        T = np.zeros((n + 1, m + 1))
        T[:,:] = np.NaN
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
            replace = self.recursive_table_search(i-1, j-1) + self.does_match(i, j)
            s1_miss = self.recursive_table_search(i-1, j) + self.gap_score
            s2_miss = self.recursive_table_search(i, j-1) + self.gap_score
            limit = 0 if self.local else float('-inf')
            self.T[i, j] = max(replace, s1_miss, s2_miss, limit)
            return self.T[i, j]

    def does_match(self, i: int, j: int):
        # helper function to return score of some (i, j) match attempt
        if self.s1[i-1] == self.s2[j-1]:
            return self.match_score
        else:
            return self.mismatch_score

    def backtracking(self,i: int, j: int):
        # backtracks score matrix to find the ideal path, building two char arrays along the way
        xseq = []
        yseq = []
        if self.local:
            while i > 0 or j > 0:
                options = [self.T[i - 1, j], self.T[i, j - 1], self.T[i - 1, j - 1]]
                if (i < 1):
                    options.remove(self.T[i - 1, j])
                    options.remove(self.T[i - 1, j - 1])
                if (j < 1):
                    options.remove(self.T[i, j - 1])
                    options.remove(self.T[i - 1, j - 1])
                choice = max(options)
                # choice = np.argmax(options)
                if choice == self.T[i - 1, j - 1]:
                    # diag
                    xseq.append(self.s1[i - 1])
                    yseq.append(self.s2[j - 1])
                    i -= 1
                    j -= 1
                elif choice == self.T[i, j - 1]:
                    # left
                    xseq.append('-')
                    yseq.append(self.s2[j - 1])
                    j -= 1
                elif choice == self.T[i - 1, j]:
                    # up
                    xseq.append(self.s1[i - 1])
                    yseq.append('-')
                    i -= 1
                if choice == 0:
                    break

        else:
            i = len(self.s1)
            j = len(self.s2)
            while (i > 0 or j > 0):
                options = [self.T[i - 1, j], self.T[i, j - 1], self.T[i - 1, j - 1]]
                if (i < 1):
                    options.remove(self.T[i - 1, j])
                    options.remove(self.T[i - 1, j - 1])
                if (j < 1):
                    options.remove(self.T[i, j - 1])
                    options.remove(self.T[i - 1, j - 1])
                choice = max(options)
                if choice == self.T[i - 1, j - 1]:
                    # diag
                    xseq.append(self.s1[i - 1])
                    yseq.append(self.s2[j - 1])
                    i -= 1
                    j -= 1
                elif choice == self.T[i - 1, j]:
                    # up
                    xseq.append(self.s1[i - 1])
                    yseq.append('-')
                    i -= 1
                elif choice == self.T[i, j - 1]:
                    # left
                    xseq.append('-')
                    yseq.append(self.s2[j - 1])
                    j -= 1
        xseq.reverse()
        yseq.reverse()
        return xseq, yseq


    def get_num_paths(self, i, j):
        # recrusively gets the number of optimal alignments
        if (self.local and self.T[i, j] == 0) or (i==0) or (j==0):
            return 1
        alignments = 0
        if self.s1[i-1] == self.s2[j-1]:
            alignments += self.get_num_paths(i - 1, j - 1)
        if self.T[i-1, j-1] + self.mismatch_score == self.T[i,j]:
            alignments += self.get_num_paths(i - 1, j - 1)
        if self.T[i, j-1] + self.gap_score == self.T[i,j]:
            alignments += self.get_num_paths(i, j - 1)
        if self.T[i-1,j] + self.gap_score == self.T[i,j]:
            alignments += self.get_num_paths(i - 1, j)
        return alignments



