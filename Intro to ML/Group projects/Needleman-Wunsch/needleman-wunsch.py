from Bio import pairwise2
import random
import time
from Bio.pairwise2 import format_alignment
class SequenceAlignment(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.solution = []

    delta = lambda self, x, y, i, j: 1 if x[i] != y[j] else 0

    def find_solution(self, OPT, m, n):
        if m == 0 and n == 0:
            return
        insert = OPT[m][n - 1] + 1 if n != 0 else float("inf")
        align = (
            OPT[m - 1][n - 1] + self.delta(self.x, self.y, m - 1, n - 1)
            if m != 0 and n != 0
            else float("inf")
        )
        delete = OPT[m - 1][n] + 1 if m != 0 else float("inf")

        best_choice = min(insert, align, delete)

        if best_choice == insert:
            self.solution.append("insert_" + str(self.y[n - 1]))
            return self.find_solution(OPT, m, n - 1)

        elif best_choice == align:
            self.solution.append("align_" + str(self.y[n - 1]))
            return self.find_solution(OPT, m - 1, n - 1)

        elif best_choice == delete:
            self.solution.append("remove_" + str(self.x[m - 1]))
            return self.find_solution(OPT, m - 1, n)

    def alignment(self):
        n = len(self.y)
        m = len(self.x)
        OPT = [[0 for i in range(n + 1)] for j in range(m + 1)]

        for i in range(1, m + 1):
            OPT[i][0] = i

        for j in range(1, n + 1):
            OPT[0][j] = j

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                OPT[i][j] = min(
                    OPT[i - 1][j - 1] + self.delta(self.x, self.y, i - 1, j - 1),
                    OPT[i - 1][j] + 1,
                    OPT[i][j - 1] + 1,
                )

        self.find_solution(OPT, m, n)

        return (OPT[m][n], self.solution[::-1])
def fill_dyn_matrix(x, y):
    L = [[0]*(len(y)+1) for _ in range(len(x)+1)]
    for x_i,x_elem in enumerate(x):
        for y_i,y_elem in enumerate(y):
            if x_elem == y_elem:
                L[x_i][y_i] = L[x_i-1][y_i-1] + 1
            else:
                L[x_i][y_i] = max((L[x_i][y_i-1],L[x_i-1][y_i]))
    return L
def LCS_DYN(x, y):
    L = fill_dyn_matrix(x, y)
    LCS = []
    x_i,y_i = len(x)-1,len(y)-1
    while x_i >= 0 and y_i >= 0:
        if x[x_i] == y[y_i]:
            LCS.append(x[x_i])
            x_i, y_i = x_i-1, y_i-1
        elif L[x_i-1][y_i] > L[x_i][y_i-1]:
            x_i -= 1
        else:
            y_i -= 1
    LCS.reverse()
    return LCS
if __name__ == '__main__':
    symbols = ["A", "G", "T", "C"]
    x = 'AGCTCT'
    x = "".join([random.choice(symbols) for i in range(900)])
    y = "".join([random.choice(symbols) for i in range(900)])
    print('We we want to transform: ' + x + ' to: ' + y)
    print(LCS_DYN(x, y))
    start = time.time()
    sqalign = SequenceAlignment(x, y)
    print("SequenceAlignment:", start - time.time())
    start = time.time()
    min_edit, steps = sqalign.alignment()
    print('Minimum amount of edit steps are: ' + str(min_edit))
    print('And the way to do it is: ' + str(steps))
    alignments = pairwise2.align.localxx(x, y)
    for a in alignments:
        print(format_alignment(*a))
    print("pairwise2",start - time.time())