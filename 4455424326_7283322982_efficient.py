import pandas as pd
import time
import sys
import math
import os, psutil

class SequenceAlignment:

    def __init__(self, filename='input.txt'):
        'Function to parse the input file, generate base strings and set penalties'

        with open(filename) as f:
            lines = f.readlines()

        lines = list(map(str.strip, lines))

        # Getting Base Strings
        base_X = lines[0]
        for i in range(1, len(lines)):
            if lines[i].isalpha():
                base_Y = lines[i]
                j = lines[1:i]
                k = lines[i + 1:] if i + 1 < len(lines) else []
                break

        # Generating Full Strings
        self.string_X = self.generateString(base_X, j)
        self.string_Y = self.generateString(base_Y, k)

        # String Length Validation
        if len(self.string_X) != (2 ** len(j) * len(base_X)) or len(self.string_Y) != (2 ** len(k) * len(base_Y)):
            print("Error: Error in String Generation")
            exit()

        # Initializing Penalty Values
        self.penalty_gap = 30
        self.penalty_alpha = pd.DataFrame([[0, 110, 48, 94],
                                           [110, 0, 118, 48],
                                           [48, 118, 0, 110],
                                           [94, 48, 110, 0]],
                                          index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])

    @staticmethod
    def generateString(base, params):
        'Function to generate full string from base string and parameters'

        for x in params:
            # Check for validity of parameter? (x could be more than length of string)
            x = int(x)
            left = base[:x + 1]
            right = base[x + 1:] if x + 1 < len(base) else ''
            base = left + base + right
        return base

    def align(self, string_X = None, string_Y = None):
        'Function that performs the sequence alignment'

        if string_X == None and string_Y == None:
            string_X = self.string_X
            string_Y = self.string_Y
            
        m = len(string_X)
        n = len(string_Y)

        M = [[float('inf') for _ in range(n + 1)] for _ in range(m + 1)]

        # Base Cases
        for i in range(m + 1):
            M[i][0] = i * self.penalty_gap

        for j in range(n + 1):
            M[0][j] = j * self.penalty_gap

        # Recurrence
        for row in range(1, m + 1):
            for col in range(1, n + 1):
                M[row][col] = min(
                    M[row - 1][col - 1] + self.penalty_alpha[string_X[row - 1]][string_Y[col - 1]],
                    M[row - 1][col] + self.penalty_gap,
                    M[row][col - 1] + self.penalty_gap)

        # Trace backwards from M[m][n]
        output_X = ''
        output_Y = ''
        x, y = m, n
        while (x != 0 or y != 0):
            current = M[x][y]

            # Finding how the minimum was obtained
            # 1. Mismatch
            if (x - 1 > -1 and y - 1 > -1) and (
                    M[x - 1][y - 1] + self.penalty_alpha[string_X[x - 1]][string_Y[y - 1]] == current):
                output_X = string_X[x - 1] + output_X
                output_Y = string_Y[y - 1] + output_Y
                x -= 1
                y -= 1
                continue

            # 2. Gap
            if x - 1 > -1 and (M[x - 1][y] + self.penalty_gap == current):
                output_X = string_X[x-1] + output_X
                output_Y = "_" + output_Y
                x -= 1
            elif y - 1 > -1 and (M[x][y - 1] + self.penalty_gap == current):
                output_X = "_" + output_X
                output_Y = string_Y[y-1] + output_Y
                y -= 1

        cost = M[m][n]
        return output_X, output_Y, cost

    def align_cost(self, string_X = None, string_Y = None):
        'Function that performs the sequence alignment'

        if string_X == None and string_Y == None:
            string_X = self.string_X
            string_Y = self.string_Y

        m = len(string_X)
        n = len(string_Y)

        if m == 0:
            return n * self.penalty_gap
        if n == 0:
            return m * self.penalty_gap

        M = [[float('inf') for _ in range(2)] for _ in range(m + 1)]

        #M dimension i = m, j = 2
        # Base Cases
        for i in range(m + 1):
            M[i][0] = i * self.penalty_gap

        for j in range(2):
            M[0][j] = j * self.penalty_gap

        j = 1

        # Recurrence
        for col in range(1, n + 1):
            M[0][1] = col * self.penalty_gap
            for row in range(1, m + 1):
                M[row][1] = min(
                    M[row - 1][0] + self.penalty_alpha[string_X[row - 1]][string_Y[col - 1]],
                    M[row - 1][1] + self.penalty_gap,
                    M[row][0] + self.penalty_gap)
            for row in range(0, m + 1):
                M[row][0] = M[row][1]

        return(M)

    def divide_and_conquer(self, string_X, string_Y):

        m = len(string_X)
        n = len(string_Y)

        if m == 0 or n == 0 or m == 1:
            output_X, output_Y, _ = self.align(string_X, string_Y)    # align of last two rows
            return output_X, output_Y

        min = float('inf')
        j = 0

        X_r = string_X[math.floor(m/2):]
        X_r = X_r[::-1]
        A = self.align_cost(string_Y = string_X[:math.floor(m/2)],string_X= string_Y)
        B = self.align_cost(string_Y =X_r, string_X=string_Y[::-1])
        B = B[::-1]
        for i in range(n):
            split_cost = A[i][1] + B[i][1]
            if split_cost < min:
                min = split_cost
                j = i

        if j < n:
            Y_to_push = string_Y[j]
        else:
            Y_to_push = ''


        if string_X[math.floor(m/2)] == Y_to_push:
            left_half_X, left_half_Y= self.divide_and_conquer( string_X[ :math.floor(m/2)] , string_Y[:j])
            right_half_X, right_half_Y = self.divide_and_conquer(string_X[ math.floor(m/2)+1:] , string_Y[j+1:])

            output_X = left_half_X + string_X[math.floor(m/2)] + right_half_X
            output_Y = left_half_Y + string_Y[  j   ] + right_half_Y
            return output_X, output_Y


        else:
            if j==0:
                left_half_X, left_half_Y = self.divide_and_conquer(string_X[:math.floor(m/2)], '')
                right_half_X, right_half_Y = self.divide_and_conquer(string_X[math.floor(m/2):], string_Y[j:])
            else:
                left_half_X, left_half_Y = self.divide_and_conquer(string_X[:math.floor(m/2)], string_Y[:j])
                right_half_X, right_half_Y = self.divide_and_conquer(string_X[math.floor(m/2):], string_Y[j:])

            output_X = left_half_X + right_half_X
            output_Y = left_half_Y + right_half_Y

            return output_X, output_Y


    def generateOutput(self, output_X, output_Y, cost, total_time, memory):
        'Function to generate the output file'

        lines = []
        if len(output_X) <= 50:
            lines += [output_X + ' ' + output_X]
        else:
            lines += [output_X[:50] + ' ' + output_X[len(output_X)-50:]]

        if len(output_Y) <= 50:
            lines += [output_Y + ' ' + output_Y]
        else:
            lines += [output_Y[:50] + ' ' + output_Y[len(output_Y)-50:]]

        lines += [str(cost), str(total_time), str(memory)]

        with open('output.txt', 'w') as f:
            f.write('\n'.join(lines))

if __name__ == '__main__':
    try:
        _, filename = sys.argv
        s = SequenceAlignment(filename)
    except:
        s = SequenceAlignment()

    start_time = time.time()
    start_memory = (psutil.Process(os.getpid()).memory_info().rss) // 1024
    
    M = s.align_cost(s.string_X, s.string_Y)
    cost = M[len(s.string_X)][1]
    output_X, output_Y = s.divide_and_conquer(s.string_X, s.string_Y)

    memory = (psutil.Process(os.getpid()).memory_info().rss) // 1024 - start_memory
    total_time = time.time() - start_time
    
    s.generateOutput(output_X, output_Y, cost, total_time, memory)