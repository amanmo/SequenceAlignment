import pandas as pd
import time
import sys
import os, psutil
from memory_profiler import memory_usage

class SequenceAlignment:

    def __init__(self, filename='input.txt'):
        'Function to parse the input file, generate base strings and set penalties'

        with open(filename) as f:
            lines = f.readlines()

        lines = list(map(str.strip, lines))

        #Getting Base Strings
        base_X = lines[0]
        for i in range(1, len(lines)):
            if lines[i].isalpha():
                base_Y = lines[i]
                j = lines[1:i]
                k = lines[i+1:] if i+1 < len(lines) else []
                break

        #Generating Full Strings
        self.string_X = self.generateString(base_X, j)
        self.string_Y = self.generateString(base_Y, k)

        #String Length Validation
        if len(self.string_X) != (2**len(j) * len(base_X)) or len(self.string_Y) != (2**len(k) * len(base_Y)):
            print("Error: Error in String Generation")
            exit()

        #Initializing Penalty Values
        self.penalty_gap = 30
        self.penalty_alpha = pd.DataFrame([[0, 110, 48, 94], [110, 0, 118, 48], [48, 118, 0, 110], [94, 48, 110, 0]], index=['A','C','G','T'], columns=['A','C','G','T'])

    @staticmethod
    def generateString(base, params):
        'Function to generate full string from base string and parameters'

        for x in params:
            x = int(x)
            left = base[:x+1]
            right = base[x+1:] if x+1 < len(base) else ''
            base = left + base + right
        return base

    def align(self):
        'Function that performs the sequence alignment'

        start_time = time.time()

        m = len(self.string_X)
        n = len(self.string_Y)

        self.M = [[float('inf') for _ in range(n+1)] for _ in range(m+1)]
        
        #Base Cases
        for i in range(m+1):
            self.M[i][0] = i * self.penalty_gap

        for j in range(n+1):
            self.M[0][j] = j * self.penalty_gap

        #Recurrence
        for row in range(1, m+1):
            for col in range(1, n+1):
                self.M[row][col] = min(self.M[row-1][col-1] + self.penalty_alpha[self.string_X[row-1]][self.string_Y[col-1]],
                                  self.M[row-1][col] + self.penalty_gap,
                                  self.M[row][col-1] + self.penalty_gap)

        #Trace backwards from self.M[m][n]
        output_X = ''
        output_Y = ''
        x, y = m, n
        while(x!=0 or y!=0):
            current = self.M[x][y]

            #Finding how the minimum was obtained
            #1. Mismatch
            if (x-1 > -1 and y-1 > -1) and (self.M[x-1][y-1] + self.penalty_alpha[self.string_X[x-1]][self.string_Y[y-1]] == current):
                output_X = self.string_X[x-1] + output_X
                output_Y = self.string_Y[y-1] + output_Y
                x-=1
                y-=1
                continue

            #2. Gap
            if x-1 > -1 and (self.M[x-1][y] + self.penalty_gap == current):
                output_X = self.string_X[x-1] + output_X
                output_Y = '_' + output_Y
                x-=1
            elif y-1 > -1 and (self.M[x][y-1] + self.penalty_gap == current):
                output_X = '_' + output_X
                output_Y = self.string_Y[y-1] + output_Y
                y-=1

        cost = self.M[m][n]
        total_time = time.time() - start_time
        return output_X, output_Y, cost, total_time

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

    #Option 1: Subtraction  (Better Results)
    start_memory = (psutil.Process(os.getpid()).memory_info().rss)//1024
    output_X, output_Y, cost, total_time = s.align()
    memory = (psutil.Process(os.getpid()).memory_info().rss)//1024 - start_memory

    #Option 2: Memory Profiler - Memory Usage
    # memories, retVal = memory_usage(s.align, retval=True)
    # output_X, output_Y, cost, total_time = retVal
    # memory = max(memories)

    s.generateOutput(output_X, output_Y, cost, total_time, memory)