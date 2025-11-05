"""
Code that defines various probability distribution functions to be used in random walks
"""

import random
import matplotlib.pyplot as plt
import math


def pmf(prob_dict):
    total = 0
    for value in prob_dict:
        total += prob_dict[value]
    
    prob_dict = {value: prob_dict[value] / total for value in prob_dict.keys()}
    values = list(prob_dict.keys())
    cumulative = [0]
    
    for value in prob_dict:
        cumulative.append(prob_dict[value] + cumulative[-1])
    
    def sample():
        r = random.random()
        for i in range(len(values)):
            if cumulative[i] < r < cumulative[i+1]:
                return values[i]
        return values[-1]
    return sample

def cauchy(location, spread):
    def sample():
        r = random.random()
        return spread * math.tan( math.pi * (r-0.5) ) + location
    return sample

def uniform(start, end):
    def sample():
        r = random.random()
        return start + (end - start) * r
    return sample

y = [cauchy(0,1)() for _ in range(500)]

y = sorted(y)

interval = [-6,6]
bin_number = 50

bin_width = (interval[1] - interval[0]) / bin_number

bins = [(interval[0], interval[0] + bin_width)]
bin_mid = []
counts = []
for i in range(bin_number - 1):
    bin_start = bins[i][1]
    bin_end = bin_start + bin_width
    bins.append((bin_start, bin_end))
    

for b in bins:
    count = 0
    for num in y:
        if b[0] < num < b[1]:
            count += 1
    counts.append(count)
    bin_mid.append(0.5 * (b[0] + b[1]))

#plt.plot(bin_mid, counts)
