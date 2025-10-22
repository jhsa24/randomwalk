"""
Code that defines various probability distribution functions to be used in random walks
"""

import random
import matplotlib.pyplot as plt
import math

y = [random.random() for i in range(50000)]


y = [math.tan( math.pi*(i-0.5)) for i in y]

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

plt.plot(bin_mid, counts)