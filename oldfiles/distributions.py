"""
Code that defines various probability distribution functions to be used in random walks
All functions should behave in the same way:
    function(params)()
Should output a random number each time, accorading to some probability distribution
fixed by function(params)

Most functions work on the principle of inverse sampling, ie, take a uniform variable r = (0,1),
and some transformation on r yields the required pdf
"""

import random
import matplotlib.pyplot as plt
import math

#defines a pmf distribution for finite sample spaces
#eg pmf({"H":1, "T"":2}) produces a weighted coin toss
def pmf(prob_dict):
    #calculate sum of dictionary values, ready to renormalise
    total = 0
    for value in prob_dict:
        total += prob_dict[value]
    
    #renormalise probability dictionary
    prob_dict = {value: prob_dict[value] / total for value in prob_dict.keys()}
    
    values = list(prob_dict.keys())
    cumulative = [0]
    #values and cumulative lists together can be used as the cumulative distribution function
    for value in prob_dict:
        cumulative.append(prob_dict[value] + cumulative[-1])
    
    def sample():
        #sample from unif(0,1)
        r = random.random()
        for i in range(len(values)):
            if cumulative[i] < r < cumulative[i+1]:
                return values[i]
        return values[-1]
    #return sample function so that behaviour is same as random.random
    return sample

#defines a cauchy pdf by inverse sampling
def cauchy(location, spread):
    def sample():
        r = random.random()
        return spread * math.tan( math.pi * (r-0.5) ) + location
    return sample

#defines an exponential pdf by inverse sampling
def exponential(scale):
    def sample():
        r = random.random()
        return - (math.log(r)) / scale
    return sample

#I had to cheat on this one :(
def normal(mean, variance):
    return lambda : random.normalvariate(mean, variance)

#stretches random.random to work on any given interval
def uniform(start, end):
    def sample():
        r = random.random()
        return start + (end - start) * r
    return sample


#graphing function for testing purposes, takes a sample size n from a pdf,
#and graphs the distribution of points. Good sanity check.
def graph_pdf_from_sample(pdf, sample_size, axis_min, axis_max, bin_number, name = None):
    sample = [pdf() for _ in range(sample_size)]
    
    
    bin_width = (axis_max - axis_min) / bin_number
    bins = [(axis_min, axis_min + bin_width)]
    
    for i in range(bin_number - 1):
        bin_start = bins[i][1]
        bin_end = bin_start + bin_width
        bins.append((bin_start, bin_end))
    
    count_dict = {}
    
    for b in bins:
        count = 0
        for num in sample:
            if b[0] <= num < b[1]:
                count += 1
            count_dict[0.5 * (b[0]+b[1])] = count
    
    plt.plot(count_dict.keys(), count_dict.values())
    plt.xlabel("Sample value")
    plt.ylabel("Value frequency")
    plt.title(f"Distribution of {sample_size} sample points")
    
    if name:
        plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        plt.show()
    else: plt.show()
    
#graph_pdf_from_sample(cauchy(0,1), 100000, -5, 5, 50)