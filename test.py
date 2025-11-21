"""
Testing ground for random walks
"""

from walk import RandomWalk
from distributions import cauchy, pmf, uniform, exponential, normal
from brw import BranchingRandomWalk

from time import time
import random
import math


RW = RandomWalk(
    step_dist = lambda: 1, 
    angle_dist = pmf({i * math.pi/2 : 1 for i in range(4)}) )

#RW.graph_walk(200, lw=0.75)
#RW.graph_walks(100, 3, lw=0.8, name = "04 - 1d skewed walks")
#RW.graph_MSD(200, 200)
#RW.graph_distribution_1d(500, 1000, -50, 50, 50)

        
BRW = BranchingRandomWalk(
    step_dist = lambda: 1,
    angle_dist = normal(0,1/3), 
    branch_angle_dist = lambda : math.pi/3,
    branch_waiting_dist = exponential(1/3)
    )

"""
#Honeycomb walk
BRW = BranchingRandomWalk(
    step_dist = lambda: 1,
    angle_dist = normal(0,1/8), 
    branch_angle_dist = pmf({-math.pi/3:1, math.pi/3:1}),
    branch_waiting_dist = lambda : 3
    )
"""

n,r,num_sim = 300, 1.5, 1

"""
t0 = time()
for _ in range(num_sim):
    test = BRW.get_barw(n,r)
BRW.graph_walk(0, walk = test)
t1 = time()
print(f"V1 Total time taken: {t1-t0}")
print(" ")


t0 = time()
for _ in range(num_sim):
    test = BRW.get_barw_v2(n,r)
BRW.graph_walk(0, walk = test)
t1 = time()
print(f"V2 Total time taken: {t1-t0}")
print(" ")
"""
t0 = time()
for _ in range(num_sim):
    test = BRW.get_barw_v3(n,r)
BRW.graph_walk(0, walk = test)
t1 = time()
print(f"V3 Total time taken: {t1-t0}")
print(" ")