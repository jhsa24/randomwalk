"""
Testing ground for random walks
"""
from time import time
import random
import math as maths

from aux import graph_walks
from rw import RandomWalk
from barw import BranchingRandomWalk
from distributions import cauchy, pmf, uniform, exponential, normal

RW = RandomWalk(
    step_dist = lambda: 1, 
    angle_dist = pmf({i * maths.pi/2 : 1 for i in range(4)}) )

#RW.graph_walk(200, lw=0.75)
#RW.graph_walks(100, 3, lw=0.8, name = "04 - 1d skewed walks")
#RW.graph_MSD(200, 200)
#RW.graph_distribution_1d(500, 1000, -50, 50, 50)

BRW = BranchingRandomWalk(
    step_dist = lambda: 1,
    angle_dist = normal(0,1/3), 
    branch_angle_dist = lambda : maths.pi/3,
    branch_waiting_dist = exponential(1/3)
    )


#Honeycomb walk
BRW = BranchingRandomWalk(
    step_dist = lambda: 1,
    angle_dist = normal(0,1/8), 
    branch_angle_dist = pmf({-maths.pi/3:1, maths.pi/3:1}),
    branch_waiting_dist = lambda : 3
    )

graph_walks(BRW.get_barw(100, 1.5))