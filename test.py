"""
Testing ground for random walks
"""

from walk import RandomWalk
from distributions import cauchy, pmf, uniform, exponential, normal
from brw import BranchingRandomWalk

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
    step_dist = uniform(0,3),
    angle_dist = normal(0,1/5), 
    branch_angle_dist = uniform(math.pi/5, math.pi/3),
    branch_waiting_dist = exponential(1/15)
    )


BRW.graph_walk(50, name = "Branching Random Walk")


