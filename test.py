from walk import RandomWalk
from distributions import cauchy, pmf, uniform, exponential, normal
from brw_test import BranchingRandomWalk

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
    step_dist = cauchy(0,1),
    angle_dist = pmf({i * math.pi/2 : 1 for i in range(4)}), 
    branch_angle_dist = lambda : math.pi/2,
    branch_waiting_dist = exponential(1/10)
    )
BRW = BranchingRandomWalk()

BRW.graph_walk(50, name = "Branching Random Walk")

#branching_walk_v3 = BRW.get_branching_walk_v3(20)
