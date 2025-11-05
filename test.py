from walk import RandomWalk
from distributions import cauchy, pmf, uniform

import math


RW = RandomWalk(1, step_dist = pmf({-1:1, 1:1}), angle_dist = uniform(0, 2*math.pi))#pmf({i * math.pi/2 : 1 for i in range(4)}) )

RW.graph_walk(500, lw=0.75, name = "2d walk")
#RW.graph_walks(100, 3, lw=0.8, name = "04 - 1d skewed walks")
RW.graph_MD(500, 200, name = "2d MSD")
#RW.graph_distribution_1d(100, 1000, -15, 65, name = "05 - 1d skewed position distrubution")