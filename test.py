from walk import RandomWalk
from distributions import cauchy, pmf, uniform

import math


RW = RandomWalk(step_dist = lambda: 1, angle_dist = uniform(- math.pi/10, math.pi / 8) )#pmf({i * math.pi/2 : 1 for i in range(4)}) )

RW.graph_walk(1000, lw=0.75)
#RW.graph_walks(100, 3, lw=0.8, name = "04 - 1d skewed walks")
RW.graph_MSD(1000, 200)
#RW.graph_distribution_1d(500, 1000, -50, 50, 50)