from walk import RandomWalk
from distributions import cauchy


RW = RandomWalk(2, step_dist = lambda : 1)#, angle_dist = cauchy(0,1) )

RW.graph_walk(200, lw=0.75, name = "2d walk")
#RW.graph_walks(100, 3, lw=0.8, name = "04 - 1d skewed walks")
RW.graph_MSD(500, 200, name = "2d MSD")
#RW.graph_distribution_1d(100, 1000, -15, 65, name = "05 - 1d skewed position distrubution")