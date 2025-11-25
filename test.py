"""
Testing ground for random walks
"""
from time import time
import math as maths

from aux import graph_walks
from rw import RandomWalk
from barw import BranchingRandomWalk
from distributions import cauchy, pmf, uniform, exponential, normal, list_seq

RW = RandomWalk(
    step_dist = lambda: 0.25, 
    angle_dist = uniform(-maths.pi/8, maths.pi/8),
    initial_pos_dist = list_seq([(-3,0), (3,0), (0,3)]),
    initial_angle_dist = list_seq([0, maths.pi, maths.pi*3/2]),
    guidance_strength = 3,
    guidance_angle = 0)

#RW.graph_walk(200, lw=0.75)
#RW.graph_walks(1000, 1000, lw=0.5, name = "03")
#RW.graph_MSD(200, 200)
#RW.graph_distribution_1d(500, 1000, -50, 50, 50)

BRW = BranchingRandomWalk(
    step_dist = lambda: 1,
    angle_dist = normal(0,1/3), 
    branch_angle_dist = lambda : maths.pi/3,
    branch_waiting_dist = exponential(1/10),
    initial_pos_dist = list_seq([(-10,0), (10,0), (0,10)]),
    initial_angle_dist = list_seq([maths.pi, 0, maths.pi/2]),
    guidance_strength = 2,
    guidance_angle = 0
    )

"""
#Honeycomb walk
BRW = BranchingRandomWalk(
    step_dist = lambda: 1,
    angle_dist = normal(0,1/8), 
    branch_angle_dist = pmf({-maths.pi/3:1, maths.pi/3:1}),
    branch_waiting_dist = lambda : 3
    )
"""

walk, somas = BRW.get_multi_barw(100, 3, 1)

graph_walks(walk, somas)
