"""
Testing ground for random walks
"""
from time import time
import math as maths

from aux import swirl, polygon, angle, sample
from rw import RandomWalk
from barw import BranchingRandomWalk
from analysis import WalkData, Analysis
from distributions import cauchy, pmf, uniform, exponential, normal, list_seq

n = 3
t0 = time()

RW = RandomWalk( 
    angle_dist = uniform(-maths.pi/8, maths.pi/8),
    initial_pos_dist = list_seq(polygon(n, n)),
    initial_angle_dist = list_seq(angle(polygon(n))),
    guidance_strength = 1,
    guidance_angle = swirl(0)
    )

BRW = BranchingRandomWalk(
    angle_dist = normal(0,1/8), 
    branch_angle_dist = uniform(maths.pi/5, maths.pi/2),
    branch_waiting_dist = exponential(1/20),
    initial_pos_dist = list_seq(polygon(n, n)),
    initial_angle_dist = list_seq(angle(polygon(n))),
    guidance_strength = 2,
    guidance_angle = maths.pi
    )


w = sample(1000, BRW.get_barw, 150, 1.5)#, initial_angle = maths.pi)

data = WalkData(w)
analysis = Analysis(data)

#analysis.graph_walk()
analysis.graph_angles(bin_number = 1000)

t1= time()
print(t1-t0)