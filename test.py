"""
Testing ground for random walks
"""
from time import time
import math as maths

from aux import swirl, graph_walks, polygon, angle, sample
from rw import RandomWalk
from barw import BranchingRandomWalk
from graph import WalkData, Analysis
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
    branch_waiting_dist = exponential(1/15),
    initial_pos_dist = list_seq(polygon(n, n)),
    initial_angle_dist = list_seq(angle(polygon(n))),
    guidance_strength = 1,
    guidance_angle = maths.pi
    )

#walk, somas = RW.get_multi_arw(150, n, 1.3)
#walk, somas = BRW.get_multi_barw(200, n, 1.6)

#graph_walks(walk, somas, name = "test", lw = 0.8)

w = sample(500, BRW.get_barw, 100, 1.5, initial_angle = maths.pi)

data = WalkData(w)
analysis = Analysis(data)

analysis.graph_walk()
analysis.graph_angles(0, 2*maths.pi, 60)

t1= time()
print(t1-t0)