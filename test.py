"""
Testing ground for random walks
"""
from time import time
import math as maths
import matplotlib.pyplot as plt

from aux import swirl, polygon, angle, sample
from rw import RandomWalk
from barw import BranchingRandomWalk
from analysis import WalkData, Analysis
from distributions import cauchy, pmf, uniform, exponential, normal, list_seq
t0 = time()

n = 9
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
    branch_waiting_dist = lambda : 3,
    initial_pos_dist = list_seq(polygon(n, n)),
    initial_angle_dist = list_seq(angle(polygon(n))),
    guidance_strength = 0.1,
    guidance_angle = maths.pi
    )

"""
for num in [1/100, 1/50, 1/20, 1/10, 1/5]:
    
    BRW = BranchingRandomWalk(
        angle_dist = normal(0,1/8), 
        branch_angle_dist = uniform(maths.pi/5, maths.pi/2),
        branch_waiting_dist = exponential(num),
        guidance_strength = 0.1,
        guidance_angle = maths.pi
        )

    
    sim_data = sample(100, BRW.get_barw, 200, 1.5)#, initial_angle = maths.pi)
    data = WalkData(sim_data)
    study = Analysis(data)
    study.graph_walk()
    study.graph_num_walkers()
    print(time() - t0)
"""

sim_data = sample(100, BRW.get_barw, 9, 1.5)#, initial_angle = maths.pi)
data = WalkData(sim_data)
study = Analysis(data)
study.graph_walk()
  
"""
interval = [70, 97]
study.graph_walk(sample = 0, iteration = interval, show = False)
study.graph_walkers_at_t(sample = 0, iteration = interval[1])
"""


"""
import scipy.special as sp
# Evaluate the modified Bessel function of the first kind, zeroth order, at 3.1
value = sp.iv(0, 3.1)
print(value)
"""