"""
Testing ground for random walks
"""
from time import time
t0 = time()
import math as maths

from aux import swirl, polygon, angle, sample, vonmis
from rw import RandomWalk
from barw import BranchingRandomWalk
from analysis import WalkData, Analysis
from distributions import cauchy, pmf, uniform, exponential, normal, list_seq

n = 100

RW = RandomWalk(
    angle_dist=normal(0, 1/3),
    initial_pos_dist=list_seq(polygon(n, n)),
    initial_angle_dist=list_seq(angle(polygon(n))),
    guidance_strength=1,
    guidance_angle=1
)

BRW = BranchingRandomWalk(
    angle_dist=uniform(-maths.pi/10, maths.pi/10),
    branch_angle_dist=uniform(maths.pi/5, maths.pi/3),
    branch_prob = 0.1,
    guidance_strength=0.1,
    guidance_angle= maths.pi
)

#sim_data = sample(200, BRW.get_barw, 150, 1.5, initial_angle = maths.pi)
#sim_data = sample(30, RW.get_rw, 150)
#data = WalkData(sim_data)
#data.save("BARW-pb-0.1")

study = Analysis("BARW-pb-0.1")


#study.graph_MSD(x_log = True, y_log = True)
#study.graph_branch_lengths(y_log = True, bin_number = 50)

study.graph_angles(bin_number = 72, show = False, sample = None)

import matplotlib.pyplot as plt
import numpy as np

pb = 0.1
ae = maths.pi/10
ab = maths.pi/3
fc = 0.1


D = 1/6 * (pb * (ab**2 + ab*ae) + ae**2)
mu = 1/2 * (pb * ab + ae)
v = (mu * fc) / D

print(f"Diffusion: {D}, mobility: {mu}")

x = np.linspace(0, 2*maths.pi, num = 50)
y = vonmis(x - np.pi, v)
plt.plot(x, y, color = "red")

print(time() - t0)