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
t0 = time()

n = 100
RW = RandomWalk(
    angle_dist=normal(0, 1/3),
    initial_pos_dist=list_seq(polygon(n, n)),
    initial_angle_dist=list_seq(angle(polygon(n))),
    guidance_strength=0,
    guidance_angle=swirl(0)
)

BRW = BranchingRandomWalk(
    angle_dist=normal(0, 1/8),
    branch_angle_dist=uniform(maths.pi/5, maths.pi/2),
    branch_waiting_dist=exponential(1/55),
    initial_pos_dist=list_seq(polygon(n, n)),
    initial_angle_dist=list_seq(angle(polygon(n))),
    guidance_strength=0.1,
    guidance_angle=0.5
)

#sim_data = sample(10, BRW.get_barw, 500, 1.6)
#sim_data = sample(30, RW.get_rw, 500)
#data = WalkData(sim_data)
#data.save("TEST4")

study = Analysis("TEST4.pkl")
study.graph_walk(sample=0)



study.animate_walk(sample = 17,
                   name = "TEST",
                   show_line = False,
                   m_size = 50,
                   lw = 0.7
                   )

print(time() - t0)

"""
#an attempt at animation!!!
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
#constants
sample = 0
max_iter = 120

positions = study.positions
x_min, x_max = positions[:,0].min(), positions[:,0].max()
y_min, y_max = positions[:,1].min(), positions[:,1].max()

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_facecolor('white')

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

branches = np.unique(study.branch_id[study.sample_id == sample])
lines = []

for b in branches:
    line, = ax.plot([], [], color = 'black')
    lines.append(line)

active_walkers_scatter = ax.scatter([], [], color = 'red', marker='o', s = 50)
    

def update(frame):
    m_sample = study.sample_id == sample
    m_iter = study.iters <= frame
    
    for line, b in zip(lines, branches):
        pts = study.positions[(study.branch_id == b) & (m_sample) & (m_iter)]
        if len(pts) > 0:
            line.set_data(pts[:,0], pts[:,1])
        else:
            line.set_data([], [])
        
        active_walkers = study.positions[(study.sample_id == sample) & (study.iters == frame)]
        active_walkers_scatter.set_offsets(active_walkers)
    
    return lines + [active_walkers_scatter]

ani = FuncAnimation(fig, 
                    update, 
                    frames = range(0, study.iters[study.sample_id == sample].max()), 
                    blit=True
                    )

ani.save('test_animation.mp4', fps = 30, dpi = 300)
"""


"""
import scipy.special as sp
# Evaluate the modified Bessel function of the first kind, zeroth order, at 3.1
value = sp.iv(0, 3.1)
print(value)
"""
