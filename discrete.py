"""
Code for testing discrete random walks, mainly in 2 dimensions

Idea is to use a particle class to keep track of position and rotation
We can then move the particle using its own internal commands,
and append its position to a list each iteration, and graph the list
"""

import random
import math
import matplotlib.pyplot as plt
from graph import simple_graph
from distributions import pmf

class Particle:
    def __init__(self, position, angle):
        self.position = position
        self.angle = angle
    
    def move(self, distance):
        new_x = distance * math.sin(self.angle) + self.position[0]
        new_y = distance * math.cos(self.angle) + self.position[1]
        self.position = (new_x, new_y)
        
    def rotate(self, angle):
        self.angle += angle
        
def random_walk_1d(length, step_dist):
    position = 0
    positions = [0]
    for _ in range(length - 1):
        position += step_dist()
        positions.append(position)
    return positions
    
def random_walk_2d(length, step_dist, angle_dist):
    walker = Particle((0,0), 0)
    positions = [(0,0)]
    for _ in range(length):
        walker.move(step_dist())
        walker.rotate(angle_dist())
        positions.append(walker.position)
    return positions

"""
basic_walk = random_walk_2d(1000, lambda: 1, lambda: random.choice([math.pi * i/2 for i in range(4)]))
simple_graph(basic_walk)

distances = []
for position in basic_walk:
    distances.append((position[0]**2 + position[1]**2) ** 1/2)
"""

walk_1d = random_walk_1d(100, pmf({-1:3,1:4}))
#plt.plot(walk_1d)

num_walks = 500000

final_positions = []
for _ in range(num_walks):
    #plt.plot(random_walk_1d(8000, pmf({-1:1, 1:1})), linewidth=0.2)
    final_positions.append(random_walk_1d(100, pmf({-1:1, 1:1}))[-1])
    final_positions.append(random_walk_1d(101, pmf({-1:1, 1:1}))[-1])
#plt.title(f"{num_walks} 1d Random Walks")
#plt.xlabel("Time")
#plt.ylabel("Position of walker")
#plt.savefig(f"plots/{num_walks} walks.png", dpi=300, bbox_inches='tight')
    
d = { n:0 for n in range(-40,41,1)}
for pos in final_positions:
    if pos in d:
        d[pos] += 1

counts = []
for key in d:
    counts.append(d[key])
    
plt.plot(counts)






