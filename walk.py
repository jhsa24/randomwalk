"""
Code for testing discrete random walks, mainly in 2 dimensions

Idea is to use a particle class to keep track of position and rotation
We can then move the particle using its own internal commands,
and append its position to a list each iteration, and graph the list
"""

import random
import math
import matplotlib.pyplot as plt

from aux import mean
from graph import simple_graph_1d, simple_graph_2d
from distributions import pmf, cauchy

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
        
class RandomWalk:
    def __init__(self, dimension = 2, step_dist = lambda: 1, angle_dist = lambda: random.choice([math.pi * i/2 for i in range(4)])):
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
        self.particle = Particle((0,0), 0)
        
    def get_walk(self, num_steps):
        positions = []
        if self.dimension == 1:
            positions.append(0)
            for _ in range(num_steps):
                positions.append( positions[-1] + self.step_dist() )
        
        if self.dimension == 2:
            self.particle.position = (0,0)
            positions.append(self.particle.position)
            for _ in range(num_steps):
                self.particle.rotate(self.angle_dist())
                self.particle.move(self.step_dist())
                positions.append(self.particle.position)
        
        return positions    
    
    def graph_walk(self, num_steps, name = None, lw = 0.5):
        if self.dimension == 1:
            simple_graph_1d(self.get_walk(num_steps), name, lw)
        if self.dimension == 2:
            simple_graph_2d(self.get_walk(num_steps), name, lw)
    
    def graph_walks(self, num_steps, num_walks, name = None, lw = 0.5):
        if self.dimension == 1:
            for _ in range(num_walks):
                plt.plot(self.get_walk(num_steps), linewidth = lw)
            plt.xlabel("Number of steps")
            plt.ylabel("Position of walker")
            plt.title(f"{num_walks} 1d random walks")
            if name:
                plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            else: plt.show()
    
    def graph_distribution_1d(self, num_steps, num_walks, axis_min, axis_max, name = None, lw = 1):
        final_positions = []
        num_iter = int(num_walks / 2)
        frequency_dict = {n:0 for n in range(axis_min, axis_max)}
        
        for i in range(num_iter):
            final_positions.append(self.get_walk(num_steps)[-1])
            final_positions.append(self.get_walk(num_steps + 1)[-1])
        
        for position in final_positions:
            if position in frequency_dict:
                frequency_dict[position] += 1
                
        plt.plot(frequency_dict.keys(), frequency_dict.values(), linewidth = lw)
        plt.xlabel(f"Position of walker after {num_steps} steps")
        plt.ylabel("Frequency")
        plt.title(f"Distribution of positions reached by {num_walks} 1d random walks")
    
    def __graph_MSD_1d(self, num_steps, num_walks, name = None, lw = 1):
        list_of_walks, list_of_averages = [], []
        for _ in range(num_walks):
            list_of_walks.append(self.get_walk(num_steps))
        
        for t in range(num_steps):
            positions_squared = [abs(list_of_walks[i][t])**2 for i in range(num_walks)]
            list_of_averages.append(mean(positions_squared))
        
        plt.plot(list_of_averages, linewidth = lw)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Squared Distance")
        plt.title(f"MSD of {num_walks} 1d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        else: plt.show()
    
    def __graph_MSD_2d(self, num_steps, num_walks, name = None, lw = 1):
        list_of_walks, list_of_averages = [], []
        for _ in range(num_walks):
            list_of_walks.append(self.get_walk(num_steps))
        
        for t in range(num_steps):
            positions = [list_of_walks[i][t] for i in range(num_walks)]
            distances_squared = [pos[0]**2 + pos[1]**2 for pos in positions]
            list_of_averages.append(mean(distances_squared))
            """
            if t%3 == 0:
                print(f" TIME: {t}")
                print(f"positions: {positions}")
                print(f"distances: {distances_squared}")
                print(f"mean: {mean(distances_squared)}")
                print()
            """
        
        plt.plot(list_of_averages, linewidth = lw)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Squared Distance")
        plt.title(f"MSD of {num_walks} 2d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        else: plt.show()
    
    def graph_MSD(self, num_steps, num_walks, name = None, lw = 1):
        if self.dimension == 1:
            self.__graph_MSD_1d(num_steps, num_walks, name, lw)
        if self.dimension == 2:
            self.__graph_MSD_2d(num_steps, num_walks, name, lw)
            
    
    
walker = RandomWalk(1, step_dist = pmf({-1:1, 1:1}))
 
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
RW = RandomWalk(2, step_dist = lambda : 1, angle_dist = cauchy(0,1))

RW.graph_walk(500, lw=0.5)
#RW.graph_MSD(100, 200)


