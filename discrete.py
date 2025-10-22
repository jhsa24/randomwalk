"""
Code for testing discrete random walks, mainly in 2 dimensions

Idea is to use a particle class to keep track of position and rotation
We can then move the particle using its own internal commands,
and append its position to a list each iteration, and graph the list
"""

import random
import numpy as np
import math

from graph import simple_graph

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
        

def random_walk(length, step_dist, angle_dist):
    walker = Particle((0,0), 0)
    positions = []
    for _ in range(length):
        positions.append(walker.position)
        walker.move(step_dist())
        walker.rotate(angle_dist())
    return positions

#simple_graph(random_walk(100, random.random, random.random))