#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 15:19:12 2025

@author: msp25jha
"""

"""
Tests for interacting random walks

The annhilating part of the BARW code isn't working, so here is a testing ground
for simplified ideas
"""

from walk import Particle
from distributions import uniform, exponential
import matplotlib.pyplot as plt
import math


class InteractingWalk:
    def __init__(self, 
                 dimension = 2, 
                 step_dist = lambda : 1, 
                 angle_dist = uniform(-math.pi/5, math.pi/5), 
                 branch_waiting_dist = exponential(1/15),
                 branch_angle_dist = lambda : math.pi / 3):
                
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
        self.branch_waiting_dist = branch_waiting_dist
        self.branch_angle_dist = branch_angle_dist

    def get_arw(self, num_steps, radius):
        w = Particle((0,0), 0)
        positions = [w.position]
        
        for _ in range(num_steps):
            dead = False
            
            for pos in positions:
                dx = pos[0] - w.position[0]
                dy = pos[1] - w.position[1]
                distance = (dx**2 + dy**2) ** 0.5
                if w.position != pos and distance < radius:
                    dead = True
                    break
            if dead:
                break
        
            w.rotate(self.angle_dist())
            w.move(self.step_dist())
            positions.append(w.position)
        
        return positions
            
    def graph_walk_v1(self, walk, name = None, lw = 0.75):
        #get_walk outputs tuples, so we need to unpack them to separate x&y
        x, y = zip(*walk)
        plt.plot(x, y, linewidth = lw)
        plt.title("2d walk")
        plt.gca().set_aspect('equal')

        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
        
    def graph_walk_v2(self, walk, name = None, lw = 0.75):
        
        for sub_walk in walk:
            x, y = zip(*sub_walk)
            plt.plot(x, y, color = "black", linewidth = lw)
        plt.title("2d walk")
        plt.gca().set_aspect('equal')
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()

IW = InteractingWalk()
walk = IW.get_arw(200, 0.9)

IW.graph_walk_v1(walk) 

"""
lw = 1
line = [(x,5) for x in range(-5,5)]
x, y = zip(*line)
plt.plot(x, y, color = "red", linewidth = lw)
"""




"""
for branch in branching_walk:
    x, y = zip(*branch)
    plt.plot(x, y, color = "black", linewidth = lw)
"""