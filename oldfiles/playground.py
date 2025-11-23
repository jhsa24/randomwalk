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

    def get_arw(self, num_steps, radius = 1):
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
    
    def get_multi_arw(self, num_steps, radius = 1, num_walkers = 2):
        wd = {}
        
        for i in range(num_walkers):
            rx = uniform(-5, 5)()
            ry = uniform(-5, 5)()
            ra = uniform(0, math.pi*2)()
            #walker_list.append([Particle((rx,ry), ra), [(rx, ry)] ])
            wd[i] = {"walker":Particle((rx,ry), ra),
                     "dead": False,
                     "positions": [(rx, ry)]}
        print(wd)
        print()
        iteration = 0
        while iteration < num_steps:
            current_walkers = wd.copy()
            for i in current_walkers:
                if wd[i]["dead"]:
                    continue
                #test distance from other paths
                for j in wd:
                    for pos in wd[j]["positions"]:
                        dx = wd[i]["walker"].position[0] - pos[0]
                        dy = wd[i]["walker"].position[1] - pos[1]
                        distance = (dx**2 + dy**2)**0.5
                        if 0 < distance < radius:
                            wd[i]["dead"] = True
                            break
                    if wd[i]["dead"]:
                        break
                if wd[i]["dead"]:
                    continue
                
                w = wd[i]["walker"]
                w.rotate(self.angle_dist())
                w.move(self.step_dist())
                w.iteration += 1
                wd[i]["positions"].append(w.position)
            if iteration in [0,1]:
                print(iteration)
                print(wd)
                print()
            
            iteration += 1
        
        positions = [wd[i]["positions"] for i in wd]
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


IW = InteractingWalk(step_dist = lambda : 1, 
                     angle_dist = uniform(-math.pi/5, math.pi/5), 
                     branch_waiting_dist = exponential(1/15), 
                     branch_angle_dist = lambda : math.pi / 3)


walk = IW.get_multi_arw(100, radius=0.5, num_walkers=2)

IW.graph_walk_v2(walk) 




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