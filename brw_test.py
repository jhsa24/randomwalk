#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 10:21:41 2025

@author: msp25jha
"""

from walk import Particle
from distributions import uniform, exponential, cauchy
import math
import matplotlib.pyplot as plt
import random

class BranchingRandomWalk:
    def __init__(self, 
                 dimension = 2, 
                 step_dist = lambda : 1, 
                 angle_dist = uniform(-math.pi/5, math.pi/5), 
                 branch_waiting_dist = exponential(1/8),
                 branch_angle_dist = lambda : math.pi / 3):
                
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
        self.branch_waiting_dist = branch_waiting_dist
        self.branch_angle_dist = branch_angle_dist
        self.list_of_walkers = []
        
    def get_walk(self, walker, num_steps):
        positions = [walker.position]
        for _ in range(num_steps):
            walker.rotate(self.angle_dist())
            walker.move(self.step_dist())
            walker.iteration += 1
            positions.append(walker.position)
        return positions
    
    def get_branching_walk_v1(self, num_steps):
        list_of_branches = []
        self.list_of_walkers = [Particle((0,0), 0)]
        walker = self.list_of_walkers[0]
        branch_time = round(self.branch_waiting_dist())
        list_of_branches.append(self.get_walk(walker, branch_time))
        
        list_of_walkers = []
        list_of_walkers.append(Particle(walker.position, 
                                        walker.angle - math.pi / 3, 
                                        walker.iteration))
        walker.rotate(math.pi / 3)
        list_of_walkers.append(walker)
        
        for w in list_of_walkers:
            branch_time = round(self.branch_waiting_dist())
            list_of_branches.append(self.get_walk(w, branch_time))
        
        return list_of_branches
        
    def get_branching_walk_v2(self, num_steps):
        list_of_branches = []
        self.list_of_walkers = [Particle((0,0), 0)]
        
        while self.list_of_walkers:
            current_walkers = self.list_of_walkers.copy()
            for w in current_walkers:
                branch_time = round(self.branch_waiting_dist())
                #Debug!
                print("--- New Branch Starting ---")
                print(f"Current number of branches: {len(self.list_of_walkers)}")
                print(f"Current walker at position {self.list_of_walkers.index(w)} in list")
                w.describe()
                print(f"Branch time: {branch_time}")
                
                if w.iteration + branch_time >= num_steps:
                    print("=============================")
                    print(f"END OF ITERATIONS ALONG BRANCH. Final branch {num_steps - w.iteration} long")
                    list_of_branches.append(self.get_walk(w, num_steps - w.iteration))
                    self.list_of_walkers.remove(w)
                    print(f"final number of iterations: {w.iteration}")
                    print("=============================")
                    continue
            
                list_of_branches.append(self.get_walk(w, branch_time))
                angle = self.branch_angle_dist()
                self.list_of_walkers.remove(w)
                self.list_of_walkers.append(Particle(w.position, w.angle - angle, w.iteration))
                self.list_of_walkers.append(Particle(w.position, w.angle + angle, w.iteration))
                
                print("--- Branch Ending ---")
                w.describe()
        
        return list_of_branches
    
    def graph_walk(self, num_steps, name = None, lw = 0.75):
        branching_walk = self.get_branching_walk_v2(num_steps)
        
        for branch in branching_walk:
            x, y = zip(*branch)
            plt.plot(x, y, color = "black", linewidth = lw)
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
        
        
BRW = BranchingRandomWalk(
    angle_dist = lambda : random.choice([i * math.pi / 2 for i in range(4)]), 
    branch_angle_dist = lambda : math.pi/2,
    branch_waiting_dist = exponential(1/20)
    )


branching_walk_v1 = BRW.get_branching_walk_v1(20)
branching_walk_v2 = BRW.get_branching_walk_v2(20)

BRW.graph_walk(20)
        

