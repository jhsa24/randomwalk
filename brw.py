#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 10:21:41 2025

@author: msp25jha
"""

from walk import Particle
from distributions import uniform, exponential
import math
import matplotlib.pyplot as plt



class BranchingRandomWalk:
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
        self.list_of_walkers = []
        
    #performs simple random walk, returns list of positions
    def get_walk(self, walker, num_steps):
        positions = [walker.position]
        for _ in range(num_steps):
            walker.rotate(self.angle_dist())
            walker.move(self.step_dist())
            walker.iteration += 1
            positions.append(walker.position)
        return positions
 
    #performs branching walk on a per-branch level, not iteratively    
    def get_branching_walk_v1(self, num_steps):
        list_of_branches = []
        self.list_of_walkers = [Particle((0,0), 0)]
        
        while self.list_of_walkers:
            #run for loop over a copy of walkers list
            current_walkers = self.list_of_walkers.copy()
            for w in current_walkers:
                #calculate how long before branch happens
                branch_time = round(self.branch_waiting_dist())
                #perform a shorter walk if we hit the max allowed number of steps
                if w.iteration + branch_time >= num_steps:
                    list_of_branches.append(self.get_walk(w, num_steps - w.iteration))
                    self.list_of_walkers.remove(w)
                    continue
                #otherwise, perform full walk, and then branch
                list_of_branches.append(self.get_walk(w, branch_time))
                #prevent this walker being re-iterated over
                self.list_of_walkers.remove(w)
                angle = self.branch_angle_dist()
                #create two new walkers to iterate over in next for loop
                self.list_of_walkers.append(Particle(w.position, w.angle - angle, w.iteration))
                self.list_of_walkers.append(Particle(w.position, w.angle + angle, w.iteration))
                
              
        return list_of_branches
    
    #performs the same walk as above, but now with an iterative approach
    def get_branching_walk_v2(self, num_steps):
        positions = []
        #dictionary keeps track of all active tips, how long they travel before branching,
        #and the trajectory it has taken so far
        walker_dict = {Particle((0,0), 0): 
                      [round(self.branch_waiting_dist()), [(0,0)] ]}
        
        iteration = 0
        
        while iteration < num_steps:
            #create copy of active walkers, ready to loop over
            current_walkers = walker_dict.copy()
            print("===========================")
            print(f"Current iteration: {iteration}")
            print(f"Current walkers: {current_walkers}")
            print("===========================")
            
            for w in current_walkers:
                #if branching point not yet reached, jump one step
                if w.iteration < walker_dict[w][0]:
                    w.rotate(self.angle_dist())
                    w.move(self.step_dist())
                    w.iteration += 1
                    walker_dict[w][1].append(w.position)
                    continue
                #otherwise, append branch's entire trajectory to list
                positions.append(walker_dict[w][1])
                #remove walker from dictionary since not active
                walker_dict.pop(w)
                angle = self.branch_angle_dist()
                #add two new walkers, these are the child walkers of branching process
                walker_dict[Particle(w.position, w.angle + angle)] = [
                    round(self.branch_waiting_dist()), [w.position] ]
                walker_dict[Particle(w.position, w.angle - angle)] = [
                    round(self.branch_waiting_dist()), [w.position] ]
            
            iteration += 1
        
        #append all 'unfinished' branches to list
        for w in walker_dict:
            positions.append(walker_dict[w][1])
        
        return positions
    
    def graph_walk(self, num_steps, name = None, lw = 0.75):
        branching_walk = self.get_branching_walk_v3(num_steps)
        
        for branch in branching_walk:
            x, y = zip(*branch)
            plt.plot(x, y, color = "black", linewidth = lw)
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
        
        
        

