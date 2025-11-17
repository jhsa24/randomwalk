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
                #remove walker from dictionary since it is no longer active
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
    
    def get_anhilating_walk(self, num_steps, anhilation_radius):
        positions = []
        walker_list = [[Particle((0,0), 0), 
                        round(self.branch_waiting_dist()), 
                        [(0,0)]]]
        iteration = 0
        
        while iteration < num_steps:
            #create copy of active walkers, ready to loop over
            current_walkers = walker_list.copy()
            dead_walkers = []
            print("===========================")
            print(f"Current iteration: {iteration}")
            print(f"Current walkers: {current_walkers}")
            print("===========================")
            
            for w in current_walkers:
                dead = False
                #if walker is too close to a branch, make inactive
                for branch in positions:
                    for duct in branch:
                        delta_x = w[0].position[0] - duct[0]
                        delta_y = w[0].position[1] - duct[1]
                        distance = round((delta_x**2 + delta_y**2) ** 0.5, 9)
                        if 0 < distance < anhilation_radius:
                            w[0].describe()
                            positions.append(w[2])
                            print(f"Walker List: {walker_list}")
                            dead_walkers.append(w)
                            dead = True
                            break
                    if dead:
                        break
                if dead:
                    continue
                
                #if branching point not yet reached, jump one step
                if w[0].iteration < w[1]:
                    w[0].rotate(self.angle_dist())
                    w[0].move(self.step_dist())
                    w[0].iteration += 1
                    w[2].append(w[0].position)
                    continue
                #otherwise, append branch's entire trajectory to list
                positions.append(w[2])
                #remove walker from dictionary since it is no longer active
                dead_walkers.append(w)
                angle = self.branch_angle_dist()
                #add two new walkers, these are the child walkers of branching process
                
                walker_list.append([Particle(w[0].position, w[0].angle + angle),
                                    round(self.branch_waiting_dist()),
                                    [w[0].position]])
                walker_list.append([Particle(w[0].position, w[0].angle - angle),
                                    round(self.branch_waiting_dist()),
                                    [w[0].position]])
            
            for w in dead_walkers:
                if w in walker_list:
                    walker_list.remove(w)
            
            iteration += 1
        
        #append all 'unfinished' branches to list
        for w in walker_list:
            positions.append(w[2])
        
        return positions
    
    def get_barw(self, num_steps, radius = 1):
        walker_dict = {0:{"walker" : Particle((0,0),0),
                          "branch_time" : self.branch_waiting_dist(),
                          "dead" : False,
                          "positions" : [(0,0)]
                          }}
        iteration = 1
        while iteration < num_steps:
            #print(f"Iteration number: {iteration}")
            current_walkers = walker_dict.copy()
            for i in current_walkers:
                w = walker_dict[i]["walker"]
                #skip inactive particles
                if walker_dict[i]["dead"]:
                    continue
                #print(f"Looking at walker with index {i}")
                #w.describe()
                #print("=== Distance Checks ===")
                
                #test neighbor distance
                for j in walker_dict:
                    for pos in walker_dict[j]["positions"]:
                        dx = w.position[0] - pos[0]
                        dy = w.position[1] - pos[1]
                        distance = round((dx**2 + dy**2) ** 0.5, 9)
                        #print(w.position, pos, distance)
                        
                        #if too close, kill the walker and exit loop
                        if 0 < distance < radius:
                            walker_dict[i]["dead"] = True
                            break
                    if walker_dict[i]["dead"]:
                        break
                if walker_dict[i]["dead"]:
                    continue
            
                #if branching doesn't occur, walk one step
                if w.iteration < walker_dict[i]["branch_time"]:
                    w.rotate(self.angle_dist())
                    w.move(self.step_dist())
                    w.iteration += 1
                    walker_dict[i]["positions"].append(w.position)
                    #print(f"Taking one more step to {w.position}")
                    continue
                #otherwise kill the walker and create two new ones
                walker_dict[i]["dead"] = True
                #calculate size of dictionary ready for new indices
                size = len(walker_dict)
                angle = self.branch_angle_dist()
                walker_dict[size] = {"walker" : Particle(w.position, w.angle + angle),
                                     "branch_time" : self.branch_waiting_dist(),
                                     "dead" : False,
                                     "positions" : [w.position]}
                
                walker_dict[size+1] = {"walker" : Particle(w.position, w.angle - angle),
                                     "branch_time" : self.branch_waiting_dist(),
                                     "dead" : False,
                                     "positions" : [w.position]}
                #print("Branching occuring, new dictionary:")
                #print(walker_dict)
            #print("================================")
            iteration += 1
            
        positions = [walker_dict[i]["positions"] for i in walker_dict]
        return positions 
            
    
    def graph_walk(self, num_steps, name = None, lw = 0.75, walk = None):
        if walk:
            branching_walk = walk
        else:    branching_walk = self.get_branching_walk_v2(num_steps)
        
        for branch in branching_walk:
            x, y = zip(*branch)
            plt.plot(x, y, color = "black", linewidth = lw)
        plt.gca().set_aspect('equal')
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
        
        
        

