"""
Code for Branching, Annihilating, Arresting, Interacting etc.... Random Walks
"""

import math as maths
import matplotlib.pyplot as plt
import numpy as np

from aux import Particle, combine_masks, mask_past, mask_relatives
from distributions import uniform, exponential


class BranchingRandomWalk:
    def __init__(self, 
                 dimension = 2, 
                 step_dist = lambda : 1, 
                 angle_dist = uniform(-maths.pi/5, maths.pi/5), 
                 branch_waiting_dist = exponential(1/15),
                 branch_angle_dist = lambda : maths.pi / 3):
                
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
        self.branch_waiting_dist = branch_waiting_dist
        self.branch_angle_dist = branch_angle_dist
        self.list_of_walkers = []
        
  
    #performs a branching random walk, outputting a nested list [ [(x1,y1), (x2,y2),  ], [...], ...]
    #each sublist contains the trajectory of a single branch
    def get_brw(self, num_steps):
        #dictionary keeps track of all walkers, how long they travel before branching,
        #whether they are active and the trajectory taken so far
        walker_dict = {0:{"walker" : Particle((0,0),0),
                          "branch_time" : self.branch_waiting_dist(),
                          "dead" : False,
                          "positions" : [(0,0)],
                          "parent" : None,
                          "sibling" : None
                          }}

        iteration = 1
        
        while iteration < num_steps:
            #initialise list of walkers that branch this iteration
            new_branches = []
            #loop over ids in walker_dict
            for i in walker_dict.keys():
                #pre-calculate walker's sub-dictionary and its particle
                walker_i = walker_dict[i]
                w = walker_i["walker"]
                #skip inactive walkers
                if walker_i["dead"]:
                    continue
                
                #if branching point not yet reached, jump one step
                if w.iteration < walker_i["branch_time"]:
                    w.rotate(self.angle_dist())
                    w.move(self.step_dist())
                    w.iteration += 1
                    walker_i["positions"].append(w.position)
                    continue
                
                #otherwise, kill walker and add index to list of new branches
                walker_i["dead"] = True
                new_branches.append(i)
                
            #add new branches for this iteration after for loop
            for parent_id in new_branches:
                #calculate size of dictionary ready for new indices
                size = len(walker_dict)
                angle = self.branch_angle_dist()
                parent = walker_dict[parent_id]["walker"]
                walker_dict[size] = {"walker" : Particle(parent.position, parent.angle + angle),
                                     "branch_time" : self.branch_waiting_dist(),
                                     "dead" : False,
                                     "positions" : [parent.position],
                                     "parent" : parent_id,
                                     "sibling" : size + 1}
            
                walker_dict[size+1] = {"walker" : Particle(parent.position, parent.angle - angle),
                                     "branch_time" : self.branch_waiting_dist(),
                                     "dead" : False,
                                     "positions" : [parent.position],
                                     "parent" : parent_id,
                                     "sibling" : size}        
            iteration += 1
        
        positions = [walker_dict[i]["positions"] for i in walker_dict]
        return positions
    
   #performs a branching annihilating random walk, outputting a nested list [ [(x1,y1), (x2,y2),  ], [...], ...]
   #each sublist contains the trajectory of a single branch
   #if a particle gets too close to an existing duct, it stops all movement (annihilation)
    def get_barw(self, num_steps, radius = 1):
        walker_dict = {0:{"walker" : Particle((0,0),0),
                          "branch_time" : self.branch_waiting_dist(),
                          "dead" : False,
                          "positions" : [(0,0)],
                          "parent" : None,
                          "sibling" : None
                          }}
        
        radius_squared = radius*radius
        iteration = 1
        
        pos_np = np.array([(0,0)])
        index_np = np.array([0])
        iteration_np = np.array([0])
        
        while iteration < num_steps:
            new_branches = []
            #print(f"Iteration number: {iteration}")
            for i in list(walker_dict.keys()):
                walker_i = walker_dict[i]
                w = walker_i["walker"]
                #skip inactive particles
                if walker_i["dead"]:
                    continue
    
                #calculate distance of all points from walker position
                diff = pos_np - np.array(w.position)
                distance_squared = np.einsum('ij,ij->i', diff, diff)
                #distance_squared = (diff * diff).sum(axis=1)

                #create masks to ignore certain positions
                #masks return true if they should be ignored
                m_self = distance_squared < 1e-9
                m_past = mask_past(index_np, iteration_np, i, w.iteration)
                m_rel = mask_relatives(index_np, iteration_np, i, walker_dict)
                
                
                #mask = mask_self | mask_past | mask_parent | mask_sibling
                mask = combine_masks([m_self, m_past, m_rel])
                too_close = (distance_squared < radius_squared) & ~mask
                
                if np.any(too_close):
                    walker_i["dead"] = True
                    continue
                
                #if branching doesn't occur, walk one step
                if w.iteration < walker_i["branch_time"]:
                    w.rotate(self.angle_dist())
                    w.move(self.step_dist())
                    w.iteration += 1
                    
                    walker_i["positions"].append(w.position)
                    
                    pos_np = np.vstack([pos_np, w.position])
                    index_np = np.append(index_np, i)
                    iteration_np = np.append(iteration_np, w.iteration)
                    continue
                
                #otherwise kill the walker and restart loop
                walker_i["dead"] = True
                new_branches.append(i)
            
            #after one full loop, create all new branched particles together
            for parent_id in new_branches:
                #calculate size of dictionary ready for new indices
                size = len(walker_dict)
                angle = self.branch_angle_dist()
                parent = walker_dict[parent_id]["walker"]
                
                walker_dict[size] = {"walker" : Particle(parent.position, parent.angle + angle),
                                     "branch_time" : self.branch_waiting_dist(),
                                     "dead" : False,
                                     "positions" : [parent.position],
                                     "parent" : parent_id,
                                     "sibling" : size + 1}
            
                walker_dict[size+1] = {"walker" : Particle(parent.position, parent.angle - angle),
                                     "branch_time" : self.branch_waiting_dist(),
                                     "dead" : False,
                                     "positions" : [parent.position],
                                     "parent" : parent_id,
                                     "sibling" : size}
                
            iteration += 1

        positions = [walker_dict[i]["positions"] for i in walker_dict]
        return positions