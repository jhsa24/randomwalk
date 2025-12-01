"""
Code for Branching, Annihilating, Arresting, Interacting etc.... Random Walks
"""

import math as maths
import numpy as np

from aux import Particle, combine_masks, mask_past, mask_relatives, constant
from distributions import uniform, exponential

"""
some constants:
"""
k = 0.1 #global guidance strength scaler

class BranchingRandomWalk:
    def __init__(self, 
                 dimension = 2,
                 step_dist = lambda : 1, 
                 angle_dist = uniform(-maths.pi/5, maths.pi/5), 
                 branch_waiting_dist = exponential(1/15),
                 branch_angle_dist = lambda : maths.pi / 3,
                 initial_pos_dist = uniform(-5,5),
                 initial_angle_dist = uniform(0, 2*maths.pi),
                 guidance_strength = 0,
                 guidance_angle = 0):
                
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
        self.branch_waiting_dist = branch_waiting_dist
        self.branch_angle_dist = branch_angle_dist
        self.initial_pos_dist = initial_pos_dist
        self.initial_angle_dist = initial_angle_dist
        
        if type(guidance_angle) in [float, int]:
            self.guidance_angle = constant(guidance_angle)
        else:
            self.guidance_angle = guidance_angle
        if type(guidance_strength) in [float, int]:
            self.guidance_strength = constant(guidance_strength)
        else:
            self.guidance_strength = guidance_strength
          
    #defines a bias function to nudge walker in the direction of the guiding field
    #this is used in place of the angle distribution on its own    
    def biased_angle_dist(self, walker):
        angle_difference = walker.angle - self.guidance_angle(*walker.position)
        unbiased_angle = self.angle_dist()
        biased_angle = - k * self.guidance_strength(*walker.position) * maths.sin(angle_difference)
        return biased_angle + unbiased_angle
    
    
    #performs a branching random walk, keeping track of all information with a dictionary
    #returned output is this dictionary, to be stored in WalkData and analysed in Analysis objects
    def get_brw(self, num_steps, initial_pos = (0,0), initial_angle = 0):
        #dictionary keeps track of all walkers, how long they travel before branching,
        #whether they are active and the trajectory taken so far
        walker_dict = {0:{"walker" : Particle(initial_pos, initial_angle),
                          "soma" : True,
                          "dead" : False,
                          "positions" : [initial_pos],
                          "angles" : [initial_angle],
                          "branch_time" : maths.ceil(self.branch_waiting_dist()),
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
                if w.iteration <= walker_i["branch_time"]:
                    angle = self.biased_angle_dist(w)
                    w.rotate(angle)
                    w.move(self.step_dist())
                    w.iteration += 1
                    walker_i["positions"].append(w.position)
                    walker_i["angles"].append(w.angle)
                    continue
                
                #otherwise, kill walker and add index to list of new branches
                walker_i["dead"] = True
                new_branches.append(i)
                
            #add new branches for this iteration after for loop
            for parent_id in new_branches:
                #calculate size of dictionary ready for new indices
                size = len(walker_dict)
                angle1 = abs(self.branch_angle_dist())
                angle2 = abs(self.branch_angle_dist())
                parent = walker_dict[parent_id]["walker"]
                walker_dict[size] = {"walker" : Particle(parent.position, parent.angle + angle1),
                                     "soma" : False,
                                     "dead" : False,
                                     "positions" : [parent.position],
                                     "angles" : [parent.angle + angle1],
                                     "branch_time" : maths.ceil(self.branch_waiting_dist()),
                                     "parent" : parent_id,
                                     "sibling" : size + 1}
            
                walker_dict[size+1] = {"walker" : Particle(parent.position, parent.angle - angle2),
                                       "soma" : False,
                                       "dead" : False,
                                       "positions" : [parent.position],
                                       "angles" : [parent.angle - angle2],
                                       "branch_time" : maths.ceil(self.branch_waiting_dist()),
                                       "parent" : parent_id,
                                       "sibling" : size}        
            iteration += 1
                
        return walker_dict
    
    
    #performs a branching annihilating random walk, returning the dictionary
    #if a particle gets too close to an existing duct, it stops all movement (annihilation)
    #walker movement follows global guidance
    def get_barw(self, num_steps, radius = 1, initial_pos = (0,0), initial_angle = 0):
        walker_dict = {0:{"walker" : Particle(initial_pos, initial_angle),
                          "soma" : True,
                          "dead" : False,
                          "positions" : [initial_pos],
                          "angles" : [initial_angle],
                          "branch_time" : maths.ceil(self.branch_waiting_dist()),
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
                    angle = self.biased_angle_dist(w)
                    w.rotate(angle)
                    w.move(self.step_dist())
                    w.iteration += 1
                    
                    walker_i["positions"].append(w.position)
                    walker_i["angles"].append(w.angle)
                    
                    pos_np = np.vstack([pos_np, w.position])
                    index_np = np.hstack([index_np, i])
                    iteration_np = np.hstack([iteration_np, w.iteration])
                    continue
                
                #otherwise kill the walker and restart loop
                walker_i["dead"] = True
                new_branches.append(i)
            
            #after one full loop, create all new branched particles together
            for parent_id in new_branches:
                #calculate size of dictionary ready for new indices
                size = len(walker_dict)
                angle1 = abs(self.branch_angle_dist())
                angle2 = abs(self.branch_angle_dist())
                parent = walker_dict[parent_id]["walker"]
                
                walker_dict[size] = {"walker" : Particle(parent.position, parent.angle + angle1),
                                     "soma" : False,
                                     "dead" : False,
                                     "positions" : [parent.position],
                                     "angles" : [parent.angle + angle1],
                                     "branch_time" : maths.ceil(self.branch_waiting_dist()),
                                     "parent" : parent_id,
                                     "sibling" : size + 1}
            
                walker_dict[size+1] = {"walker" : Particle(parent.position, parent.angle - angle2),
                                       "soma" : False,
                                       "dead" : False,
                                       "positions" : [parent.position],
                                       "angles" : [parent.angle - angle2],
                                       "branch_time" : maths.ceil(self.branch_waiting_dist()),
                                       "parent" : parent_id,
                                       "sibling" : size}
                
            iteration += 1

        return walker_dict
    
    def get_multi_barw(self, num_steps, num_walkers, radius = 1):
        walker_dict = {}

        pos_np = np.zeros((num_walkers, 2))
        index_np = np.zeros(num_walkers, dtype=int)
        iteration_np = np.zeros(num_walkers, dtype=int)
        
        for i in range(num_walkers):
            theta = self.initial_angle_dist()
            x = self.initial_pos_dist()
            if type(x) == tuple:
                x,y = x
            else:
                y = self.initial_pos_dist()
            
            walker_dict[i] = {"walker" : Particle((x,y), theta),
                              "soma" : True,
                              "dead" : False,  
                              "positions" : [(x,y)],
                              "angles" : [theta],
                              "branch_time" : maths.ceil(self.branch_waiting_dist()),
                              "parent" : None,
                              "sibling" : None}
            pos_np[i] = (x,y)
            index_np[i] = i
            iteration_np[i] = 0
            
        radius_squared = radius*radius
        iteration = 1
             
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
                if w.iteration <= walker_i["branch_time"]:
                    angle = self.biased_angle_dist(w)
                    w.rotate(angle)
                    w.move(self.step_dist())
                    w.iteration += 1
                    
                    walker_i["positions"].append(w.position)
                    walker_i["angles"].append(w.angle)
                    
                    pos_np = np.vstack([pos_np, w.position])
                    index_np = np.hstack([index_np, i])
                    iteration_np = np.hstack([iteration_np, w.iteration])
                    continue
                
                #otherwise kill the walker and restart loop
                walker_i["dead"] = True
                new_branches.append(i)
            
            #after one full loop, create all new branched particles together
            for parent_id in new_branches:
                #calculate size of dictionary ready for new indices
                size = len(walker_dict)
                angle1 = abs(self.branch_angle_dist())
                angle2 = abs(self.branch_angle_dist())
                parent = walker_dict[parent_id]["walker"]
                
                walker_dict[size] = {"walker" : Particle(parent.position, parent.angle + angle1),
                                     "soma" : False,
                                     "dead" : False, 
                                     "positions" : [parent.position],
                                     "angles" : [parent.angle + angle1],
                                     "branch_time" : maths.ceil(self.branch_waiting_dist()),
                                     "parent" : parent_id,
                                     "sibling" : size + 1}
            
                walker_dict[size+1] = {"walker" : Particle(parent.position, parent.angle - angle2),
                                     "soma" : False,
                                     "dead" : False, 
                                     "positions" : [parent.position],
                                     "angles" : [parent.angle - angle2],
                                     "branch_time" : maths.ceil(self.branch_waiting_dist()),
                                     "parent" : parent_id,
                                     "sibling" : size}
                
            iteration += 1

        return walker_dict