"""
Code for discrete random walks, in 1 and 2 dimensions

Originally this was used as foundations for more complicated branching algorithms
Should still be a useful testing ground to implement preliminary versions of
self avoidance, global guidance, and other effects.

Idea is to use the particle class in aux.py to track position and rotation
We can then move the particle using its own internal commands,
and append its position to a list each iteration, and graph the list
"""

import math as maths
import matplotlib.pyplot as plt
import numpy as np

from aux import Particle, mask_past, combine_masks, constant
from distributions import pmf, uniform

"""
some constants:
"""
k = 0.1 #global guidance strength scaler
     
class RandomWalk:
    def __init__(self, 
                 dimension = 2, 
                 step_dist = lambda: 1, 
                 angle_dist = pmf({maths.pi * i/2:1 for i in range(4)}),
                 initial_pos_dist = uniform(-5,5),
                 initial_angle_dist = uniform(0, 2*maths.pi),
                 guidance_strength = 0,
                 guidance_angle = 0):
        
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
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
            
        self.particle = Particle((0,0), 0)
    
    
    def biased_angle_dist(self, walker):
        angle_difference = walker.angle - self.guidance_angle(*walker.position)
        unbiased_angle = self.angle_dist()
        biased_angle = - k * self.guidance_strength(*walker.position) * maths.sin(angle_difference)
        return biased_angle + unbiased_angle


    def get_rw(self, num_steps, initial_pos = (0,0), initial_angle = 0):
        walker_dict = {0: {
                        "walker" : Particle(initial_pos, initial_angle),
                        "dead" : False,
                        "positions" : [initial_pos],
                        "angles" : [initial_angle],
                        "parent" : None
                        }}
        w = walker_dict[0]["walker"]
        
        for _ in range(num_steps):
            angle = self.biased_angle_dist(w)
            w.rotate(angle)
            w.move(self.step_dist())
            walker_dict[0]["positions"].append(w.position)
            walker_dict[0]["angles"].append(w.angle)
                    
        return walker_dict
    
    
    def get_multi_arw(self, num_steps, num_walkers, radius):
        walker_dict = {}
        pos_np = np.zeros((num_walkers, 2))
        index_np = np.zeros(num_walkers, dtype=int)
        iteration_np = np.zeros(num_walkers, dtype=int)
        
        for i in range(num_walkers):
            x = self.initial_pos_dist()
            if type(x) == tuple:
                x,y = x
            else:
                y = self.initial_pos_dist()
            
            theta = self.initial_angle_dist()
            walker_dict[i] = {"walker" : Particle((x,y),theta),
                              "dead" : False,
                              "positions" : [(x,y)],
                              "angles" : [theta],
                              "parent" : None
                              }
            pos_np[i] = (x,y)
            index_np[i] = i
            iteration_np[i] = 0
            
        radius_squared = radius*radius
        iteration = 1
       
        while iteration < num_steps:
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
                mask = m_self | m_past
                #mask = combine_masks([m_self, m_past])
                too_close = (distance_squared < radius_squared) & ~mask
               
                if np.any(too_close):
                    walker_i["dead"] = True
                    continue
                angle = self.biased_angle_dist(w)
                w.rotate(angle)
                w.move(self.step_dist())
                w.iteration += 1
                
                walker_i["positions"].append(w.position)
                walker_i["angles"].append(w.angle)
                    
                pos_np = np.vstack([pos_np, w.position])
                index_np = np.hstack([index_np, i])
                iteration_np = np.hstack([iteration_np, w.iteration])
                  
            iteration += 1

        return walker_dict