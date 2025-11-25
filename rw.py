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

from aux import Particle, mean, mask_past, combine_masks
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
        self.guidance_strength = guidance_strength
        self.guidance_angle = guidance_angle
        self.particle = Particle((0,0), 0)
    
    def biased_angle_dist(self, angle):
        angle_difference = angle - self.guidance_angle
        unbiased_angle = self.angle_dist()
        biased_angle = - k * self.guidance_strength * maths.sin(angle_difference)
        return biased_angle + unbiased_angle
    
    #perform walk of given length, output is a list 
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
                #self.particle.angle = self.angle_dist()
                self.particle.move(self.step_dist())
                positions.append(self.particle.position)
            
        return positions

    def get_guided_walk(self, num_steps):
        positions, soma = [], []
        self.particle.position = (0,0)
        positions.append(self.particle.position)
        soma.append(self.particle.position)
        for _ in range(num_steps):
            angle = self.biased_angle_dist(self.particle.angle)
            self.particle.rotate(angle)
            self.particle.move(self.step_dist())
            positions.append(self.particle.position)
        
        return [positions], soma
    
    def get_multi_arw(self, 
                      num_steps, 
                      num_walkers, 
                      radius = 1 ):
        
        walker_dict = {}
        somas = []
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
                              "positions" : [(x,y)]
                              }
            somas.append((x,y))
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
             
                w.rotate(self.angle_dist())
                w.move(self.step_dist())
                w.iteration += 1
                
                walker_i["positions"].append(w.position)
                    
                pos_np = np.vstack([pos_np, w.position])
                index_np = np.hstack([index_np, i])
                iteration_np = np.hstack([iteration_np, w.iteration])
                  
            iteration += 1

        positions = [walker_dict[i]["positions"] for i in walker_dict]
        return positions, somas
    
    def get_guided_multi_arw(self, num_steps, num_walkers, radius):
        walker_dict = {}
        somas = []
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
                              "positions" : [(x,y)]
                              }
            somas.append((x,y))
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
                angle = self.biased_angle_dist(w.angle)
                w.rotate(angle)
                w.move(self.step_dist())
                w.iteration += 1
                
                walker_i["positions"].append(w.position)
                    
                pos_np = np.vstack([pos_np, w.position])
                index_np = np.hstack([index_np, i])
                iteration_np = np.hstack([iteration_np, w.iteration])
                  
            iteration += 1

        positions = [walker_dict[i]["positions"] for i in walker_dict]
        return positions, somas
    
    #hidden helper function for graph_walk(..)
    def __graph_walk_1d(self, num_steps, name = None, lw = 0.5):
        plt.plot(self.get_walk(num_steps), linewidth = lw)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.gca().set_aspect('equal')
        plt.title(f"1d walk of length {num_steps}")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    #hidden helper function for graph_walk(..)
    def __graph_walk_2d(self, num_steps, name = None, lw = 0.5):
        #get_walk outputs tuples, so we need to unpack them to separate x&y
        x, y = zip(*self.get_walk(num_steps))
        plt.plot(x, y, linewidth = lw)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.gca().set_aspect('equal')
        plt.title(f"2d walk of length {num_steps}")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    def graph_walk(self, num_steps, name = None, lw = 0.5):
        if self.dimension == 1:
            self.__graph_walk_1d(num_steps, name, lw)
        if self.dimension == 2:
            self.__graph_walk_2d(num_steps, name, lw)
    
    #useful to plot multiple (1d) walks all together
    def graph_walks(self, num_steps, num_walks, name = None, lw = 0.5):
        if self.dimension == 1:
            for _ in range(num_walks):
                plt.plot(self.get_walk(num_steps), linewidth = lw)
            plt.xlabel("Number of steps")
            plt.ylabel("Position of walker")
            plt.title(f"{num_walks} 1d random walks")
            if name:
                plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
                plt.show()
            else: plt.show()
    
    #plots distribution of points after n steps of a 1d random walk
    def graph_distribution_1d(self, num_steps, num_walks, axis_min, axis_max, bin_number, name = None, lw = 1):
        final_positions = []
        num_iter = int(num_walks / 2)
        #simulate n walks, and append final positions to a single list
        # +1 trick is if discrete steps, we may only hit even/odd numbers only
        for i in range(num_iter):
            final_positions.append(self.get_walk(num_steps)[-1])
            final_positions.append(self.get_walk(num_steps + 1)[-1])
        
        bin_width = (axis_max - axis_min) / bin_number
        bins = [(axis_min, axis_min + bin_width)]
        
        for i in range(bin_number - 1):
            bin_start = bins[i][1]
            bin_end = bin_start + bin_width
            bins.append((bin_start, bin_end))
        
        count_dict = {}
        
        for b in bins:
            count = 0
            for num in final_positions:
                if b[0] <= num < b[1]:
                    count += 1
                count_dict[0.5 * (b[0]+b[1])] = count
                
        plt.plot(count_dict.keys(), count_dict.values(), linewidth = lw)
        plt.xlabel(f"Position of walker after {num_steps} steps")
        plt.ylabel("Frequency")
        plt.title(f"Distribution of positions reached by {num_walks} 1d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    #hidden helper function for graphing mean squared distance
    def __graph_MSD_1d(self, num_steps, num_walks, name = None, lw = 1):
        list_of_walks, list_of_averages = [], []
        for _ in range(num_walks):
            #simulate n walks, put all positions into a single nested list
            list_of_walks.append(self.get_walk(num_steps))
        
        #iterate for each time step
        for t in range(num_steps + 1):
            #calculate distance squared
            distance_squared = [abs(list_of_walks[i][t])**2 for i in range(num_walks)]
            #calculate average (using function from aux.py)
            list_of_averages.append(mean(distance_squared))
        
        plt.plot(list_of_averages, linewidth = lw)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Squared Distance")
        plt.title(f"MSD of {num_walks} 1d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    #hidden helper function for graphing mean squared distance
    def __graph_MSD_2d(self, num_steps, num_walks, name = None, lw = 1):
        list_of_walks, list_of_averages = [], []
        #simulate n walks, put all positions into a single nested list
        for _ in range(num_walks):
            list_of_walks.append(self.get_walk(num_steps))
        
        #iterate for each time step
        for t in range(num_steps + 1):
            #get list of positions at time t from all simulated walks
            positions = [list_of_walks[i][t] for i in range(num_walks)]
            #calculate distances, and average
            distance_squared = [pos[0]**2 + pos[1]**2 for pos in positions]
            list_of_averages.append(mean(distance_squared))
        
        plt.plot(list_of_averages, linewidth = lw)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Squared Distance")
        plt.title(f"MSD of {num_walks} 2d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    def graph_MSD(self, num_steps, num_walks, name = None, lw = 1):
        if self.dimension == 1:
            self.__graph_MSD_1d(num_steps, num_walks, name, lw)
        if self.dimension == 2:
            self.__graph_MSD_2d(num_steps, num_walks, name, lw)