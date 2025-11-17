"""
Code for testing discrete random walks, mainly in 2 dimensions

Idea is to use a particle class to keep track of position and rotation
We can then move the particle using its own internal commands,
and append its position to a list each iteration, and graph the list
"""

import random
import math
import matplotlib.pyplot as plt

from aux import mean

#defines a Particle class that can be moved and rotated ready for
#use in random walk functions
class Particle:
    def __init__(self, position, angle, iteration = 0):
        self.position = position
        self.angle = angle % (2 * math.pi)
        self.iteration = iteration
    
    def move(self, distance):
        new_x = distance * math.cos(self.angle) + self.position[0]
        new_y = distance * math.sin(self.angle) + self.position[1]
        self.position = (new_x, new_y)
        
    def rotate(self, angle):
        self.angle += angle % (2 * math.pi)
    
    def describe(self):
        print(f"Particle at position {self.position} with angle {self.angle} after {self.iteration} iterations")
        
class RandomWalk:
    def __init__(self, 
                 dimension = 2, 
                 step_dist = lambda: 1, 
                 angle_dist = lambda: random.choice([math.pi * i/2 for i in range(4)])):
        
        self.dimension = dimension
        self.step_dist = step_dist
        self.angle_dist = angle_dist
        self.particle = Particle((0,0), 0)
        
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
    
    #hidden helper function for graph_walk(..)
    def __graph_walk_1d(self, num_steps, name = None, lw = 0.5):
        plt.plot(self.get_walk(num_steps), linewidth = lw)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
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
    
    #old code for mean distance graphs, probably don't need it anymore            
    """def __graph_MD_1d(self, num_steps, num_walks, name = None, lw = 1):
        list_of_walks, list_of_averages = [], []
        for _ in range(num_walks):
            list_of_walks.append(self.get_walk(num_steps))
        
        for t in range(num_steps + 1):
            positions = [list_of_walks[i][t] for i in range(num_walks)]
            list_of_averages.append(mean(positions))
        
        plt.plot(list_of_averages, linewidth = lw)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Distance")
        plt.title(f"MD of {num_walks} 1d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
        
    def __graph_MD_2d(self, num_steps, num_walks, name = None, lw = 1):
        list_of_walks, list_of_averages = [], []
        for _ in range(num_walks):
            list_of_walks.append(self.get_walk(num_steps))
        
        for t in range(num_steps + 1):
            positions = [list_of_walks[i][t] for i in range(num_walks)]
            distances = [(pos[0]**2 + pos[1]**2) ** 1/2 for pos in positions]
            list_of_averages.append(mean(distances))
        
        plt.plot(list_of_averages, linewidth = lw)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Distance")
        plt.title(f"MD of {num_walks} 2d random walks")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    def graph_MD(self, num_steps, num_walks, name = None, lw = 1):
        if self.dimension == 1:
            self.__graph_MD_1d(num_steps, num_walks, name, lw)
        if self.dimension == 2:
            self.__graph_MD_2d(num_steps, num_walks, name, lw)"""