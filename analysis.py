"""
This file defines two objects, one for the storing and extracting of simulation data, 
the other for the analysis and graphing of said data

pipeline is Simulation -> WalkData -> Analysis

Many samples from a simulation can be run, and to save time, stored for later use in WalkData
then many graphs / analyses made from a single run of the simulation
"""
import matplotlib.pyplot as plt
import numpy as np

class WalkData:
    def __init__(self, data, simulation = None):
        if type(data) == dict:
            self.data = [data]
        else:
            self.data = data
        self.sample_size = len(self.data)
        self.simulation = simulation
        
        self.positions = self.get_positions()
        self.angles = self.get_angles()
        self.somas = self.get_somas()

    def get_positions(self):
        all_pos = []
        for d in self.data:
            walk_pos = [np.array(d[i]["positions"]) for i in d]
            all_pos.append(walk_pos)
        
        return all_pos

    def get_angles(self):
        all_angles = []
        for d in self.data:
            for i in d:
                for theta in d[i]["angles"]:
                    all_angles.append(theta)        
        return all_angles
    
    def get_somas(self):
        all_somas = []
        for d in self.data:
            walk_somas = []
            for i in d:
                if d[i]["soma"]:
                    walk_somas.append(d[i]["positions"][0])
            all_somas.append(np.array(walk_somas))
        
        return all_somas
    
class Analysis:
    def __init__(self, walkdata):
        self.data = walkdata
    
        
    def graph_walk(self, index = 0, lw = 0.75, col = "white", name = None):
        walk = self.data.positions[index]
        somas = self.data.somas[index]
        
        for branch in walk:
            x, y = zip(*branch)
            plt.plot(x, y, color = col, linewidth = lw)
        for soma in somas:
            x, y = zip(soma)
            plt.plot(x, y, color = 'red', marker = 'o', markersize = 2.5)
            
        plt.gca().set_aspect('equal')
        ax = plt.gca()
        ax.set_facecolor("black")
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
        
    def graph_angles(self, axis_min, axis_max, bin_number, name = None):
        bin_width = (axis_max - axis_min) / bin_number
        bins = [(axis_min, axis_min + bin_width)]
        
        for i in range(bin_number - 1):
            bin_start = bins[i][1]
            bin_end = bin_start + bin_width
            bins.append((bin_start, bin_end))
        
        count_dict = {}
        
        for b in bins:
            count = 0
            for angle in self.data.angles:
                if b[0] <= angle < b[1]:
                    count += 1
                count_dict[0.5 * (b[0]+b[1])] = count
        
        plt.plot(count_dict.keys(), count_dict.values())
        plt.xlabel("Walker angle")
        plt.ylabel("Frequency")
        plt.title("Distribution of angles")
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()