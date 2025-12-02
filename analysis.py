"""
This file defines two objects, one for the storing and extracting of simulation data, 
the other for the analysis and graphing of said data

pipeline is Simulation -> WalkData -> Analysis

Many samples from a simulation can be run, and to save time, stored for later use in WalkData
then many graphs / analyses made from a single run of the simulation
"""
import matplotlib.pyplot as plt
import numpy as np

from aux import mean

class WalkData:
    def __init__(self, data, simulation = None):
        if type(data) == dict:
            self.data = [data]
        else:
            self.data = data
        self.sample_size = len(self.data)
        self.simulation = simulation
        
        self.somas = self.get_somas()
        self.numpyify()
    
    def get_branch_start_times(self):
        all_times = []
        for walk in self.data:
            start_times = {}
            incomplete = True
            
            while incomplete:
                incomplete = False
                for index in walk:
                    if index in start_times:
                        continue
                    parent = walk[index]["parent"]
                    if parent is None:
                        start_times[index] = 0
                    elif parent in start_times:
                        start_times[index] = start_times[parent] + walk[parent]["branch_time"]
                    else:
                        incomplete = True
            
            all_times.append(start_times)
        return all_times
    
    def get_positions_by_iter(self):
        all_pos = []
        branch_start_times = self.get_branch_start_times()
        
        for i, wd in enumerate(self.data):
            positions = {}
            for index in wd:
                iteration = branch_start_times[i][index]
                for pos in wd[index]["positions"]:
                    positions.setdefault(iteration, []).append(pos)
                    iteration += 1
            all_pos.append(positions)
        return all_pos
    
    
    def numpyify(self):
        all_pos = []
        all_iters = []
        all_angles = []
        all_branch_id = []
        all_parent_id = []
        all_sample_id = []
        
        branch_start_times = self.get_branch_start_times()
        for i, wd in enumerate(self.data):
            for branch in wd:
                positions = np.array(wd[branch]["positions"])
                angles = np.array(wd[branch]["angles"])
                
                start = branch_start_times[i][branch]
                length = len(positions)
                iterations = np.arange(start, start + length)
                
                all_pos.append(positions)
                all_iters.append(iterations)
                all_angles.append(angles)
                all_branch_id.append(np.full(length, branch))
                all_parent_id.append(np.full(length, wd[branch]["parent"]))
                all_sample_id.append(np.full(length, i))
        
        self.positions = np.vstack(all_pos)
        self.iters = np.concatenate(all_iters)
        self.angles = np.concatenate(all_angles)
        self.branch_id = np.concatenate(all_branch_id)
        self.parent_id = np.concatenate(all_parent_id)
        self.sample_id = np.concatenate(all_sample_id)
            
    
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
        
        self.iters = walkdata.iters
        self.angles = walkdata.angles
        self.positions = walkdata.positions
        self.branch_id = walkdata.branch_id
        self.parent_id = walkdata.parent_id
        self.sample_id = walkdata.sample_id
        
    def tmp(self, index=0):
        pass
    
    
    def graph_walk(self, index = 0, lw = 0.75, col = "white", name = None):
        #create arrays for masking sample
        m_sample = self.sample_id == index
        
        somas = self.positions[(self.iters == 0) & m_sample]
        unique_branches = np.unique(self.branch_id[m_sample])
        
        for b in unique_branches:
            pts = self.positions[(self.branch_id == b) & m_sample]
            x = [p[0] for p in pts]
            y = [p[1] for p in pts]
            plt.plot(x, y, color = col, linewidth = lw)
        
        for soma in somas:
            x, y = soma[0], soma[1]
            plt.plot(x, y, color = 'red', marker = 'o', markersize = 2.5)
            
        plt.gca().set_aspect('equal')
        ax = plt.gca()
        ax.set_facecolor("black")
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
    
    
    def graph_MSD(self, name = None):
        
        unique_iter = np.unique(self.iters)
        list_of_averages = []
        
        for t in range(len(unique_iter)):
            pos = self.positions[self.iters == t]
            distance_squared = (pos*pos).sum(axis=1)
            
            list_of_averages.append(mean(distance_squared))
        
        plt.plot(list_of_averages)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Squared Distance")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()
            
    
        
    def graph_angles(self, indices = [], bin_number = None, name = None):
        if type(indices) in [float, int]:
            indices = [indices]
        
        angles = self.angles
        plt.hist(angles, bins = bin_number)
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
            plt.show()
        else: plt.show()