"""
This file defines two objects, one for the storing and extracting of simulation data, 
the other for the analysis and graphing of said data

pipeline is Simulation -> WalkData -> Analysis

Many samples from a simulation can be run, and to save time, stored for later use in WalkData
then many graphs / analyses made from a single run of the simulation
"""
import matplotlib.pyplot as plt
from  matplotlib.animation import FuncAnimation
import numpy as np
import pickle

class WalkData:
    def __init__(self, data, simulation = None):
        if type(data) == dict:
            self.data = [data]
        else:
            self.data = data
        self.sample_size = len(self.data)
        
        self.simulation = simulation
        self.metadata = {}
        if simulation:
            self.metadata["dimension"] = getattr(simulation, "dimension", None)
            self.metadata["step_dist"] = getattr(simulation, "step_dist", None)
            self.metadata["angle_dist"] = getattr(simulation, "angle_dist", None)
            self.metadata["branch_waiting_dist"] = getattr(simulation, "branch_waiting_dist", None)
            self.metadata["branch_angle_dist"] = getattr(simulation, "branch_angle_dist", None)
            self.metadata["initial_pos_dist"] = getattr(simulation, "initial_pos_dist", None)
            self.metadata["initial_angle_dist"] = getattr(simulation, "initial_angle_dist", None)
            self.metadata["guidance_angle"] = getattr(simulation, "guidance_angle", None)
            self.metadata["guidance_strength"] = getattr(simulation, "guidance_strength", None)
            
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
                        #start_times[index] = start_times[parent] + walk[parent]["branch_time"]
                        start_times[index] = start_times[parent] + len(walk[parent]["positions"]) - 1
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
            
    
    def save(self, filename):
        if not filename.endswith(".pkl"):
            filename += ".pkl"
        
        with open("sim/" + filename, 'wb') as file:
            pickle.dump(self, file)
        
    
    
class Analysis:
    def __init__(self, walkdata):
        
        if isinstance(walkdata, str):
            walkdata = self.load(walkdata)
                
        self.data = walkdata.data
        
        self.iters = walkdata.iters
        self.angles = walkdata.angles
        self.positions = walkdata.positions
        self.branch_id = walkdata.branch_id
        self.parent_id = walkdata.parent_id
        self.sample_id = walkdata.sample_id
        
  
    def load(self, filename):
        if not filename.endswith(".pkl"):
            filename += ".pkl"
            
        with open("sim/" + filename, 'rb') as file:
            return pickle.load(file)
    
    
    def graph_walk(self, 
                   sample = 0, 
                   iteration = None, 
                   lw = 0.75, 
                   col = "black", 
                   name = None,
                   show = True):
        
        if iteration is None:
            iteration = [0, float("inf")]
        elif isinstance(iteration, (int, float)):
            iteration = [0, iteration]
        elif len(iteration != 2):
            raise ValueError("Iteration must be None, Int, Float, or a [low, high] list")
            
        #create arrays for masking sample and iterations
        m_sample = self.sample_id == sample
        m_iter = (iteration[0] <= self.iters) & (self.iters <= iteration[1])
        
        mask = (m_sample) & (m_iter)
        branch_ids = self.branch_id[mask]
        positions = self.positions[mask]
        
        unique_branches = np.unique(branch_ids)
        
        for b in unique_branches:
            pts = positions[branch_ids == b]
            x, y = pts[:,0], pts[:,1] 
            plt.plot(x, y, color = col, linewidth = lw)
        
        
        if iteration[0] == 0:
            somas = self.positions[(self.iters == 0) & m_sample]
            for soma in somas:
                x, y = soma[0], soma[1]
                plt.plot(x, y, color = 'red', marker = 'o', markersize = 2.5)
            
        plt.gca().set_aspect('equal')
        ax = plt.gca()
        ax.set_facecolor("white")
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        if show: 
            plt.show()


    def graph_walkers_at_t(self,
                           sample = 0,
                           iteration = 0,
                           col = 'red',
                           name = None,
                           show = True):
        
        m_sample = self.sample_id == sample
        m_iter = self.iters == iteration
        
        current_pts = self.positions[m_sample & m_iter]
        x, y = current_pts[:,0], current_pts[:,1]
        plt.scatter(x, y, color = col)
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        if show: 
            plt.show()
    
            
    def graph_MSD(self, 
                  sample = None, 
                  name = None,
                  show = True,
                  x_log = False,
                  y_log = False):
        
        if sample is None:
            sample = np.unique(self.sample_id)
        elif isinstance(sample, (int, float)):
            sample = [sample]
        m_sample = np.isin(self.sample_id, sample)
        
        iters = self.iters[m_sample]
        positions = self.positions[m_sample]
        distance_squared = (positions * positions).sum(axis=1)
        
        unique_iter = np.unique(iters)
        list_of_averages = []

        for t in unique_iter:
            list_of_averages.append(np.mean(distance_squared[iters == t]))
        
        plt.plot(list_of_averages)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean Squared Distance")
        if x_log:
            plt.xscale("log")
        if y_log:
            plt.yscale("log")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        if show: 
            plt.show()
    
    
    def graph_num_walkers(self,
                             sample = None,
                             name = None,
                             show = True):
        
        if sample == None:
            sample = np.unique(self.sample_id)
        elif isinstance(sample, (float, int)):
            sample = [sample]
        
        m_sample = np.isin(self.sample_id, sample)
        
        iters = self.iters[m_sample]
        positions = self.positions[m_sample]
        samples = self.sample_id[m_sample]
            
        unique_iter = np.unique(iters)
        list_of_averages = []
        
        for t in unique_iter:
            m_iter = iters == t
            num_walkers = np.unique(positions[m_iter], axis = 0)
            num_samples = np.unique(samples[m_iter])
            list_of_averages.append(len(num_walkers)/len(num_samples))
            
        plt.plot(list_of_averages)
        plt.xlabel("Number of steps")
        plt.ylabel("Mean number of active walkers")
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        if show: 
            plt.show()
    
        
    def graph_angles(self, 
                     sample = None, 
                     bin_number = None, 
                     name = None, 
                     show = True):
        
        if sample == None:
            sample = np.unique(self.sample_id)
        elif isinstance(sample, (float, int)):
            sample = [sample]
        
        m_sample = np.isin(self.sample_id, sample)
            
        angles = self.angles[m_sample]
        plt.hist(angles, bins = bin_number, density = True)
        
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        if show: 
            plt.show()
            

    def graph_branch_lengths(self,
                             sample = None,
                             bin_number = None,
                             name = None,
                             show = True,
                             x_log = False,
                             y_log = False):
   
        if sample == None:
            sample = np.unique(self.sample_id)
        elif isinstance(sample, (float, int)):
            sample = [sample]

        m_sample = np.isin(self.sample_id, sample)
        sample_id = self.sample_id[m_sample]
        branch_id = self.branch_id[m_sample]
        
        branch_lengths = []
        for i in sample:
            current_branches = branch_id[sample_id == i]
            _, count = np.unique(current_branches, return_counts=True)
            branch_lengths.extend(count)
        
        plt.hist(branch_lengths, bins = bin_number)
        
        if x_log:
            plt.xscale("log")
        if y_log:
            plt.yscale("log")
        if name:
            plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        if show: 
            plt.show()
       
        
    def animate_walk(self, 
                     sample = 0, 
                     name = "Animation", 
                     show_dot = True,
                     show_line = True,
                     lw = 0.75,
                     m_size = 50,
                     col = 'black',
                     bg = 'white',
                     m_col = 'red'):
        
        if not name.endswith(".mp4"):
            name += ".mp4"
        buffer = 0.08
        m_sample = self.sample_id == sample
        
        positions = self.positions[m_sample]
        x_min, x_max = positions[:,0].min(), positions[:,0].max()
        y_min, y_max = positions[:,1].min(), positions[:,1].max()

        fig, ax = plt.subplots(figsize=(6,6))
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
        ax.set_aspect('equal')
        ax.set_facecolor(bg)
        
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        """
        ax.set_xlim(x_min - buffer * abs(x_min), x_max + buffer * abs(x_max))
        ax.set_ylim(y_min - buffer * abs(y_min), y_max + buffer * abs(y_max))
        """
        ax.set_xlim(min(x_min, y_min), max(x_max, y_max))
        ax.set_ylim(min(x_min, y_min), max(x_max, y_max))
        branches = np.unique(self.branch_id[self.sample_id == sample])
        lines = []

        for b in branches:
            line, = ax.plot([], [], color = col, linewidth = lw)
            lines.append(line)
        active_walkers_scatter = ax.scatter([], [], color = m_col, marker='o', s = m_size)
        

        def update(frame):
            m_iter = self.iters <= frame
            
            for line, b in zip(lines, branches):
                pts = self.positions[(self.branch_id == b) & (m_sample) & (m_iter)]
                
                if show_line:
                    if len(pts) > 0:
                        line.set_data(pts[:,0], pts[:,1])
                    else:
                        line.set_data([], [])
                
                if show_dot:
                    active_walkers = self.positions[(self.sample_id == sample) & (self.iters == frame)]
                    active_walkers_scatter.set_offsets(active_walkers)
            
            return lines + [active_walkers_scatter]

        ani = FuncAnimation(fig, 
                            update, 
                            frames = range(0, self.iters[self.sample_id == sample].max() + 30), 
                            blit=True
                            )
  
        ani.save("animations/" + name, fps = 30, dpi = 300)