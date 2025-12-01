#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Various helper functions for the random walk simulations
"""
import math as maths
import numpy as np
import matplotlib.pyplot as plt

"""
1: Particle class to be used to track walkers in random walk simulations
"""
class Particle:
    def __init__(self, position, angle, iteration = 0):
        self.position = position
        self.angle = angle % (2 * maths.pi)
        self.iteration = iteration
    
    def move(self, distance):
        new_x = distance * maths.cos(self.angle) + self.position[0]
        new_y = distance * maths.sin(self.angle) + self.position[1]
        self.position = (new_x, new_y)
        
    def rotate(self, angle):
        self.angle += angle
        self.angle = self.angle % (2*maths.pi)
    
    def describe(self):
        print(f"Particle at position {self.position} with angle {self.angle} after {self.iteration} iterations")

"""
2: A couple miscelaneous function(s) for graphing and simulating walks
""" 
#runs the function f(*args, **kwargs) n times, collecting the results into a single list
def sample(n, f, *args, **kwargs):
    output = []
    for _ in range(n):
        output.append(f(*args, **kwargs))
    return output

def mean(lst):
    return sum(lst) / len(lst)

#used for seeding initial positions of multi walks, use as list_seq(polygon(n))
def polygon(n, scale = 1):
    return [ (scale * maths.cos(2*k*maths.pi/n), 
              scale * maths.sin(2*k*maths.pi/n) ) for k in range(n)]

#used for seeding initial angles radially for multi walks
def angle(lst):
    return [maths.atan2(pos[1], pos[0]) for pos in lst]

"""
3: Several functions to help mask creation in the annihilating part of simulation
"""
#define a history constant to tweak how far back in time to look:
k = 2
def combine_masks(mask_list):
    return np.bitwise_or.reduce(mask_list)

def mask_past(index_np, iteration_np, walker_index, walker_iteration):
    return (index_np == walker_index) & (iteration_np >= walker_iteration - k)

def mask_relatives(index_np, iteration_np, walker_index, walker_dict):
    parent = walker_dict[walker_index]["parent"]
    w = walker_dict[walker_index]["walker"]
    if parent == None or w.iteration >= k+1:
        return np.zeros_like(index_np, dtype=bool)
    pbt = walker_dict[parent]["branch_time"]
    mask1 = (index_np == parent) & (iteration_np >= pbt-k)
    
    sibling = walker_dict[walker_index]["sibling"]
    mask2 = (index_np == sibling) & (iteration_np <= k)
    
    grandparent = walker_dict[parent]["parent"]
    if (grandparent != None) and (pbt-1 < 0):
        gpbt = walker_dict[grandparent]["branch_time"]
        mask3 = (index_np == grandparent) & (iteration_np >= gpbt-1)
        
        return combine_masks([mask1, mask2, mask3])
    return combine_masks([mask1, mask2])

"""
4: Some functions to be used in global guidance

Instead of phi = K for some constant K, 
suppose instead phi = f(x,y) 
i.e. guidance can vary with position
"""
def constant(c):
    def inner_func(x,y):
        return c 
    return inner_func

def f(x,y):
    x = np.asarray(x)
    y = np.asarray(y)
    return np.where(x<0, 0, np.pi)

#set offset = 0 for radial field, pi/2 for circlular
def swirl(offset):    
    def inner_func(x,y):
        x = np.asarray(x)
        y = np.asarray(y)
        angle = np.arctan2(y,x) + offset
        return np.mod(angle, 2*np.pi)
    return inner_func

"""
5: A few commonly used graph functions
"""
#for testing guidance angle field functions
def graph_angle_field(angle_field, x = (-5,5), y = (-5,5), num_arrows = 10, name = None):
    x_axis = np.linspace(x[0], x[1], num_arrows)
    y_axis = np.linspace(y[0], y[1], num_arrows)
    X, Y = np.meshgrid(x_axis, y_axis)
    
    Angle = angle_field(X,Y)
    
    U = np.cos(Angle)
    V = np.sin(Angle)
    
    plt.quiver(X, Y, U, V)
    if name:
        plt.savefig("plots/" + name + ".png", dpi=300, bbox_inches='tight')
        plt.show()
    else: plt.show()

#graph_angle_field(swirl(0.8))  

def graph_walks(nested_list, somas, lw = 0.75, col = "white", name = None):
    for lst in nested_list:
        x, y = zip(*lst)
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