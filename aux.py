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
    #old definition of branch_length, doesn't work if p_b is a float
    #p_len = walker_dict[parent]["branch_time"]
    p_len = len(walker_dict[parent]["positions"])
    mask1 = (index_np == parent) & (iteration_np >= p_len - k + w.iteration)
    
    sibling = walker_dict[walker_index]["sibling"]
    mask2 = (index_np == sibling) & (iteration_np <= k)
    
    grandparent = walker_dict[parent]["parent"]
    if (grandparent != None) and (p_len < k):
        gp_len = len(walker_dict[grandparent]["positions"])
        mask3 = (index_np == grandparent) & (iteration_np >= gp_len - k + p_len + w.iteration)
        
        uncle = walker_dict[parent]["sibling"]
        mask4 = (index_np == uncle) & (iteration_np <= k)
        
        return combine_masks([mask1, mask2, mask3, mask4])
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
    
"""
6: A few mathematical formulas for equation plotting
"""

#integral calculator for functions such as Bessel's functions
def integral(f, a, b, n = 100):
    width = (b-a)/n
    integral = 0
    
    for i in range(0,n+1):
        x = a + i * width
        prefactor = 2
        if (i == 0 or i ==n):
            prefactor = 1
        if (i -1)/2 - int((i -1)/2) == 0:
            prefactor = 4
        integral += prefactor * f(x)
    
    return width/3 * integral


def i0(v, n = 100):
    def integrand(theta):
        return maths.exp(v * maths.cos(theta))
    return 1/maths.pi * integral(integrand, 0, maths.pi, n)


def vonmis(theta, v):
    normalisation = 2 * maths.pi * i0(v)
    return normalisation**(-1) * np.exp(v * np.cos(theta))



"""
import scipy.special as sp
#Evaluate the modified Bessel function of the first kind, zeroth order, at val
val = 3
value = sp.iv(0, val)
print("SciPy value: ", value)
"""
